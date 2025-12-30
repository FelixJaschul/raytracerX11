#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define XKEYS_IMPLEMENTATION
#define XMATH_IMPLEMENTATION
#define XUTIL_IMPLEMENTATION

#include <wrapperX11/x11.h>

#define XTITLE "X11"
#define XFPS 60

#define SENS_X 0.40f
#define SENS_Y 0.15f

#define CLAMP(x, low, high) (fmaxf((low), fminf((x), (high))))

#define EPSILON          0.001f
#define AMBIENT_STRENGTH 0.2f
#define DIFFUSE_STRENGTH 0.6f
#define MAX_BOUNCES      2

typedef struct {
    bool hit;
    float t;
    Vec3 point;
    Vec3 normal;
    xMaterial mat;
} HitRecord;

typedef struct {
    Vec3 v0, v1, v2;
} Triangle;

typedef struct {
    Triangle* triangles;
    Triangle* transformed_triangles;
    int num_triangles;
    Vec3 position;
    float rot_x, rot_y, rot_z;
    float scale;
    xMaterial mat;
} Model;

#define MAX_MODELS 100
Model scene_models[MAX_MODELS];
int num_models = 0;

// Rotation helpers
Vec3 rotate_x(const Vec3 v, const float a)
{
    const float c = cosf(a), s = sinf(a);
    return vec3(v.x, v.y * c - v.z * s, v.y * s + v.z * c);
}

Vec3 rotate_y(const Vec3 v, const float a)
{
    const float c = cosf(a), s = sinf(a);
    return vec3(v.x * c - v.z * s, v.y, v.x * s + v.z * c);
}

Vec3 rotate_z(const Vec3 v, const float a)
{
    const float c = cosf(a), s = sinf(a);
    return vec3(v.x * c - v.y * s, v.x * s + v.y * c, v.z);
}

Vec3 transform_vertex(Vec3 v, const Model* m)
{
    v = mul(v, m->scale);
    if (m->rot_z) v = rotate_z(v, m->rot_z);
    if (m->rot_x) v = rotate_x(v, m->rot_x);
    if (m->rot_y) v = rotate_y(v, m->rot_y);
    return add(v, m->position);
}

Model* create_model(const Vec3 color, const float refl)
{
    if (num_models >= MAX_MODELS) return NULL;
    
    Model* m = &scene_models[num_models++];
    *m = (Model){
        .position = {0, 0, 0},
        .scale = 1.0f,
        .mat = {color, refl, 0.0f}
    };
    return m;
}

void load(Model* m, const char* path)
{
    FILE* f = fopen(path, "r");
    if (!f) {
        printf("Failed to load: %s\n", path);
        return;
    }
    
    Vec3 verts[1000];
    Triangle tris[1000];
    int nv = 0, nt = 0;
    
    char buf[256];
    while (fgets(buf, sizeof(buf), f)) 
    {
        if (buf[0] == 'v' && buf[1] == ' ') 
        {
            float x, y, z;
            sscanf(buf + 2, "%f %f %f", &x, &y, &z);
            verts[nv++] = vec3(x, y, z);
        }
        else if (buf[0] == 'f') {
            int a, b, c;
            sscanf(buf + 2, "%d %d %d", &a, &b, &c);
            tris[nt++] = (Triangle){verts[a-1], verts[b-1], verts[c-1]};
        }
    }
    fclose(f);
    
    m->triangles = malloc(nt * sizeof(Triangle));
    m->transformed_triangles = malloc(nt * sizeof(Triangle));
    memcpy(m->triangles, tris, nt * sizeof(Triangle));
    m->num_triangles = nt;
    printf("Loaded %s: %d triangles\n", path, nt);
}

void transform(Model* m, const Vec3 pos, const Vec3 rot, const float scale)
{
    m->position = pos;
    m->rot_x = rot.x;
    m->rot_y = rot.y;
    m->rot_z = rot.z;
    m->scale = scale;
}

void update_transformed_vertices()
{
    for (int i = 0; i < num_models; i++)
    {
        const Model* m = &scene_models[i];
        for (int j = 0; j < m->num_triangles; j++)
        {
            m->transformed_triangles[j].v0 = transform_vertex(m->triangles[j].v0, m);
            m->transformed_triangles[j].v1 = transform_vertex(m->triangles[j].v1, m);
            m->transformed_triangles[j].v2 = transform_vertex(m->triangles[j].v2, m);
        }
    }
}

// Scene setup
void scene_init() 
{
    Model* cube1 = create_model(vec3(1.0f, 0.5f, 0.3f), 0.3f);
    load(cube1, "cube.obj");
    transform(cube1, vec3(0.0f, 1.0f, 0.0f), vec3(0, 0, 0), 1.0f);
    
    Model* cube2 = create_model(vec3(0.3f, 0.5f, 1.0f), 0.5f);
    load(cube2, "cube.obj");
    transform(cube2, vec3(2.0f, 2.0f, 0.0f), vec3(0, M_PI/4, 0), 0.5f);
    
    Model* cube3 = create_model(vec3(0.3f, 1.0f, 0.5f), 0.2f);
    load(cube3, "cube.obj");
    transform(cube3, vec3(-2.0f, 1.5f, -1.0f), vec3(M_PI/6, M_PI/3, 0), 0.8f);

    update_transformed_vertices();
}

// Ray-triangle intersection
bool intersect_triangle(const Ray ray, const Triangle tri, const xMaterial mat, HitRecord *rec)
{
    const Vec3 v0 = tri.v0;
    const Vec3 v1 = tri.v1;
    const Vec3 v2 = tri.v2;
    
    const Vec3 edge1 = sub(v1, v0);
    const Vec3 edge2 = sub(v2, v0);
    const Vec3 h  = cross(ray.direction, edge2);
    const float a = dot(edge1, h);
    
    if (fabsf(a) < EPSILON) return false;
    
    const float f = 1.0f / a;
    const Vec3 s = sub(ray.origin, v0);
    const float u = f * dot(s, h);
    
    if (u < 0.0f || u > 1.0f) return false;
    
    const Vec3 q = cross(s, edge1);
    const float v = f * dot(ray.direction, q);
    
    if (v < 0.0f || u + v > 1.0f) return false;
    
    const float t = f * dot(edge2, q);
    
    if (t < EPSILON || t >= rec->t) return false;
    
    rec->hit = true;
    rec->t = t;
    rec->point = add(ray.origin, mul(ray.direction, t));
    rec->normal = norm(cross(edge1, edge2));
    rec->mat = mat;
    return true;
}

bool trace_scene(const Ray ray, HitRecord *closest_rec)
{
    closest_rec->hit = false;
    closest_rec->t = FLT_MAX;
    
    for (int i = 0; i < num_models; i++)
    {
        const Model* model = &scene_models[i];
        for (int j = 0; j < model->num_triangles; j++)
        {
            intersect_triangle(ray, model->transformed_triangles[j], model->mat, closest_rec);
        }
    }
    
    return closest_rec->hit;
}

Vec3 compute_lighting(const Vec3 point, const Vec3 normal, const Vec3 view_dir, const xMaterial mat)
{
    const Vec3 light_pos = vec3(0.0f, 5.0f, 3.0f);
    const Vec3 light_vec = sub(light_pos, point);
    const float light_len_sq = dot(light_vec, light_vec);
    const float light_len = sqrtf(light_len_sq);
    const Vec3 light_dir = mul(light_vec, 1.0f / light_len);

    const Vec3 ambient = mul(mat.color, AMBIENT_STRENGTH);
    const float diff = fmaxf(dot(normal, light_dir), 0.0f);
    const Vec3 diffuse = mul(mat.color, diff * DIFFUSE_STRENGTH);

    const Vec3 reflect_dir = reflect(mul(light_dir, -1.0f), normal);
    float spec = fmaxf(dot(view_dir, reflect_dir), 0.0f);
    spec *= spec; spec *= spec; spec *= spec; spec *= spec; spec *= spec;

    const Vec3 specular = vec3(spec * mat.specular, spec * mat.specular, spec * mat.specular);

    return add(add(ambient, diffuse), specular);
}

Vec3 ray_color(const Ray ray, const int depth)
{
    HitRecord rec;
    if (trace_scene(ray, &rec))
    {
        const Vec3 view_dir = mul(ray.direction, -1.0f);
        Vec3 color = compute_lighting(rec.point, rec.normal, view_dir, rec.mat);

        if (depth > 1 && rec.mat.reflectivity > 0.0f)
        {
            const Vec3 reflect_dir = reflect(ray.direction, rec.normal);
            const Ray reflect_ray = {rec.point, reflect_dir};
            const Vec3 reflect_color = ray_color(reflect_ray, depth - 1);
            color = add(mul(color, 1.0f - rec.mat.reflectivity), mul(reflect_color, rec.mat.reflectivity));
        }
        return color;
    }
    return vec3(0.0f, 0.0f, 0.0f);
}

static inline uint32_t color_to_uint32(const Vec3 color)
{
    const float r = CLAMP(color.x, 0.0f, 1.0f);
    const float g = CLAMP(color.y, 0.0f, 1.0f);
    const float b = CLAMP(color.z, 0.0f, 1.0f);
    return ((int)(r * 255) << 16) | ((int)(g * 255) << 8) | (int)(b * 255);
}

int main()
{
    xWindow win;
    xWindowInit(&win);
    win.title = XTITLE;
    win.fps = XFPS;
    xCreateWindow(&win);

    xCamera camera;
    xCameraInit(&camera);
    camera.position = vec3(2.0f, 2.0f, 0.6f);
    
    scene_init();

    const float viewport_height = 2.0f;
    const float viewport_width = (float)win.width / (float)win.height * viewport_height;

    // Precompute viewport offsets
    float* u_offsets = malloc(win.width * sizeof(float));
    float* v_offsets = malloc(win.height * sizeof(float));
    for (int x = 0; x < win.width; x++)  u_offsets[x] = ((float)x / (float)(win.width - 1) - 0.5f) * viewport_width;
    for (int y = 0; y < win.height; y++) v_offsets[y] = ((float)(win.height - 1 - y) / (float)(win.height - 1) - 0.5f) * viewport_height;

    bool mouse_grabbed = false;

    while (1)
    {
        const float move_speed = 0.03f;
        
        if (xPollEvents(win.display)) break;
        if (xIsKeyPressed(KEY_ESCAPE)) break;

        // Toggle mouse grab with mouse click
        if (xIsMousePressed(MOUSE_LEFT) && !mouse_grabbed)
        {
            xGrabMouse(win.display, win.window, win.width, win.height);
            mouse_grabbed = true;
        }

        // Mouse camera control
        if (mouse_grabbed)
        {
            int dx, dy;
            xGetMouseDelta(&dx, &dy);
            xCameraRotate(&camera, dx * SENS_X, -dy * SENS_Y);
        }

        // Keyboard movement
        if (xIsKeyDown(KEY_W)) xCameraMove(&camera, camera.front, move_speed);
        if (xIsKeyDown(KEY_S)) xCameraMove(&camera, mul(camera.front, -1), move_speed);
        if (xIsKeyDown(KEY_A)) xCameraMove(&camera, mul(camera.right, -1), move_speed);
        if (xIsKeyDown(KEY_D)) xCameraMove(&camera, camera.right, move_speed);
        if (xIsKeyDown(KEY_Q)) xCameraMove(&camera, vec3(0, -1, 0), move_speed);
        if (xIsKeyDown(KEY_E)) xCameraMove(&camera, vec3(0, 1, 0), move_speed);

        // Render
        #pragma omp parallel for schedule(dynamic) default(none) shared(win, camera, u_offsets, v_offsets)
        for (int y = 0; y < win.height; y++)
        {
            const Vec3 row_offset = mul(camera.up, v_offsets[y]);
            for (int x = 0; x < win.width; x++)
            {
                Vec3 rd = add(camera.front, add(row_offset, mul(camera.right, u_offsets[x])));
                rd = norm(rd);

                const Ray ray = {camera.position, rd};
                const Vec3 color = ray_color(ray, MAX_BOUNCES);
                win.buffer[y * win.width + x] = color_to_uint32(color);
            }
        }

        xUpdateFramebuffer(&win);
        xUpdateFrame(&win);
        xUpdateInput();
    }

    free(u_offsets);
    free(v_offsets);
    
    // Cleanup models
    for (int i = 0; i < num_models; i++) {
        if (scene_models[i].triangles) free(scene_models[i].triangles);
        if (scene_models[i].transformed_triangles) free(scene_models[i].transformed_triangles);
    }

    xDestroyWindow(&win);
    return 0;
}
