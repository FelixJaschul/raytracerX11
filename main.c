#include <sys/types.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <x11.h>
#include <xMath.h>

#define XTITLE "X11"
#define XFPS 60

#define CLAMP(x, low, high) (fmaxf((low), fminf((x), (high))))

#define SPECULAR_VALUE   0.0f
#define EPSILON          0.001f
#define AMBIENT_STRENGTH 0.2f
#define DIFFUSE_STRENGTH 0.6f

typedef struct 
{
    Vec3 color;
    float reflectivity;
    float specular;
} 
Material;

typedef struct 
{
    Vec3 v0, v1, v2;
} 
Triangle;

typedef struct
{
    Triangle* triangles;
    int num_triangles;
    Vec3 position;
    float rot_x, rot_y, rot_z;
    float scale;
    Material mat;
}
Model;

typedef struct 
{
    bool hit;
    float t;
    Vec3 point;
    Vec3 normal;
    Material mat;
} 
HitRecord;

#define MAX_MODELS 100
#define MAX_TRIANGLES 1000
#define MAX_BOUNCES 2

Model scene_models[MAX_MODELS];
int num_models = 0;

// Rotation matrices
Vec3 rotate_x(Vec3 v, float a) 
{
    float c = cosf(a), s = sinf(a);
    return vec3(v.x, v.y * c - v.z * s, v.y * s + v.z * c);
}

Vec3 rotate_y(Vec3 v, float a) 
{
    float c = cosf(a), s = sinf(a);
    return vec3(v.x * c - v.z * s, v.y, v.x * s + v.z * c);
}

Vec3 rotate_z(Vec3 v, float a) 
{
    float c = cosf(a), s = sinf(a);
    return vec3(v.x * c - v.y * s, v.x * s + v.y * c, v.z);
}

// Transform vertex with model matrix
Vec3 transform_vertex(Vec3 v, Model* m) 
{
    v = mul(v, m->scale);
    if (m->rot_z) v = rotate_z(v, m->rot_z);
    if (m->rot_x) v = rotate_x(v, m->rot_x);
    if (m->rot_y) v = rotate_y(v, m->rot_y);
    return add(v, m->position);
}

Model* create(Vec3 color, float refl) 
{
    if (num_models >= MAX_MODELS) return NULL;
    
    Model* m = &scene_models[num_models++];
    *m = (Model){
        .position = {0, 0, 0},
        .scale = 1.0f,
        .mat = {color, refl, SPECULAR_VALUE}};
    return m;
}

void load(Model* m, const char* path) 
{
    FILE* f = fopen(path, "r");
    if (!f) return;
    
    Vec3 verts[MAX_TRIANGLES];
    Triangle tris[MAX_TRIANGLES];
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
    memcpy(m->triangles, tris, nt * sizeof(Triangle));
    m->num_triangles = nt;
}

void transform(Model* m, Vec3 pos, Vec3 rot, float scale) 
{
    m->position = pos;
    m->rot_x = rot.x;
    m->rot_y = rot.y;
    m->rot_z = rot.z;
    m->scale = scale;
}

void sInit() 
{
    Model* cube1 = create(vec3(1.0f, 0.5f, 0.3f), 0.3f);
    load(cube1, "cube.obj");
    transform(cube1, vec3(0, 1, 0), vec3(0, 0, 0), 1.0f);
    
    Model* cube2 = create(vec3(0.3f, 0.5f, 1.0f), 0.5f);
    load(cube2, "cube.obj");
    transform(cube2, vec3(2, 2, 0), vec3(0, M_PI/4, 0), 0.5f);
}

// Forward declarations
Vec3 ray_color(Ray ray, int depth);
Vec3 compute_lighting(Vec3 point, Vec3 normal, Vec3 view_dir, Material mat);
bool trace_scene(Ray ray, HitRecord *closest_rec);
static inline uint32_t color_to_uint32(Vec3 color);

int main()
{
    xWindow win;
    xInit(&win);
    win.title = XTITLE;
    win.fps = XFPS;
    xCreateWindow(&win);

    Vec3 camera_pos = vec3(2, 2, 2);
    float yaw = -135.0f;
    float pit =    0.0f;
    sInit();

    const float viewport_height = 2.0f;
    const float viewport_width = (float)win.width / (float)win.height * viewport_height;

    float* u_offsets = malloc(win.width * sizeof(float));
    float* v_offsets = malloc(win.height * sizeof(float));
    for (int x = 0; x < win.width; x++) u_offsets[x] = ((float)x / (float)(win.width - 1) - 0.5f) * viewport_width;
    for (int y = 0; y < win.height; y++) v_offsets[y] = ((float)(win.height - 1 - y) / (float)(win.height - 1) - 0.5f) * viewport_height;

    while (1)
    {
        const float rspeed = 1.9f;
        const float mspeed = 0.03f;
        if (xPollEvents(win.display)) break;
        if (xIsKeyPressed(Escape))    break;

        if (xIsKeyDown(Left))  yaw -= rspeed;
        if (xIsKeyDown(Right)) yaw += rspeed;
        if (xIsKeyDown(Up))    pit += rspeed;
        if (xIsKeyDown(Down))  pit -= rspeed;

        if (pit >  89.0f) pit =  89.0f;
        if (pit < -89.0f) pit = -89.0f;

        const float sin_pit = sinf(pit * M_PI / 180.0f);
        const float cos_pit = cosf(pit * M_PI / 180.0f);
        const float sin_yaw = sinf(yaw * M_PI / 180.0f);
        const float cos_yaw = cosf(yaw * M_PI / 180.0f);

        Vec3 front = {cos_yaw * cos_pit, sin_pit, sin_yaw * cos_pit};
        front = norm(front);

        const Vec3 right = norm(cross(front, vec3(0, 1, 0)));
        const Vec3 up    = cross(right, front);

        if (xIsKeyDown(W)) camera_pos = add(camera_pos, mul(front, mspeed));
        if (xIsKeyDown(S)) camera_pos = sub(camera_pos, mul(front, mspeed));
        if (xIsKeyDown(A)) camera_pos = sub(camera_pos, mul(right, mspeed));
        if (xIsKeyDown(D)) camera_pos = add(camera_pos, mul(right, mspeed));

        #pragma omp parallel for schedule(dynamic)
        for (int y = 0; y < win.height; y++)
        {
            const Vec3 row_offset = mul(up, v_offsets[y]);
            for (int x = 0; x < win.width; x++)
            {
                Vec3 rd = add(front, add(row_offset, mul(right, u_offsets[x])));
                rd = norm(rd);

                const Ray ray = {camera_pos, rd};
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
    xDestroyWindow(&win);
    return 0;
}

static inline uint32_t color_to_uint32(Vec3 color)
{
    const float r = CLAMP(color.x, 0.0f, 1.0f);
    const float g = CLAMP(color.y, 0.0f, 1.0f);
    const float b = CLAMP(color.z, 0.0f, 1.0f);
    return ((int)(r * 255) << 16) | ((int)(g * 255) << 8) | (int)(b * 255);
}

bool intersect_triangle(Ray ray, Triangle tri, Model* model, HitRecord *rec)
{
    // Transform triangle vertices
    Vec3 v0 = transform_vertex(tri.v0, model);
    Vec3 v1 = transform_vertex(tri.v1, model);
    Vec3 v2 = transform_vertex(tri.v2, model);
    
    const Vec3 edge1 = sub(v1, v0);
    const Vec3 edge2 = sub(v2, v0);
    const Vec3 h = cross(ray.direction, edge2);
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
    rec->mat = model->mat;
    return true;
}

bool trace_scene(Ray ray, HitRecord *closest_rec)
{
    closest_rec->hit = false;
    closest_rec->t = FLT_MAX;
    
    for (int i = 0; i < num_models; i++)
    {
        Model* model = &scene_models[i];
        for (int j = 0; j < model->num_triangles; j++)
        {
            intersect_triangle(ray, model->triangles[j], model, closest_rec);
        }
    }
    
    return closest_rec->hit;
}

Vec3 compute_lighting(Vec3 point, Vec3 normal, Vec3 view_dir, Material mat)
{
    const Vec3 light_pos = vec3(0.0f, 4.0f, 0.0f);
    const Vec3 light_vec = sub(light_pos, point);
    const float light_dist = len(light_vec);
    const Vec3 light_dir = mul(light_vec, 1.0f / light_dist);

    const Vec3 ambient = mul(mat.color, AMBIENT_STRENGTH);

    const float dot_nl = dot(normal, light_dir);
    const float diff = fmaxf(dot_nl, 0.0f);
    const Vec3 diffuse = mul(mat.color, diff * DIFFUSE_STRENGTH);

    const Vec3 reflect_dir = reflect(mul(light_dir, -1.0f), normal);
    float dot_vr = fmaxf(dot(view_dir, reflect_dir), 0.0f);

    float spec = dot_vr;
    spec *= spec;
    spec *= spec;
    spec *= spec;
    spec *= spec;
    spec *= spec;

    const Vec3 specular = vec3(spec * mat.specular, spec * mat.specular, spec * mat.specular);

    return add(add(ambient, diffuse), specular);
}

Vec3 ray_color(Ray ray, int depth)
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
    return vec3(0, 0, 0);
}

