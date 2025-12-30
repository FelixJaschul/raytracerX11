#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <omp.h>

#define XKEYS_IMPLEMENTATION
#define XMATH_IMPLEMENTATION
#define XCAMERA_IMPLEMENTATION
#define XMODEL_IMPLEMENTATION
#include <assert.h>

#include "wrapperX11/x11.h"

#define XTITLE "Raytracer X11"
#define XFPS 60
#define SENS_X 0.40f
#define SENS_Y 0.15f

#define CLAMP(x, low, high) (fmaxf((low), fminf((x), (high))))
#define EPSILON 0.001f
#define AMBIENT_STRENGTH 0.2f
#define DIFFUSE_STRENGTH 0.6f
#define MAX_BOUNCES 3
#define MAX_MODELS 100

typedef struct {
    bool hit;
    float t;
    Vec3 point;
    Vec3 normal;
    xMaterial mat;
} HitRecord;

xModel scene_models[MAX_MODELS];
int num_models = 0;

// Scene setup
void scene_init()
{
    #define TOGGLE_REFLECTIVITY 1
    xModel *cube_red = xModelCreate(scene_models, &num_models, MAX_MODELS, vec3(1.0f, 0.0f, 0.0f), 0.3f * TOGGLE_REFLECTIVITY);
    xModelLoad(cube_red, "res/cube.obj");
    xModelTransform(cube_red, vec3(0.0f, 1.0f, 0.0f), vec3(0, 0, 0), vec3(1.0f, 1.0f, 1.0f));

    xModel *cube_blue = xModelCreate(scene_models, &num_models, MAX_MODELS, vec3(0.0f, 0.0f, 1.0f), 0.5f  * TOGGLE_REFLECTIVITY);
    xModelLoad(cube_blue, "res/cube.obj");
    xModelTransform(cube_blue, vec3(2.0f, 2.0f, 0.0f), vec3(0, M_PI/4, 0), vec3(0.5f, 0.5f, 0.5f));

    xModel *cube_green = xModelCreate(scene_models, &num_models, MAX_MODELS, vec3(0.0f, 1.0f, 0.0f), 0.2f  * TOGGLE_REFLECTIVITY);
    xModelLoad(cube_green, "res/cube.obj");
    xModelTransform(cube_green, vec3(-2.0f, 1.5f, -1.0f), vec3(M_PI/6, M_PI/3, 0), vec3(0.8f, 0.8f, 0.8f));

    #define WALL_REFLECTIVITY (0.0f * TOGGLE_REFLECTIVITY)
    xModel *floor = xModelCreate(scene_models, &num_models, MAX_MODELS, vec3(1.0f, 1.0f, 1.0f), WALL_REFLECTIVITY);
    xModelLoad(floor, "res/rect.obj");
    xModelTransform(floor, vec3(-5.0f, 0.0f, 5.0f), vec3(-M_PI/2, 0, 0), vec3(10.0f, 10.0f, 100.0f));

    xModel *ceil = xModelCreate(scene_models, &num_models, MAX_MODELS, vec3(1.0f, 1.0f, 1.0f), WALL_REFLECTIVITY);
    xModelLoad(ceil, "res/rect.obj");
    xModelTransform(ceil, vec3(-5.0f, 3.0f, 5.0f), vec3(-M_PI/2, 0, 0), vec3(10.0f, 10.0f, 100.0f));

    xModel *left_wall = xModelCreate(scene_models, &num_models, MAX_MODELS, vec3(1.0f, 1.0f, 1.0f), WALL_REFLECTIVITY);
    xModelLoad(left_wall, "res/rect.obj");
    xModelTransform(left_wall, vec3(-5.0f, 0.0f, -5.0f), vec3(0, M_PI/2, 0), vec3(10.0f, 3.0f, 1.0f));

    xModel *right_wall = xModelCreate(scene_models, &num_models, MAX_MODELS, vec3(1.0f, 1.0f, 1.0f), WALL_REFLECTIVITY);
    xModelLoad(right_wall, "res/rect.obj");
    xModelTransform(right_wall, vec3(5.0f, 0.0f, -5.0f), vec3(0, M_PI/2, 0), vec3(10.0f, 3.0f, 1.0f));

    xModel *front_wall = xModelCreate(scene_models, &num_models, MAX_MODELS, vec3(1.0f, 1.0f, 1.0f), WALL_REFLECTIVITY);
    xModelLoad(front_wall, "res/rect.obj");
    xModelTransform(front_wall, vec3(-5.0f, 0.0f, -5.0f), vec3(0, 0, 0), vec3(10.0f, 3.0f, 1.0f));

    xModel *back_wall = xModelCreate(scene_models, &num_models, MAX_MODELS, vec3(1.0f, 1.0f, 1.0f), WALL_REFLECTIVITY);
    xModelLoad(back_wall, "res/rect.obj");
    xModelTransform(back_wall, vec3(-5.0f, 0.0f, 5.0f), vec3(0, 0, 0), vec3(10.0f, 3.0f, 1.0f));

    xModelUpdate(scene_models, num_models);
}

// Ray-triangle intersection
bool intersect_triangle(const Ray ray, const xTriangle tri, const xMaterial mat, HitRecord *rec)
{
    const Vec3 v0 = tri.v0;
    const Vec3 v1 = tri.v1;
    const Vec3 v2 = tri.v2;
    const Vec3 edge1 = sub(v1, v0);
    const Vec3 edge2 = sub(v2, v0);
    const Vec3 h = cross(ray.direction, edge2);
    const float a = dot(edge1, h);

    if (a < EPSILON) return false;

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
        const xModel* model = &scene_models[i];
        for (int j = 0; j < model->num_triangles; j++)
        {
            intersect_triangle(ray, model->transformed_triangles[j], model->mat, closest_rec);
        }
    }

    return closest_rec->hit;
}

Vec3 calculate_lighting(const Vec3 point, const Vec3 normal, const Vec3 view_dir, const xMaterial mat)
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

Vec3 calculate_ray_color(const Ray ray, const int depth)
{
    HitRecord rec;
    if (trace_scene(ray, &rec))
    {
        const Vec3 view_dir = mul(ray.direction, -1.0f);
        Vec3 color = calculate_lighting(rec.point, rec.normal, view_dir, rec.mat);

        if (depth > 1 && rec.mat.reflectivity > 0.0f)
        {
            const Vec3 reflect_dir = reflect(ray.direction, rec.normal);
            const Ray reflect_ray = {rec.point, reflect_dir};
            const Vec3 reflect_color = calculate_ray_color(reflect_ray, depth - 1);
            color = add(mul(color, 1.0f - rec.mat.reflectivity), mul(reflect_color, rec.mat.reflectivity));
        }

        return color;
    }

    return vec3(0.0f, 0.0f, 0.0f);
}

static inline uint32_t uint32(const Vec3 color)
{
    const float r = CLAMP(color.x, 0.0f, 1.0f);
    const float g = CLAMP(color.y, 0.0f, 1.0f);
    const float b = CLAMP(color.z, 0.0f, 1.0f);
    return ((int)(r * 255) << 16) | ((int)(g * 255) << 8) | (int)(b * 255);
}

int main()
{
    // Initialize window
    xWindow win;
    xWindowInit(&win);
    win.title = XTITLE;
    win.fps = XFPS;

    if (!xCreateWindow(&win)) return 1;

    // Initialize input state
    xInputState input;
    xInputInit(&input);

    // Initialize camera
    xCamera camera;
    xCameraInit(&camera);
    camera.position = vec3(2.0f, 2.0f, 0.6f);
    camera.yaw = -90.0f;

    // Load scene
    scene_init();

    // Precompute viewport offsets (once at startup, not per frame!)
    const float viewport_height = 2.0f;
    const float viewport_width = (float)win.width / (float)win.height * viewport_height;

    float* u_offsets = (float*)malloc(win.width * sizeof(float));
    float* v_offsets = (float*)malloc(win.height * sizeof(float));
    assert(u_offsets && v_offsets && "Failed to allocate viewport offset buffers");

    for (int x = 0; x < win.width; x++)  u_offsets[x] = ((float)x / (float)(win.width - 1) - 0.5f) * viewport_width;
    for (int y = 0; y < win.height; y++) v_offsets[y] = ((float)(win.height - 1 - y) / (float)(win.height - 1) - 0.5f) * viewport_height;

    // Main loop
    while (1)
    {
        const float move_speed = 0.03f;

        // Poll events with explicit input state
        if (xPollEvents(win.display, &input)) break;
        if (xIsKeyPressed(&input, KEY_ESCAPE)) break;

        // Mouse look
        int dx, dy;
        xGetMouseDelta(&input, &dx, &dy);
        xGrabMouse(win.display, win.window, win.width, win.height, &input);
        xCameraRotate(&camera, dx * SENS_X, -dy * SENS_Y);

        // Keyboard movement
        if (xIsKeyDown(&input, KEY_W)) xCameraMove(&camera, camera.front, move_speed);
        if (xIsKeyDown(&input, KEY_S)) xCameraMove(&camera, mul(camera.front, -1), move_speed);
        if (xIsKeyDown(&input, KEY_A)) xCameraMove(&camera, mul(camera.right, -1), move_speed);
        if (xIsKeyDown(&input, KEY_D)) xCameraMove(&camera, camera.right, move_speed);
        if (xIsKeyDown(&input, KEY_Q)) xCameraMove(&camera, vec3(0, -1, 0), move_speed);
        if (xIsKeyDown(&input, KEY_E)) xCameraMove(&camera, vec3(0, 1, 0), move_speed);

        // Render
        #pragma omp parallel for schedule(dynamic) default(none) shared(win, camera, u_offsets, v_offsets)
        for (int y = 0; y < win.height; y++)
        {
            for (int x = 0; x < win.width; x++)
            {
                const Ray ray = xCameraGetRay(&camera, u_offsets[x], v_offsets[y]);
                const Vec3 color = calculate_ray_color(ray, MAX_BOUNCES);
                win.buffer[y * win.width + x] = uint32(color);
            }
        }


        xUpdateFramebuffer(&win);
        xUpdateFrame(&win);
        xUpdateInput(&input);
    }

    // Cleanup
    free(u_offsets);
    free(v_offsets);

    for (int i = 0; i < num_models; i++) xModelFree(&scene_models[i]);

    xDestroyWindow(&win);
    return 0;
}