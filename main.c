#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <stdint.h>

#define XKEYS_IMPLEMENTATION
#define XMATH_IMPLEMENTATION
#define XCAMERA_IMPLEMENTATION
#define XMODEL_IMPLEMENTATION
#include <wrapperX11/x11.h>

#define XTITLE "X11"
#define XFPS 120
#define SENS_X 0.40f
#define SENS_Y 0.15f

#define CLAMP(x, low, high) (fmaxf((low), fminf((x), (high))))

#define EPSILON          0.001f
#define RAY_OFFSET       0.0002f
#define AMBIENT_STRENGTH 0.2f
#define DIFFUSE_STRENGTH 0.6f

#define MAX_BOUNCES 2
#define MAX_MODELS  10
#define TILE_SIZE   16

typedef struct {
    bool hit;
    float t;
    Vec3 point;
    Vec3 normal;
    xMaterial mat;
} HitRecord;

#define UTIL_IMPLEMENTATION
#include "util.h"

xModel scene_models[MAX_MODELS];
int num_models = 0;
BVHNode *bvh_root = NULL;

#define TOGGLE_REFLECTIVITY 1
#define WALL_REFLECTIVITY (0.3f * TOGGLE_REFLECTIVITY)

xModel* model(const char* path, xModel* storage, int* count, const Vec3 color, const float refl)
{
    xModel *f = xModelCreate(storage, count, MAX_MODELS, color, WALL_REFLECTIVITY * refl);
    xModelLoad(f, path);
    return f;
}

// Scene setup
void scene_init()
{
    model("res/rect.obj",  scene_models, &num_models, vec3(0.0f, 0.0f, 1.0f), 1.0f);
    model("res/rect.obj",  scene_models, &num_models, vec3(0.0f, 1.0f, 0.0f), 1.0f);
    model("res/rect.obj",  scene_models, &num_models, vec3(1.0f, 0.0f, 0.0f), 1.0f);
    model("res/bunni.obj", scene_models, &num_models, vec3(1.0f, 1.0f, 1.0f), 0.3f);
    xModelTransform(&scene_models[0], vec3(-2.0f,  0.0f,  2.0f), vec3(-M_PI/2, 0, 0), vec3(4.0f, 4.0f, 1.0f));
    xModelTransform(&scene_models[1], vec3(-2.0f,  0.0f,  2.0f), vec3(0, -M_PI/2, 0), vec3(4.0f, 4.0f, 1.0f));
    xModelTransform(&scene_models[2], vec3(-2.0f,  0.0f, -2.0f), vec3(0,       0, 0), vec3(4.0f, 4.0f, 1.0f));

    xModelUpdate(scene_models, num_models);

    // Build BVH after all models are loaded and transformed
    bvh_build(&bvh_root, scene_models, num_models);
    printf("BVH built: %s\n", bvh_root ? "YES" : "NO");
}

#undef TOGGLE_REFLECTIVITY
#undef WALL_REFLECTIVITY

bool trace_scene(const Ray ray, HitRecord *rec)
{
    rec->hit = false;
    rec->t   = FLT_MAX;
    return bvh_intersect(bvh_root, ray, rec);
}

Vec3 calculate_lighting(const Vec3 point, const Vec3 normal, const Vec3 view_dir, const xMaterial mat)
{
    const Vec3 light_pos = vec3(0.0f, 3.8f, 0.0f);
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

        // Add Fresnel-like effect (edges more reflective)
        const float facing = dot(view_dir, rec.normal);
        const float adjusted_reflectivity = rec.mat.reflectivity * (1.0f - facing * 0.5f);

        if (depth > 1 && adjusted_reflectivity > 0.05f)
        {
            const Vec3 reflect_dir = reflect(ray.direction, rec.normal);
            const Ray reflect_ray = {rec.point, reflect_dir};
            const Vec3 reflect_color = calculate_ray_color(reflect_ray, depth - 1);
            color = add(mul(color, 1.0f - adjusted_reflectivity), mul(reflect_color, adjusted_reflectivity));
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

    // Initialize input
    xInput input;
    xInputInit(&input);

    // Initialize camera
    xCamera camera;
    xCameraInit(&camera);
    camera.position = vec3(0.0f, 0.7f, 2.0f);
    camera.yaw = -90.0f;

    // Load scene
    scene_init();

    // Precompute viewport offsets
    const float viewport_height = 2.0f;
    const float viewport_width = (float)win.width / (float)win.height * viewport_height;

    float* u_offsets = malloc(win.width * sizeof(float));
    float* v_offsets = malloc(win.height * sizeof(float));
    assert(u_offsets && v_offsets && "Failed to allocate viewport offset buffers");

    for (int x = 0; x < win.width; x++)  u_offsets[x] = ((float)x / (float)(win.width - 1) - 0.5f) * viewport_width;
    for (int y = 0; y < win.height; y++) v_offsets[y] = ((float)(win.height - 1 - y) / (float)(win.height - 1) - 0.5f) * viewport_height;

    int x = 0;
    // Main loop
    while (1)
    {
        if (x == 360) x = 0;
        x++;
        printf("FPS: %.2f\n", xGetFPS(&win));
        const float move_speed = 0.03f;

        // Update bunni rotation
        xModelTransform(&scene_models[3], vec3(0.0f, -0.33f, 0.0f), vec3(0, x * 0.01f, 0), vec3(10.0f, 10.0f, 10.0f));
        xModelUpdate(scene_models, num_models); bvh_free(bvh_root); // Free previous BVH
        bvh_build(&bvh_root, scene_models, num_models);

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
        #pragma omp parallel for schedule(dynamic) collapse(2) default(none) shared(win, camera, u_offsets, v_offsets)
        for (int ty = 0; ty < win.height; ty += TILE_SIZE)
        {
            for (int tx = 0; tx < win.width; tx += TILE_SIZE)
            {
                // Process tile
                const int y_end = (ty + TILE_SIZE < win.height) ? ty + TILE_SIZE : win.height;
                const int x_end = (tx + TILE_SIZE < win.width ) ? tx + TILE_SIZE : win.width;

                for (int y = ty; y < y_end; y++)
                {
                    uint32_t* restrict row = &win.buffer[y * win.width];
                    for (int x = tx; x < x_end; x++)
                    {
                        const Ray ray = xCameraGetRay(&camera, u_offsets[x], v_offsets[y]);
                        row[x] = uint32(calculate_ray_color(ray, MAX_BOUNCES));
                    }
                }
            }
        }


        xUpdateFramebuffer(&win);
        xUpdateFrame(&win);
        xUpdateInput(&input);
    }

    // Cleanup
    bvh_free(bvh_root);

    free(u_offsets);
    free(v_offsets);

    for (int i = 0; i < num_models; i++) xModelFree(&scene_models[i]);

    xDestroyWindow(&win);
    return 0;
}