#include <sys/types.h>
#include <stdbool.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
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
    Vec3 center;
    float radius;
    Material mat;
} 
Sphere;

typedef struct 
{
    Vec3 point;
    Vec3 normal;
    Vec3 u;
    Vec3 v;
    float width;
    float height;
    Material mat;
} 
Rect;

typedef struct 
{
    bool hit;
    float t;
    Vec3 point;
    Vec3 normal;
    Material mat;
} 
HitRecord;

#define MAX_SPHERES 1
#define MAX_RECTS 12
#define MAX_BOUNCES 3

Sphere scene_spheres[MAX_SPHERES];
Rect   scene_rects[MAX_RECTS];

int num_spheres = 0;
int num_rects   = 0;

void add_sphere(
    Vec3 center,
    float radius,
    Vec3 color,
    float refl)
{
    if (num_spheres < MAX_SPHERES) 
        scene_spheres[num_spheres++] = (Sphere) {
                center, 
                radius, 
                {color, refl, SPECULAR_VALUE}
        };
}

void add_rect(
    Vec3 point,
    Vec3 normal,
    Vec3 u,
    Vec3 v,
    float width,
    float height,
    Vec3 color,
    float refl)
{
    if (num_rects < MAX_RECTS) 
        scene_rects[num_rects++] = (Rect) {
            point, 
            norm(normal), 
            norm(u), 
            norm(v), 
            width, 
            height, 
            {color, refl, SPECULAR_VALUE}
        };
}

void add_cube(
    Vec3 center,
    float x,
    float y,
    float z,
    Vec3 color,
    float refl)
{
    const float half_x = x / 2.0f;
    const float half_y = y / 2.0f;
    const float half_z = z / 2.0f;
    
    add_rect(vec3(center.x,          center.y - half_y, center.z),          vec3( 0, -1,  0), vec3(1, 0, 0), vec3(0, 0, 1), x, z, color, refl);
    add_rect(vec3(center.x,          center.y + half_y, center.z),          vec3( 0,  1,  0), vec3(1, 0, 0), vec3(0, 0, 1), x, z, color, refl);
    add_rect(vec3(center.x - half_x, center.y,          center.z),          vec3(-1,  0,  0), vec3(0, 0, 1), vec3(0, 1, 0), z, y, color, refl);
    add_rect(vec3(center.x + half_x, center.y,          center.z),          vec3( 1,  0,  0), vec3(0, 0, 1), vec3(0, 1, 0), z, y, color, refl);
    add_rect(vec3(center.x,          center.y,          center.z - half_z), vec3( 0,  0, -1), vec3(1, 0, 0), vec3(0, 1, 0), x, y, color, refl);
    add_rect(vec3(center.x,          center.y,          center.z + half_z), vec3( 0,  0,  1), vec3(1, 0, 0), vec3(0, 1, 0), x, y, color, refl);
}


void sInit()
{
    /*
    add_rect(
             Vec3 point,            -> Center point 
             Vec3 normal,           -> Direction the rectangle faces
             Vec3 u,                -> 1. edge direction (normalized)
             Vec3 v,                -> 2. edge direction (normalized, perpendicular to u)
             float width,           -> Size along u direction
             float height,          -> Size along v direction  
             Vec3 color,            -> Color RGB (0.0 to 1.0)
             float refl,            -> Reflectivity (0.0 = no reflection, 1.0 = mirror)
    )
    */

    add_rect(
             vec3( 0.0f,  0.0f,  0.0f), 
             vec3( 0.0f,  1.0f,  0.0f), 
             vec3( 1, 0, 0), 
             vec3( 0, 0, 1), 
             6.0f, 
             6.0f, 
             vec3( 0.2f,  0.2f,  1.0f), 
             0.30f 
    );

    add_rect(
             vec3( 0.0f,  4.0f,  0.0f), 
             vec3( 0.0f, -1.0f,  0.0f), 
             vec3( 1, 0, 0), 
             vec3( 0, 0, 1), 
             6.0f, 
             6.0f, 
             vec3( 1.0f,  1.0f,  1.0f), 
             0.08f
    );

    add_rect(
             vec3(-3.0f,  2.0f,  0.0f), 
             vec3( 1.0f,  0.0f,  0.0f), 
             vec3( 0, 0, 1), 
             vec3( 0, 1, 0), 
             6.0f, 
             6.0f, 
             vec3( 1.0f,  0.2f,  0.2f), 
             0.08f
    );

    add_rect(
             vec3( 3.0f,  2.0f,  0.0f), 
             vec3(-1.0f,  0.0f,  0.0f), 
             vec3( 0, 0, 1), 
             vec3( 0, 1, 0), 
             6.0f, 
             6.0f, 
             vec3( 0.2f,  1.0f,  0.2f), 
             0.08f
    );

    add_rect(
             vec3( 0.0f,  2.0f, -3.0f), 
             vec3( 0.0f,  0.0f,  1.0f), 
             vec3( 1, 0, 0), 
             vec3( 0, 1, 0), 
             6.0f, 
             6.0f, 
             vec3( 1.0f,  1.0f,  1.0f), 
             0.08f
    );

    add_rect(
             vec3( 0.0f,  2.0f, +3.0f), 
             vec3( 0.0f,  0.0f,  1.0f), 
             vec3( 1, 0, 0), 
             vec3( 0, 1, 0), 
             6.0f, 
             6.0f, 
             vec3( 0.5f,  0.0f,  0.8f), 
             0.08f
    );
    
    add_cube(
             vec3( 0.0f,  2.0f,  0.0f),
             1.0f,
             1.0f,
             1.0f,
             vec3( 1.0f,  1.0f,  1.0f),
             0.80f
    );
    
    /*
    add_sphere(
             Vec3 center,           -> Pos in 3d space
             float radius,          -> Size
             Vec3 color,            -> Color RGB (0.0 to 1.0)
             float refl,            -> Reflectivity (0.0 = no reflection, 1.0 = mirror)
    )
    */ /*
    
    add_sphere(
             vec3( 0.0f,  2.0f,  0.0f), 
             1.0f, 
             vec3( 1.0f,  1.0f,  1.0f), 
             0.8f
    ); */
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
    float pit = 0.0f;
    sInit();

    const float viewport_height = 2.0f;
    const float viewport_width = (float)win.width / (float)win.height * viewport_height;

    // Precompute u/v offsets
    float* u_offsets = malloc(win.width * sizeof(float));
    float* v_offsets = malloc(win.height * sizeof(float));
    for (int x = 0; x < win.width; x++) u_offsets[x] = ((float)x / (float)(win.width - 1) - 0.5f) * viewport_width;
    for (int y = 0; y < win.height; y++) v_offsets[y] = ((float)(win.height - 1 - y) / (float)(win.height - 1) - 0.5f) * viewport_height;

    while (1)
    {
        const float rspeed = 0.4f;
        const float mspeed = 0.02f;
        if (xPollEvents(win.display)) break;
        if (xIsKeyPressed(Escape))    break;

        if (xIsKeyDown(Left))  yaw -= rspeed;
        if (xIsKeyDown(Right)) yaw += rspeed;
        if (xIsKeyDown(Up))    pit += rspeed;
        if (xIsKeyDown(Down))  pit -= rspeed;

        if (pit >  89.0f) pit =  89.0f;
        if (pit < -89.0f) pit = -89.0f;

        // Camera front vector
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

// Forward declarations
bool intersect_sphere(Ray ray, Sphere sphere, HitRecord *rec);
bool intersect_rect(Ray ray, Rect rect, HitRecord *rec);

bool trace_scene(Ray ray, HitRecord *closest_rec)
{
    closest_rec->hit = false;
    closest_rec->t = FLT_MAX;
    for (int i = 0; i < num_spheres; i++) intersect_sphere(ray, scene_spheres[i], closest_rec);
    for (int i = 0; i < num_rects;   i++) intersect_rect(ray, scene_rects[i], closest_rec);
    return closest_rec->hit;
}

Vec3 compute_lighting(
    Vec3 point,
    Vec3 normal,
    Vec3 view_dir,
    Material mat)
{
    const Vec3 light_pos = vec3(0.0f, 4.0f, 0.0f);
    const Vec3 light_vec = sub(light_pos, point);
    const float light_dist = len(light_vec);
    const Vec3 light_dir = vdiv(light_vec, light_dist);

    const Vec3 ambient = mul(mat.color, AMBIENT_STRENGTH);

    const float dot_nl = dot(normal, light_dir);
    const float diff = fmaxf(dot_nl, 0.0f);
    const Vec3 diffuse = mul(mat.color, diff * DIFFUSE_STRENGTH);

    const Vec3 reflect_dir = reflect(mul(light_dir, -1.0f), normal);
    float dot_vr = fmaxf(dot(view_dir, reflect_dir), 0.0f);

    // Fast pow^32 using repeated squaring
    float spec = dot_vr;
    spec *= spec;  // 2
    spec *= spec;  // 4
    spec *= spec;  // 8
    spec *= spec;  // 16
    spec *= spec;  // 32

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

bool intersect_sphere(
    Ray ray,
    Sphere sphere,
    HitRecord *rec)
{
    const Vec3 oc = sub(ray.origin, sphere.center);
    const float b = dot(oc, ray.direction);
    const float c = dot(oc, oc) - sphere.radius * sphere.radius;
    const float discriminant = b * b - c;
    if (discriminant < 0) return false;

    const float sqrt_d = sqrtf(discriminant);
    float t = -b - sqrt_d;
    if (t < EPSILON)
    {
        t = -b + sqrt_d;
        if (t < EPSILON) return false;
    }

    if (t >= rec->t) return false;

    rec->hit = true;
    rec->t = t;
    rec->point = add(ray.origin, mul(ray.direction, t));
    rec->normal = vdiv(sub(rec->point, sphere.center), sphere.radius);
    rec->mat = sphere.mat;
    return true;
}

bool intersect_rect(
    Ray ray,
    Rect rect,
    HitRecord *rec)
{
    const float denom = dot(rect.normal, ray.direction);
    if (fabsf(denom) < 0.0001f) return false;

    const float t = dot(sub(rect.point, ray.origin), rect.normal) / denom;
    if (t < EPSILON || t >= rec->t) return false;

    const Vec3 p = add(ray.origin, mul(ray.direction, t));
    const Vec3 v_hit = sub(p, rect.point);

    const float u = dot(v_hit, rect.u);
    const float v = dot(v_hit, rect.v);
    if (u < -rect.width / 2.0f || u > rect.width / 2.0f) return false;
    if (v < -rect.height / 2.0f || v > rect.height / 2.0f) return false;

    rec->hit = true;
    rec->t = t;
    rec->point = p;
    rec->normal = rect.normal;
    rec->mat = rect.mat;
    return true;
}

