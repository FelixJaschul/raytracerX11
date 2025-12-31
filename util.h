// ============================================================================
// xBvd.h - Bounding Volume Hierarchy for fast ray-triangle intersection
// ============================================================================
#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <math.h>

#ifndef EPSILON
#define EPSILON 0.0001f
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Axis-Aligned Bounding Box
typedef struct {
    Vec3 min, max;
} AABB;

// BVH Tree Node
typedef struct BVHNode {
    AABB bounds;
    struct BVHNode *left, *right;
    xTriangle *tris;
    xMaterial *mats;
    int count;
} BVHNode;

// Build BVH from array of models
/*  -> Example:
 *  BVHNode* root = bvh_build(models, num_models);
 */
BVHNode* bvh_build(const xModel *models, int num);

// Free all resources used by BVH
/*  -> Example:
 *  bvh_free(root);
 */
void bvh_free(BVHNode *n);

// Intersect ray with BVH, returns true if hit found
/*  -> Example:
 *  HitRecord rec = { .hit = false, .t = 1e30f };
 *  if (bvh_intersect(root, ray, &rec)) { ... }
 */
bool bvh_intersect(BVHNode *root, Ray ray, HitRecord *rec);

#ifdef __cplusplus
}
#endif

#ifdef UTIL_IMPLEMENTATION

typedef struct {
    xTriangle tri; xMaterial mat; Vec3 center;
} Item;

// Internal helpers
static inline AABB box_tri(const xTriangle t)
{
    AABB b;
    b.min.x = fminf(fminf(t.v0.x,t.v1.x),t.v2.x);
    b.min.y = fminf(fminf(t.v0.y,t.v1.y),t.v2.y);
    b.min.z = fminf(fminf(t.v0.z,t.v1.z),t.v2.z);
    b.max.x = fmaxf(fmaxf(t.v0.x,t.v1.x),t.v2.x);
    b.max.y = fmaxf(fmaxf(t.v0.y,t.v1.y),t.v2.y);
    b.max.z = fmaxf(fmaxf(t.v0.z,t.v1.z),t.v2.z);
    return b;
}

static inline AABB box_merge(const AABB a, const AABB b)
{
    AABB r;
    r.min.x=fminf(a.min.x,b.min.x); r.min.y=fminf(a.min.y,b.min.y); r.min.z=fminf(a.min.z,b.min.z);
    r.max.x=fmaxf(a.max.x,b.max.x); r.max.y=fmaxf(a.max.y,b.max.y); r.max.z=fmaxf(a.max.z,b.max.z);
    return r;
}

#ifndef __ARM_NEON__
static inline bool box_hit(AABB box, Ray ray, float tmin, float tmax)
{
    for (int i = 0; i < 3; i++)
    {
        const float inv = 1.0f / ((float*)&ray.direction)[i];
        float t0 = (((float*)&box.min)[i] - ((float*)&ray.origin)[i]) * inv;
        float t1 = (((float*)&box.max)[i] - ((float*)&ray.origin)[i]) * inv;
        if (inv < 0)
        {
            const float tmp = t0;
            t0 = t1;
            t1 = tmp;
        }

        tmin = t0 > tmin ? t0 : tmin;
        tmax = t1 < tmax ? t1 : tmax;
        if (tmax <= tmin) return false;
    }
    return true;
}
#endif

#ifdef __ARM_NEON__
#include <arm_neon.h>

static inline bool box_hit(const AABB box, const Ray ray, const float tmin, const float tmax)
{
    float ray_o[4]   = {ray.origin.x, ray.origin.y, ray.origin.z, 0};
    float ray_d[4]   = {ray.direction.x, ray.direction.y, ray.direction.z, 1};
    float box_min[4] = {box.min.x, box.min.y, box.min.z, 0};
    float box_max[4] = {box.max.x, box.max.y, box.max.z, 0};

    const float32x4_t ray_orig = vld1q_f32(ray_o);
    const float32x4_t ray_dir = vld1q_f32(ray_d);
    const float32x4_t bmin = vld1q_f32(box_min);
    const float32x4_t bmax = vld1q_f32(box_max);

    float32x4_t inv_dir = vrecpeq_f32(ray_dir);
    inv_dir = vmulq_f32(vrecpsq_f32(ray_dir, inv_dir), inv_dir);

    const float32x4_t t0 = vmulq_f32(vsubq_f32(bmin, ray_orig), inv_dir);
    const float32x4_t t1 = vmulq_f32(vsubq_f32(bmax, ray_orig), inv_dir);

    const float32x4_t t_near = vminq_f32(t0, t1);
    const float32x4_t t_far = vmaxq_f32(t0, t1);

    const float t_entry = fmaxf(fmaxf(vgetq_lane_f32(t_near,0), vgetq_lane_f32(t_near,1)), fmaxf(vgetq_lane_f32(t_near,2), tmin));
    const float t_exit  = fminf(fminf(vgetq_lane_f32(t_far ,0), vgetq_lane_f32(t_far ,1)), fminf(vgetq_lane_f32(t_far ,2), tmax));

    return t_exit > t_entry;
}
#endif

static int cmp_x(const void *a, const void *b)
{
    return (((Item*)a)->center.x > ((Item*)b)->center.x) - (((Item*)a)->center.x < ((Item*)b)->center.x);
}
static int cmp_y(const void *a, const void *b)
{
    return (((Item*)a)->center.y > ((Item*)b)->center.y) - (((Item*)a)->center.y < ((Item*)b)->center.y);
}
static int cmp_z(const void *a, const void *b)
{
    return (((Item*)a)->center.z > ((Item*)b)->center.z) - (((Item*)a)->center.z < ((Item*)b)->center.z);
}

static BVHNode* build(Item *items, const int n)
{
    BVHNode *node = (BVHNode*)malloc(sizeof(BVHNode));
    node->bounds = box_tri(items[0].tri);
    for (int i = 1; i < n; i++) node->bounds = box_merge(node->bounds, box_tri(items[i].tri));

    if (n <= 4)
    {
        node->count = n;
        node->tris = (xTriangle*)malloc(n * sizeof(xTriangle));
        node->mats = (xMaterial*)malloc(n * sizeof(xMaterial));
        for (int i = 0; i < n; i++)
        {
            node->tris[i] = items[i].tri;
            node->mats[i] = items[i].mat;
        }
        node->left = node->right = NULL;
        return node;
    }

    Vec3 extent = sub(node->bounds.max, node->bounds.min);
    int axis = (extent.y > extent.x) ? 1 : 0;
    if (extent.z > ((float*)&extent)[axis]) axis = 2;

    qsort(items, n, sizeof(Item), axis==0?cmp_x:axis==1?cmp_y:cmp_z);

    const int mid = n / 2;
    node->count = 0; node->tris = node->mats = NULL;
    node->left  = build(items, mid);
    node->right = build(items + mid, n - mid);
    return node;
}

inline BVHNode* bvh_build(const xModel *models, const int num)
{
    int total = 0;
    for (int i = 0; i < num; i++) total += models[i].num_triangles;
    if (!total) return NULL;

    Item *items = (Item*)malloc(total * sizeof(Item));
    int idx = 0;
    for (int i = 0; i < num; i++)
    {
        const xModel *m = &models[i];
        for (int j = 0; j < m->num_triangles; j++)
        {
            items[idx].tri = m->transformed_triangles[j];
            items[idx].mat = m->mat;
            xTriangle t = m->transformed_triangles[j];
            items[idx].center = vec3((t.v0.x+t.v1.x+t.v2.x)/3.0f, (t.v0.y+t.v1.y+t.v2.y)/3.0f, (t.v0.z+t.v1.z+t.v2.z)/3.0f);
            idx++;
        }
    }

    BVHNode *root = build(items, total);
    free(items);
    return root;
}

inline void bvh_free(BVHNode *n)
{
    if (!n) return;
    if (n->count)
    {
        free(n->tris);
        free(n->mats);
    }
    else
    {
        bvh_free(n->left);
        bvh_free(n->right);
    }
    free(n);
}

static inline bool intersect_triangle(const Ray ray, const xTriangle tri, const xMaterial mat, HitRecord *rec)
{
    const Vec3 v0 = tri.v0;
    const Vec3 v1 = tri.v1;
    const Vec3 v2 = tri.v2;

    const Vec3 edge1 = sub(v1, v0);
    const Vec3 edge2 = sub(v2, v0);

    const Vec3 h = cross(ray.direction, edge2);
    const float a = dot(edge1, h);

    // if (fabsf(a) < EPSILON) return false;

    const float f = 1.0f / a;
    const Vec3 s = sub(ray.origin, v0);
    const float u = f * dot(s, h);

    const Vec3 q = cross(s, edge1);
    const float v = f * dot(ray.direction, q);

    if (u < 0.0f || u > 1.0f)     return false;
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

#define STACK_SIZE 64
inline bool bvh_intersect(BVHNode *root, Ray ray, HitRecord *rec)
{
    if (!root) return false;
    BVHNode *stack[STACK_SIZE];
    int sp = 0;
    stack[sp++] = root;
    bool hit = false;

    while (sp > 0)
    {
        BVHNode *n = stack[--sp];
        if (!box_hit(n->bounds, ray, 0.001f, rec->t)) continue;

        if (n->count)
        {
            for (int i = 0; i < n->count; i++)
                if (intersect_triangle(ray, n->tris[i], n->mats[i], rec)) hit = true;
        }
        else
        {
            if (sp + 2 < STACK_SIZE)
            {
                stack[sp++] = n->left;
                stack[sp++] = n->right;
            }
        }
    }
    return hit;
}

#endif // UTIL_IMPLEMENTATION
#endif // UTIL_H