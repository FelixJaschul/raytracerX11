// ============================================================================
// xBvd.h -> BVH implementation
// ============================================================================

#ifndef XBVH_H
#define XBVH_H

#include <stdlib.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    Vec3 min, max;
} AABB;

typedef struct BVHNode {
    AABB bounds;
    struct BVHNode *left, *right;
    xTriangle *tris;
    xMaterial *mats;
    int count;
} BVHNode;

// Inline helpers
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

typedef struct {
    xTriangle tri; xMaterial mat; Vec3 center;
} Item;

// Straight up ChatGPTed this lol
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
    BVHNode *node = malloc(sizeof(BVHNode));
    node->bounds = box_tri(items[0].tri);
    for (int i = 1; i < n; i++) node->bounds = box_merge(node->bounds, box_tri(items[i].tri));

    if (n <= 4)
    {
        node->count = n;
        node->tris = malloc(n * sizeof(xTriangle));
        node->mats = malloc(n * sizeof(xMaterial));
        for (int i = 0; i < n; i++) { node->tris[i] = items[i].tri; node->mats[i] = items[i].mat; }
        node->left = node->right = NULL;
        return node;
    }

    Vec3 extent = sub(node->bounds.max, node->bounds.min);
    int axis = (extent.y > extent.x) ? 1 : 0;
    if (extent.z > ((float*)&extent)[axis]) axis = 2;

    qsort(items, n, sizeof(Item), axis==0?cmp_x:axis==1?cmp_y:cmp_z);

    const int mid = n / 2;
    node->count = 0; node->tris = node->mats = NULL;
    node->left = build(items, mid);
    node->right = build(items + mid, n - mid);
    return node;
}

static BVHNode* bvh_build(const xModel *models, const int num)
{
    int total = 0;
    for (int i = 0; i < num; i++) total += models[i].num_triangles;
    if (!total) return NULL;

    Item *items = malloc(total * sizeof(Item));
    int idx = 0;
    for (int i = 0; i < num; i++)
    {
        const xModel *m = &models[i];
        for (int j = 0; j < m->num_triangles; j++)
        {
            items[idx].tri = m->transformed_triangles[j];
            items[idx].mat = m->mat;
            xTriangle t = m->transformed_triangles[j];
            items[idx].center = vec3((t.v0.x+t.v1.x+t.v2.x)/3, (t.v0.y+t.v1.y+t.v2.y)/3, (t.v0.z+t.v1.z+t.v2.z)/3);
            idx++;
        }
    }

    BVHNode *root = build(items, total);
    free(items);
    return root;
}

static void bvh_free(BVHNode *n)
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

bool intersect_triangle(const Ray ray, const xTriangle tri, const xMaterial mat, HitRecord *rec)
{
    const Vec3 v0 = tri.v0;
    const Vec3 v1 = tri.v1;
    const Vec3 v2 = tri.v2;

    // Compute edges
    const Vec3 edge1 = sub(v1, v0);
    const Vec3 edge2 = sub(v2, v0);

    // Begin Möller-Trumbore algorithm
    const Vec3 h = cross(ray.direction, edge2);
    const float a = dot(edge1, h);

    // Backface culling: skip back-facing triangles (2x speedup!)
    if (a < EPSILON) return false;

    const float f = 1.0f / a;
    const Vec3 s = sub(ray.origin, v0);
    const float u = f * dot(s, h);

    // Check barycentric coordinate u
    if (u < 0.0f || u > 1.0f) return false;

    const Vec3 q = cross(s, edge1);
    const float v = f * dot(ray.direction, q);

    // Check barycentric coordinate v
    if (v < 0.0f || u + v > 1.0f) return false;

    // Compute intersection distance
    const float t = f * dot(edge2, q);

    // Check if valid and closer than previous hits
    if (t < EPSILON || t >= rec->t) return false;

    // Valid hit! Update hit record
    rec->hit = true;
    rec->t = t;
    rec->point = add(ray.origin, mul(ray.direction, t));
    rec->normal = norm(cross(edge1, edge2));
    rec->mat = mat;

    return true;
}

// Iterative traversal (faster than recursion)
#define STACK_SIZE 64
static inline bool bvh_intersect(BVHNode *root, Ray ray, HitRecord *rec)
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
            {
                if (intersect_triangle(ray, n->tris[i], n->mats[i], rec))
                    hit = true;
            }
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

#ifdef __cplusplus
}
#endif
#endif