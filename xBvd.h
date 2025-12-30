#ifndef XBVD_H
#define XBVD_H

#include <stdlib.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// AABB (Axis-Aligned Bounding Box)

typedef struct {
    Vec3 min;
    Vec3 max;
} AABB;

// Create AABB from a triangle
static inline AABB aabb_from_triangle(const xTriangle tri)
{
    AABB box;
    box.min.x = fminf(fminf(tri.v0.x, tri.v1.x), tri.v2.x);
    box.min.y = fminf(fminf(tri.v0.y, tri.v1.y), tri.v2.y);
    box.min.z = fminf(fminf(tri.v0.z, tri.v1.z), tri.v2.z);
    box.max.x = fmaxf(fmaxf(tri.v0.x, tri.v1.x), tri.v2.x);
    box.max.y = fmaxf(fmaxf(tri.v0.y, tri.v1.y), tri.v2.y);
    box.max.z = fmaxf(fmaxf(tri.v0.z, tri.v1.z), tri.v2.z);
    return box;
}

// Merge two AABBs
static inline AABB aabb_merge(const AABB a, const AABB b)
{
    AABB result;
    result.min.x = fminf(a.min.x, b.min.x);
    result.min.y = fminf(a.min.y, b.min.y);
    result.min.z = fminf(a.min.z, b.min.z);
    result.max.x = fmaxf(a.max.x, b.max.x);
    result.max.y = fmaxf(a.max.y, b.max.y);
    result.max.z = fmaxf(a.max.z, b.max.z);
    return result;
}

// Ray-AABB intersection test (slab method)
static inline bool aabb_intersect(const AABB box, const Ray ray, float t_min, float t_max)
{
    for (int axis = 0; axis < 3; axis++)
    {
        const float inv_d = 1.0f / ((float*)&ray.direction)[axis];
        float t0 = (((float*)&box.min)[axis] - ((float*)&ray.origin)[axis]) * inv_d;
        float t1 = (((float*)&box.max)[axis] - ((float*)&ray.origin)[axis]) * inv_d;

        if (inv_d < 0.0f) {
            float temp = t0;
            t0 = t1;
            t1 = temp;
        }

        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;

        if (t_max <= t_min) return false;
    }
    return true;
}

// Get AABB center
static inline Vec3 aabb_center(const AABB box)
{
    return vec3(
        (box.min.x + box.max.x) * 0.5f,
        (box.min.y + box.max.y) * 0.5f,
        (box.min.z + box.max.z) * 0.5f
    );
}

// BVH (Bounding Volume Hierarchy)

#define BVH_MAX_LEAF_TRIANGLES 4

typedef struct BVHNode {
    AABB bounds;
    struct BVHNode *left;
    struct BVHNode *right;
    xTriangle *triangles;
    xMaterial *materials;
    int num_triangles;
    bool is_leaf;
} BVHNode;

typedef struct {
    xTriangle triangle;
    xMaterial material;
    AABB bounds;
    Vec3 center;
} BVHTriangle;

// Comparison functions for qsort
static int compare_x(const BVHTriangle *a, const BVHTriangle *b)
{
    const BVHTriangle *ta = a;
    const BVHTriangle *tb = b;
    return (ta->center.x > tb->center.x) - (ta->center.x < tb->center.x);
}

static int compare_y(const BVHTriangle *a, const BVHTriangle *b)
{
    const BVHTriangle *ta = a;
    const BVHTriangle *tb = b;
    return (ta->center.y > tb->center.y) - (ta->center.y < tb->center.y);
}

static int compare_z(const BVHTriangle *a, const BVHTriangle *b)
{
    const BVHTriangle *ta = a;
    const BVHTriangle *tb = b;
    return (ta->center.z > tb->center.z) - (ta->center.z < tb->center.z);
}

// Build BVH recursively
static BVHNode* bvh_build_recursive(BVHTriangle *tris, const int num_tris)
{
    BVHNode *node = malloc(sizeof(BVHNode));

    // Compute bounds for this node
    node->bounds = tris[0].bounds;
    for (int i = 1; i < num_tris; i++)
    {
        node->bounds = aabb_merge(node->bounds, tris[i].bounds);
    }

    // Leaf node condition
    if (num_tris <= BVH_MAX_LEAF_TRIANGLES)
    {
        node->is_leaf = true;
        node->num_triangles = num_tris;
        node->triangles = (xTriangle*)malloc(num_tris * sizeof(xTriangle));
        node->materials = (xMaterial*)malloc(num_tris * sizeof(xMaterial));

        for (int i = 0; i < num_tris; i++)
        {
            node->triangles[i] = tris[i].triangle;
            node->materials[i] = tris[i].material;
        }

        node->left = NULL;
        node->right = NULL;
        return node;
    }

    // Interior node - split along longest axis
    Vec3 extent = sub(node->bounds.max, node->bounds.min);
    int axis = 0;
    if (extent.y > extent.x) axis = 1;
    if (extent.z > ((float*)&extent)[axis]) axis = 2;

    // Sort triangles along chosen axis
    if (axis == 0) qsort(tris, num_tris, sizeof(BVHTriangle), compare_x);
    else if (axis == 1) qsort(tris, num_tris, sizeof(BVHTriangle), compare_y);
    else qsort(tris, num_tris, sizeof(BVHTriangle), compare_z);

    // Split at median
    int mid = num_tris / 2;

    node->is_leaf = false;
    node->num_triangles = 0;
    node->triangles = NULL;
    node->materials = NULL;
    node->left = bvh_build_recursive(tris, mid);
    node->right = bvh_build_recursive(tris + mid, num_tris - mid);

    return node;
}

// Build BVH from scene models
static BVHNode* bvh_build(const xModel *models, const int num_models)
{
    // Count total triangles
    int total_triangles = 0;
    for (int i = 0; i < num_models; i++)
    {
        total_triangles += models[i].num_triangles;
    }

    if (total_triangles == 0) return NULL;

    // Collect all triangles with their materials and bounds
    BVHTriangle *tris = malloc(total_triangles * sizeof(BVHTriangle));
    int tri_idx = 0;

    for (int i = 0; i < num_models; i++)
    {
        const xModel *model = &models[i];
        for (int j = 0; j < model->num_triangles; j++)
        {
            tris[tri_idx].triangle = model->transformed_triangles[j];
            tris[tri_idx].material = model->mat;
            tris[tri_idx].bounds = aabb_from_triangle(model->transformed_triangles[j]);
            tris[tri_idx].center = aabb_center(tris[tri_idx].bounds);
            tri_idx++;
        }
    }

    BVHNode *root = bvh_build_recursive(tris, total_triangles);
    free(tris);
    return root;
}

// Free BVH tree
static void bvh_free(BVHNode *node)
{
    if (!node) return;

    if (node->is_leaf)
    {
        free(node->triangles);
        free(node->materials);
    }
    else
    {
        bvh_free(node->left);
        bvh_free(node->right);
    }

    free(node);
}

#ifdef __cplusplus
}
#endif

#endif // XBVD_H