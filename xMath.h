#ifndef XMATH_H
#define XMATH_H

#include <math.h>

typedef struct { float x, y, z; } Vec3;
typedef struct { Vec3 origin, direction; } Ray;

// Constructor
static inline Vec3 vec3(const float x, const float y, const float z) {
    const Vec3 v = {x, y, z};
    return v;
}

// Arithmetic
static inline Vec3 add(Vec3 a, const Vec3 b) {
    a.x += b.x; a.y += b.y; a.z += b.z;
    return a;
}

static inline Vec3 sub(Vec3 a, const Vec3 b) {
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
    return a;
}

static inline Vec3 mul(Vec3 v, const float t) {
    v.x *= t; v.y *= t; v.z *= t;
    return v;
}

static inline Vec3 vdiv(Vec3 v, const float t) {
    float inv = 1.0f / t;
    v.x *= inv; v.y *= inv; v.z *= inv;
    return v;
}

// Dot & Cross
static inline float dot(const Vec3 a, const Vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

static inline Vec3 cross(const Vec3 a, const Vec3 b) {
    return (Vec3){
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

// Length & Normalize
static inline float len(const Vec3 v) {
    return sqrtf(dot(v,v));
}

static inline Vec3 norm(Vec3 v) {
    const float l2 = dot(v,v);
    if (l2 > 0.0f) {
        const float inv = 1.0f / sqrtf(l2);
        v.x *= inv; v.y *= inv; v.z *= inv;
    }
    return v;
}

// Reflect
static inline Vec3 reflect(Vec3 v, const Vec3 n) {
    const float d = 2.0f * dot(v,n);
    v.x -= d * n.x;
    v.y -= d * n.y;
    v.z -= d * n.z;
    return v;
}

#endif // XMATH_H
