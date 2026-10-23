/* bench/ref.cpp -- 严格 IEEE 语义参考实现。见 ref.h 顶部说明。 */
#include "ref.h"

#include <math.h>

#define REF_EPS ((float)1e-6)

void ref_add(const float *ax, const float *ay, const float *bx, const float *by,
             float *ox, float *oy, unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        ox[i] = ax[i] + bx[i];
        oy[i] = ay[i] + by[i];
    }
}

void ref_scale(const float *ax, const float *ay, float s, float *ox, float *oy, unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        ox[i] = ax[i] * s;
        oy[i] = ay[i] * s;
    }
}

float ref_dot_sum(const float *ax, const float *ay, const float *bx, const float *by, unsigned n)
{
    float s = 0.0f;
    for (unsigned i = 0; i < n; ++i) {
        s = s + (ax[i] * bx[i] + ay[i] * by[i]);
    }
    return s;
}

void ref_dot_elem(const float *ax, const float *ay, const float *bx, const float *by,
                  float *ox, float *oy, unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        const float d = ax[i] * bx[i] + ay[i] * by[i];
        ox[i] = d;
        oy[i] = d;
    }
}

void ref_len(const float *ax, const float *ay, float *ox, float *oy, unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        const float l2 = ax[i] * ax[i] + ay[i] * ay[i];
        float l = sqrtf(l2);
        if (l <= REF_EPS) {
            l = 0.0f;
        }
        ox[i] = l;
        oy[i] = 0.0f;
    }
}

void ref_normalize(const float *ax, const float *ay, float *ox, float *oy, unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        const float l2 = ax[i] * ax[i] + ay[i] * ay[i];
        const float l = sqrtf(l2);
        if (l <= REF_EPS) {
            ox[i] = 0.0f;
            oy[i] = 0.0f;
        } else {
            const float inv = 1.0f / l;
            ox[i] = ax[i] * inv;
            oy[i] = ay[i] * inv;
        }
    }
}

void ref_collision(const float *px, const float *py, const float *qx, const float *qy,
                   const float *vx, const float *vy, float *ox, float *oy, unsigned n)
{
    for (unsigned i = 0; i < n; ++i) {
        const float dx = px[i] - qx[i];
        const float dy = py[i] - qy[i];
        const float l2 = dx * dx + dy * dy;
        const float l = sqrtf(l2);
        float nx = 0.0f;
        float ny = 0.0f;
        if (l > REF_EPS) {
            const float inv = 1.0f / l;
            nx = dx * inv;
            ny = dy * inv;
        }
        const float dn = vx[i] * nx + vy[i] * ny;
        const float k = 2.0f * dn;
        ox[i] = vx[i] - k * nx;
        oy[i] = vy[i] - k * ny;
    }
}
