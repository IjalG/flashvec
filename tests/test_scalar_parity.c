/* Scalar API bit-level parity: default (header static inline) vs FV_USE_EXTERN (library symbols).
 *
 * The same source is compiled twice; every call's raw output bytes are written to a file, and the two
 * files are compared with cmp. This is the P0 counterpart of tests/test_parity.c (which covers batch).
 * Inputs include zero / tiny / singular / non-finite cases so the guards are exercised too.
 *
 * Usage: test_scalar_parity <output.bin>
 */
#include "flashvec.h"

#include <stdio.h>
#include <string.h>

static FILE *g_out;

static void wf(float x)
{
    fwrite(&x, sizeof x, 1, g_out);
}

static void wv2(fv_vec2 v)
{
    wf(v.x); wf(v.y);
}

static void wv3(fv_vec3 v)
{
    wf(v.x); wf(v.y); wf(v.z);
}

static void wv4(fv_vec4 v)
{
    wf(v.x); wf(v.y); wf(v.z); wf(v.w);
}

static void wq(fv_quat q)
{
    wf(q.x); wf(q.y); wf(q.z); wf(q.w);
}

static void wm3(fv_mat3 m)
{
    int i;
    for (i = 0; i < 9; ++i) {
        wf(m.m[i]);
    }
}

static void wm4(fv_mat4 m)
{
    int i;
    for (i = 0; i < 16; ++i) {
        wf(m.m[i]);
    }
}

int main(int argc, char **argv)
{
    fv_vec2 a2, b2, zero2, tiny2, onto2;
    fv_vec3 a3, b3, axis3;
    fv_vec4 a4, b4;
    fv_quat q1, q2, qz;
    fv_mat3 m3, m3b, m3s, m3q;
    fv_mat4 m4, m4b, m4s;

    if (argc != 2) {
        fprintf(stderr, "usage: %s <output.bin>\n", argv[0]);
        return 2;
    }
    g_out = fopen(argv[1], "wb");
    if (g_out == NULL) {
        perror(argv[1]);
        return 2;
    }

    a2 = fv_vec2_make(1.5f, -2.25f);
    b2 = fv_vec2_make(0.125f, 3.5f);
    zero2 = fv_vec2_make(0.0f, 0.0f);
    tiny2 = fv_vec2_make(1e-7f, -2e-7f);
    onto2 = fv_vec2_make(3.0f, -4.0f);

    wv2(fv_vec2_add(a2, b2));
    wv2(fv_vec2_sub(a2, b2));
    wv2(fv_vec2_neg(a2));
    wv2(fv_vec2_scale(a2, -3.25f));
    wv2(fv_vec2_lerp(a2, b2, 0.0f));
    wv2(fv_vec2_lerp(a2, b2, 1.0f));
    wv2(fv_vec2_lerp(a2, b2, -0.5f));
    wv2(fv_vec2_min(a2, b2));
    wv2(fv_vec2_max(a2, b2));
    wv2(fv_vec2_clamp(a2, fv_vec2_make(-1.0f, -1.0f), fv_vec2_make(1.0f, 1.0f)));
    wv2(fv_vec2_madd(a2, b2, 0.75f));
    wf(fv_vec2_dot(a2, b2));
    wf(fv_vec2_len_sqr(a2));
    wf(fv_vec2_len(a2));
    wf(fv_vec2_dist_sqr(a2, b2));
    wf(fv_vec2_dist(a2, b2));
    wv2(fv_vec2_normalize(a2));
    wv2(fv_vec2_normalize(zero2));
    wv2(fv_vec2_normalize(tiny2));
    wv2(fv_vec2_project(a2, onto2));
    wv2(fv_vec2_project(a2, zero2));
    wf((float)fv_vec2_near(a2, a2, 0.0f));
    wf((float)fv_vec2_near(a2, b2, 0.0f));
    wv2(fv_vec2_addmul(a2, b2, 0.5f));
    wv2(fv_vec2_maddv(a2, b2, a2));

    a3 = fv_vec3_make(1.5f, -2.25f, 0.75f);
    b3 = fv_vec3_make(-0.5f, 2.0f, 3.0f);
    wv3(fv_vec3_add(a3, b3));
    wv3(fv_vec3_sub(a3, b3));
    wv3(fv_vec3_scale(a3, 2.5f));
    wv3(fv_vec3_lerp(a3, b3, 0.25f));
    wf(fv_vec3_dot(a3, b3));
    wf(fv_vec3_len(a3));
    wv3(fv_vec3_normalize(a3));
    wv3(fv_vec3_normalize(fv_vec3_make(0.0f, 0.0f, 0.0f)));
    wv3(fv_vec3_cross(a3, b3));
    wv3(fv_vec3_project(a3, b3));
    wv3(fv_vec3_addmul(a3, b3, 0.5f));

    a4 = fv_vec4_make(1.5f, -2.25f, 0.75f, 0.5f);
    b4 = fv_vec4_make(-0.5f, 2.0f, 3.0f, -1.0f);
    wv4(fv_vec4_add(a4, b4));
    wv4(fv_vec4_sub(a4, b4));
    wv4(fv_vec4_scale(a4, -1.25f));
    wv4(fv_vec4_lerp(a4, b4, 0.75f));
    wf(fv_vec4_dot(a4, b4));
    wf(fv_vec4_len(a4));
    wv4(fv_vec4_normalize(a4));
    wv4(fv_vec4_normalize(fv_vec4_make(0.0f, 0.0f, 0.0f, 0.0f)));
    wv4(fv_vec4_project(a4, b4));
    wv4(fv_vec4_addmul(a4, b4, 0.5f));

    axis3 = fv_vec3_make(0.0f, 0.0f, 2.0f);
    q1 = fv_quat_from_axis_angle(axis3, 0.75f);
    q2 = fv_quat_from_euler_xyz(0.25f, -0.5f, 0.75f);
    qz = fv_quat_make(0.0f, 0.0f, 0.0f, 0.0f);
    wq(fv_quat_add(q1, q2));
    wq(fv_quat_sub(q1, q2));
    wq(fv_quat_scale(q1, -2.0f));
    wq(fv_quat_mul(q1, q2));
    wf(fv_quat_dot(q1, q2));
    wf(fv_quat_len(q1));
    wq(fv_quat_normalize(q1));
    wq(fv_quat_normalize(qz));
    wq(fv_quat_conjugate(q1));
    wq(fv_quat_inverse(q1));
    wq(fv_quat_inverse(qz));
    wv3(fv_quat_rotate_vec3(q1, a3));
    wq(fv_quat_nlerp(q1, q2, 0.25f));
    wq(fv_quat_slerp(q1, q2, 0.25f));
    wq(fv_quat_slerp(q1, q1, 0.5f));

    m3 = fv_mat3_from_rows(fv_vec3_make(1.0f, 0.5f, 0.25f),
                           fv_vec3_make(-0.5f, 2.0f, 0.75f),
                           fv_vec3_make(0.25f, -0.75f, 1.5f));
    m3b = fv_mat3_from_rows(fv_vec3_make(-1.0f, 0.25f, 0.5f),
                            fv_vec3_make(0.75f, 0.5f, -0.25f),
                            fv_vec3_make(0.5f, -1.0f, 2.0f));
    m3s = fv_mat3_from_rows(fv_vec3_make(1.0f, 2.0f, 3.0f),
                            fv_vec3_make(2.0f, 4.0f, 6.0f),
                            fv_vec3_make(3.0f, 6.0f, 9.0f));
    m3q = fv_mat3_from_quat(q1);
    wv3(fv_mat3_transform_vec3(&m3, a3));
    wm3(fv_mat3_mul(&m3, &m3b));
    wm3(fv_mat3_transpose(&m3));
    wm3(m3q);
    wf(fv_mat3_determinant(&m3));
    wf(fv_mat3_determinant(&m3s));
    wm3(fv_mat3_inverse(&m3));
    wm3(fv_mat3_inverse(&m3s));
    wq(fv_quat_from_mat3(&m3q));

    m4 = fv_mat4_identity();
    m4.m[3] = 1.5f; m4.m[7] = -2.0f; m4.m[11] = 0.75f;
    m4.m[12] = 0.25f; m4.m[13] = -0.5f; m4.m[14] = 2.0f; m4.m[15] = 1.0f;
    m4b = fv_mat4_identity();
    m4b.m[0] = 2.0f; m4b.m[5] = 0.5f; m4b.m[10] = -1.5f; m4b.m[15] = 1.0f;
    m4s = fv_mat4_identity();
    m4s.m[0] = 0.0f; m4s.m[5] = 0.0f; m4s.m[10] = 0.0f; m4s.m[15] = 0.0f;
    wm4(fv_mat4_mul(&m4, &m4b));
    wm4(fv_mat4_inverse(&m4));
    wm4(fv_mat4_inverse(&m4s));

    wf(fv_atan_approx(-1.0f));
    wf(fv_atan_approx(-0.5f));
    wf(fv_atan_approx(0.0f));
    wf(fv_atan_approx(0.5f));
    wf(fv_atan_approx(1.0f));

    fclose(g_out);
    return 0;
}
