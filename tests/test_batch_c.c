/* 用 C 编译器编译本文件，验证批量入口的声明在纯 C11 下可用：
 * size_t / uint8_t / uint32_t 的来源（stddef.h / stdint.h）、__restrict 参数、以及语义本身。
 * 与 tests/test_c_api.c 一样，这个文件只写 C。 */
#include "flashvec.h"

#include <stdio.h>
#include <string.h>

/* 位模式比较（C 版本）。 */
static int bits_eq(float a, float b)
{
    unsigned char ba[4], bb[4];
    memcpy(ba, &a, 4);
    memcpy(bb, &b, 4);
    return memcmp(ba, bb, 4) == 0;
}

int main(void)
{
    fv_vec2 a[8], b[8], o[8], want;
    uint8_t mask[8];
    uint32_t idx[3];
    float ax[5], ay[5], bx[5], by[5], ox[5], oy[5];
    int i, bad = 0;

    for (i = 0; i < 8; ++i) {
        a[i] = fv_vec2_make((fv_real)i, (fv_real)(-2 * i));
        b[i] = fv_vec2_make(1.0f, 0.5f);
        mask[i] = (uint8_t)(i % 2);
    }

    /* n = 5：主循环 4 个 + 尾巴 1 个 */
    fv_vec2_add_x4(a, b, o, 5);
    for (i = 0; i < 5; ++i) {
        if (!bits_eq(o[i].x, a[i].x + 1.0f) || !bits_eq(o[i].y, a[i].y + 0.5f)) {
            bad = 1;
        }
    }

    /* 掩码：mask[i] != 0 做加法，否则写回 a[i] */
    memset(o, 0, sizeof(o));
    fv_vec2_add_x4_mask(a, b, o, mask, 8);
    for (i = 0; i < 8; ++i) {
        want = (mask[i] != 0) ? fv_vec2_add(a[i], b[i]) : a[i];
        if (!bits_eq(o[i].x, want.x) || !bits_eq(o[i].y, want.y)) {
            bad = 1;
        }
    }

    /* 索引：重复 + 乱序；idx 之外的位置保持不动 */
    idx[0] = 7u;
    idx[1] = 0u;
    idx[2] = 7u;
    memcpy(o, a, sizeof(o));
    fv_vec2_add_x4_idx(a, b, o, idx, 3);
    for (i = 0; i < 8; ++i) {
        want = (i == 0 || i == 7) ? fv_vec2_add(a[i], b[i]) : a[i];
        if (!bits_eq(o[i].x, want.x) || !bits_eq(o[i].y, want.y)) {
            bad = 1;
        }
    }

    /* SoA：n = 5 覆盖尾巴 */
    for (i = 0; i < 5; ++i) {
        ax[i] = (fv_real)i * 0.25f;
        ay[i] = (fv_real)(-i);
        bx[i] = 2.0f;
        by[i] = -0.5f;
    }
    fv_vec2_add_soa(ax, ay, bx, by, ox, oy, 5);
    for (i = 0; i < 5; ++i) {
        if (!bits_eq(ox[i], ax[i] + bx[i]) || !bits_eq(oy[i], ay[i] + by[i])) {
            bad = 1;
        }
    }
    {
        const float d = fv_vec2_dot_soa(ax, ay, bx, by, 5);
        double ref = 0.0;
        for (i = 0; i < 5; ++i) {
            ref += (double)ax[i] * (double)bx[i] + (double)ay[i] * (double)by[i];
        }
        if (d != (float)ref) {
            bad = 1; /* n < 8 时是纯顺序累加，必须逐位相等 */
        }
    }

    /* P0 的新符号在 C 下也要能链接 */
    o[0] = fv_vec2_clamp(fv_vec2_make(5.0f, -5.0f), fv_vec2_make(-1.0f, -1.0f), fv_vec2_make(1.0f, 1.0f));
    if (!bits_eq(o[0].x, 1.0f) || !bits_eq(o[0].y, -1.0f)) {
        bad = 1;
    }
    {
        fv_aabb2 box;
        box.lo = fv_vec2_make(0.0f, 0.0f);
        box.hi = fv_vec2_make(1.0f, 1.0f);
        box = fv_aabb2_merge(box, fv_vec2_make(3.0f, -1.0f));
        if (bits_eq(box.hi.x, 3.0f) == 0 || fv_aabb2_overlaps(box, box) != 1) {
            bad = 1;
        }
    }

    printf("C batch API ok: add_x4(5) mask idx soa dot clamp aabb  (bad=%d)\n", bad);
    return bad;
}
