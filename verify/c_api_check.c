/* verify/c_api_check.c —— 纯 C11 编译 + 链接 + 运行检查（验证 include/flashvec.h 保持纯 C）
 * 编译: cc -std=c11 -O3 -pedantic-errors -Wall -Wextra -ffp-contract=off -Iinclude c_api_check.c lib -lm
 *
 * 契约（见 include/flashvec.h 注释）：min/max 是逐分量 (a<b)?a:b 语义，**不是** IEEE minNum/fminf；
 * NaN 的吞噬行为因此不同（见文件末尾 INFO 行）。clamp 的期望一律在 lo <= hi 的合法区间上构造。
 */
#include "flashvec.h"

#include <stdio.h>
#include <string.h>

static int g_fail = 0;

static void chk(const char* what, int ok)
{
    printf("  %-46s %s\n", what, ok ? "PASS" : "FAIL");
    if (!ok) g_fail++;
}

/* 用 volatile 中转，避免编译器把期望值本身收缩成 FMA */
static float mul_then_add(float a, float b, float s)
{
    volatile float bb = b, ss = s;
    float p = bb * ss;
    volatile float pp = p;
    return a + pp;
}

int main(void)
{
    /* ---- P0 常量 ---- */
    chk("FV_EPS == 1e-6f", FV_EPS == (fv_real)1e-6f);
    chk("FV_EPS_SQR == 1e-12f 且 == FV_EPS*FV_EPS", FV_EPS_SQR == (fv_real)1e-12f && FV_EPS_SQR == FV_EPS * FV_EPS);

    /* ---- 逐分量 (a<b)?a:b 契约，全部在合法区间 / 非 NaN 数据上 ---- */
    {
        fv_vec2 a = fv_vec2_make(1.0f, 2.0f), b = fv_vec2_make(3.0f, 4.0f);
        fv_vec2 lo = fv_vec2_make(0.0f, 0.0f), hi = fv_vec2_make(10.0f, 10.0f);
        fv_vec2 mn = fv_vec2_min(a, b), mx = fv_vec2_max(a, b);
        fv_vec2 c_in = fv_vec2_clamp(a, lo, hi);                              /* 界内 → 原值 */
        fv_vec2 c_out = fv_vec2_clamp(fv_vec2_make(20.0f, -5.0f), lo, hi);    /* 界外 → 被夹 */
        fv_vec2 md = fv_vec2_madd(a, b, 2.0f);
        chk("vec2 min((1,2),(3,4)) = (1,2)", mn.x == 1.0f && mn.y == 2.0f);
        chk("vec2 max((1,2),(3,4)) = (3,4)", mx.x == 3.0f && mx.y == 4.0f);
        chk("vec2 clamp 界内 (1,2)->(1,2)", c_in.x == 1.0f && c_in.y == 2.0f);
        chk("vec2 clamp 界外 (20,-5)->(10,0)", c_out.x == 10.0f && c_out.y == 0.0f);
        chk("vec2 madd(a,b,2) = (7,10)", md.x == 7.0f && md.y == 10.0f && md.x == mul_then_add(1.0f, 3.0f, 2.0f));

        {
            fv_vec3 p = fv_vec3_make(1.0f, 2.0f, 3.0f), q = fv_vec3_make(3.0f, 2.0f, 1.0f);
            fv_vec3 lo3 = fv_vec3_min(p, q), hi3 = fv_vec3_max(p, q);   /* lo <= hi 保证成立 */
            fv_vec3 in3 = fv_vec3_clamp(p, lo3, hi3);
            fv_vec3 out3 = fv_vec3_clamp(fv_vec3_make(0.0f, 9.0f, 9.0f), lo3, hi3);
            fv_vec3 md3 = fv_vec3_madd(p, q, 0.5f);
            chk("vec3 min/max 逐分量", lo3.x == 1.0f && lo3.y == 2.0f && lo3.z == 1.0f &&
                                       hi3.x == 3.0f && hi3.y == 2.0f && hi3.z == 3.0f);
            chk("vec3 clamp 界内 = 原值", in3.x == p.x && in3.y == p.y && in3.z == p.z);
            chk("vec3 clamp 界外 (0,9,9)->(1,2,3)", out3.x == 1.0f && out3.y == 2.0f && out3.z == 3.0f);
            chk("vec3 madd(p,q,0.5) = (2.5,3,3.5)", md3.x == 2.5f && md3.y == 3.0f && md3.z == 3.5f);
        }

        {
            fv_aabb2 box = { fv_vec2_make(0.0f, 0.0f), fv_vec2_make(1.0f, 1.0f) };
            fv_aabb2 m2 = fv_aabb2_merge(box, fv_vec2_make(-1.0f, 2.0f));
            fv_aabb2 o2 = { fv_vec2_make(0.5f, 0.5f), fv_vec2_make(2.0f, 2.0f) };
            fv_aabb2 far = { fv_vec2_make(9.0f, 9.0f), fv_vec2_make(10.0f, 10.0f) };
            chk("aabb2_merge = [min(lo,p), max(hi,p)]", m2.lo.x == -1.0f && m2.lo.y == 0.0f && m2.hi.x == 1.0f && m2.hi.y == 2.0f);
            chk("aabb2_overlaps 相交=1 / 相离=0", fv_aabb2_overlaps(box, o2) == 1 && fv_aabb2_overlaps(box, far) == 0);
        }
    }

    /* ---- 批量接口：C 侧符号存在性 + 真实调用 + 一处值检查 ---- */
    {
        enum { N = 9 };
        fv_vec2 a[N], b[N], o[N], rl[N];
        float ax[N], ay[N], bx[N], by[N], ox[N], oy[N], lo[N];
        uint8_t mask[N];
        uint32_t idx[N];
        int i, ok = 1;
        for (i = 0; i < N; ++i) {
            a[i] = fv_vec2_make((fv_real)(i + 1), (fv_real)(i + 2));
            b[i] = fv_vec2_make(1.0f, 1.0f);
            rl[i] = fv_vec2_make(a[i].x + 1.0f, a[i].y + 1.0f);
            ax[i] = a[i].x; ay[i] = a[i].y; bx[i] = b[i].x; by[i] = b[i].y;
            mask[i] = (uint8_t)(i % 2);
            idx[i] = (uint32_t)(N - 1 - i);
        }
        fv_vec2_add_x4(a, b, o, N);
        for (i = 0; i < N; ++i) if (o[i].x != rl[i].x || o[i].y != rl[i].y) ok = 0;
        chk("add_x4 全 9 个元素正确（含尾巴）", ok);
        fv_vec2_sub_x4(a, b, o, N);   fv_vec2_scale_x4(a, 2.0f, o, N);
        fv_vec2_lerp_x4(a, b, 0.5f, o, N); fv_vec2_madd_x4(a, b, 2.0f, o, N);
        fv_vec2_add_x4_mask(a, b, o, mask, N);
        ok = 1;
        for (i = 0; i < N; ++i) {
            if (mask[i]) { if (o[i].x != rl[i].x || o[i].y != rl[i].y) ok = 0; }
            else if (o[i].x != a[i].x || o[i].y != a[i].y) ok = 0;   /* 掩码外必须写回 a[i] 原值 */
        }
        chk("add_x4_mask 掩码外逐位等于 a[i]", ok);
        fv_vec2_add_x4_idx(a, b, o, idx, N);
        fv_vec2_scale_x4_idx(a, 2.0f, o, idx, N);
        fv_vec2_add_soa(ax, ay, bx, by, ox, oy, N);
        fv_vec2_scale_soa(ax, ay, 2.0f, ox, oy, N);
        {
            float dot = fv_vec2_dot_soa(ax, ay, bx, by, N);
            fv_vec2_len_soa(ax, ay, lo, N);
            fv_vec2_normalize_soa(ax, ay, ox, oy, N);
            chk("SoA 入口可从纯 C 调用且产出有限值", dot == dot && lo[0] > 0.0f && ox[0] == ox[0]);
        }
    }

    /* ---- 契约差异（只记录，不判定）：min/max 是 (a<b)?a:b，与 IEEE minNum/fminf 的 NaN 语义不同 ---- */
    {
        float nan_v = 0.0f;
        fv_vec2 nn = fv_vec2_make(0.0f / nan_v, 1.0f);           /* (NaN, 1) */
        fv_vec2 zz = fv_vec2_make(0.0f, 0.0f);
        fv_vec2 mn = fv_vec2_min(nn, zz);
        printf("  %-46s INFO  min((NaN,1),(0,0)) = (%.1f,%.1f)：契约 (a<b)?a:b 给出 (0,0)；IEEE fminf 会给 (0,1)\n",
               "min 的 NaN 语义与 fminf 的差异", (double)mn.x, (double)mn.y);
    }

    printf("\nc_api_check: %s\n", g_fail ? "FAIL" : "全部 PASS");
    return g_fail ? 1 : 0;
}
