// bench_all.cpp — flashvec vs GLM 全量 API 基准
//
// 口径（三类，结论不可混用）：
//   P  = 标量 API 逐元素循环。不加 LTO 时 flashvec 标量是跨 TU 调用（真实默认形态，对它不利）；
//        加 -flto 后编译器可内联，此时才是与 GLM 模板真正的对等口径。
//   B  = batch / SoA / stream / fast 入口。flashvec 头文件内 static inline，与 GLM 同 TU 内联，天然对等。
//   S  = batch 入口在小 n 下的调用开销。
//   AA = 同一实现两个 TU 互比，给出本机噪声下限。
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>

#define GLM_ENABLE_EXPERIMENTAL
#include "flashvec.h"
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>   // nlerp / rotate(vec3,quat)
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/epsilon.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtc/type_ptr.hpp>

void aa_add4_a(const fv_vec4 *a, const fv_vec4 *b, fv_vec4 *o, size_t n);
void aa_add4_b(const fv_vec4 *a, const fv_vec4 *b, fv_vec4 *o, size_t n);

// ------------------------- 数据 -------------------------
static float *gA, *gB, *gC, *gO;                  // AoS，容量 gN 个 vec4
static float *SA[4], *SB[4], *SO[4];               // SoA 4 通道
static float *MA[9], *MB[9], *MO[9];               // mat3 SoA 9 通道
static uint8_t *gMask; static uint32_t *gIdx; static float *gAng;
static fv_quat *gQA, *gQB; static glm::quat *GQA, *GQB;
static size_t gN, gM;
static volatile float gSink;
static fv_mat4 gFvM4; static glm::mat4 gGlM4;
static fv_mat3 gFvM3; static glm::mat3 gGlM3;
static fv_quat gFvQ;  static glm::quat gGlQ;

static void *a64(size_t b) { void *p=nullptr; if(posix_memalign(&p,64,b)) {perror("alloc");exit(1);} memset(p,0,b); return p; }
static float frand(unsigned *s) { *s=*s*1103515245u+12345u; return 0.1f+(float)((*s>>8)&0xFFFFFF)/(float)0xFFFFFF; }
static void fill(float *p, size_t n) { unsigned s=987654321u; for(size_t i=0;i<n;i++) p[i]=frand(&s); }

typedef void (*Fn)();
static double now_ms(){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec*1e3+t.tv_nsec/1e6; }

// ------------------------- 形态 runner -------------------------
template<class T,class F> static void p1(F f){ const T*a=(const T*)gA; T*o=(T*)gO; for(size_t i=0;i<gN;i++) o[i]=f(a[i]); gSink+=((const float*)o)[0]; }
template<class T,class F> static void p2(F f){ const T*a=(const T*)gA; const T*b=(const T*)gB; T*o=(T*)gO; for(size_t i=0;i<gN;i++) o[i]=f(a[i],b[i]); gSink+=((const float*)o)[0]; }
template<class T,class F> static void p3(F f){ const T*a=(const T*)gA; const T*b=(const T*)gB; const T*c=(const T*)gC; T*o=(T*)gO; for(size_t i=0;i<gN;i++) o[i]=f(a[i],b[i],c[i]); gSink+=((const float*)o)[0]; }
template<class T,class F> static void p1s(F f){ const T*a=(const T*)gA; T*o=(T*)gO; for(size_t i=0;i<gN;i++) o[i]=f(a[i],1.5f); gSink+=((const float*)o)[0]; }
template<class T,class F> static void p2s(F f){ const T*a=(const T*)gA; const T*b=(const T*)gB; T*o=(T*)gO; for(size_t i=0;i<gN;i++) o[i]=f(a[i],b[i],0.25f); gSink+=((const float*)o)[0]; }
template<class T,class F> static void s1(F f){ const T*a=(const T*)gA; float s=0; for(size_t i=0;i<gN;i++) s+=(float)f(a[i]); gSink+=s; }
template<class T,class F> static void s2(F f){ const T*a=(const T*)gA; const T*b=(const T*)gB; float s=0; for(size_t i=0;i<gN;i++) s+=(float)f(a[i],b[i]); gSink+=s; }
template<class T,class F> static void si2(F f){ const T*a=(const T*)gA; const T*b=(const T*)gB; float s=0; for(size_t i=0;i<gN;i++) s+=(float)f(a[i],b[i]); gSink+=s; }
template<class TI,class TO,class F> static void x1(F f){ const TI*a=(const TI*)gA; TO*o=(TO*)gO; for(size_t i=0;i<gN;i++) o[i]=f(a[i]); gSink+=((const float*)o)[0]; }
template<class TI,class TO,class F> static void xq1(F f){ const TI*a=(const TI*)gA; TO*o=(TO*)gO; for(size_t i=0;i<gM;i++) o[i]=f(a[i]); gSink+=((const float*)o)[0]; }
template<class T,class F> static void pm1(F f){ const T*a=(const T*)gA; T*o=(T*)gO; for(size_t i=0;i<gM;i++) o[i]=f(a[i]); gSink+=((const float*)o)[0]; }
template<class T,class F> static void pm2(F f){ const T*a=(const T*)gA; const T*b=(const T*)gB; T*o=(T*)gO; for(size_t i=0;i<gM;i++) o[i]=f(a[i],b[i]); gSink+=((const float*)o)[0]; }
template<class T,class F> static void sm1(F f){ const T*a=(const T*)gA; float s=0; for(size_t i=0;i<gM;i++) s+=(float)f(a[i]); gSink+=s; }
// 常量算子 + vec 数组
template<class T,class F> static void c1(F f){ const T*a=(const T*)gA; T*o=(T*)gO; for(size_t i=0;i<gN;i++) o[i]=f(a[i]); gSink+=((const float*)o)[0]; }
// quat（独立缓冲与布局）
template<class T,class F> static void q1(F f){ const T*a=gQA; T*o=(T*)gO; for(size_t i=0;i<gN;i++) o[i]=f(a[i]); gSink+=((const float*)o)[0]; }
template<class T,class F> static void q2(F f){ const T*a=gQA; const T*b=(const T*)gQB; T*o=(T*)gO; for(size_t i=0;i<gN;i++) o[i]=f(a[i],b[i]); gSink+=((const float*)o)[0]; }
template<class T,class F> static void q2s(F f){ const T*a=gQA; const T*b=gQB; T*o=(T*)gO; for(size_t i=0;i<gN;i++) o[i]=f(a[i],b[i],0.25f); gSink+=((const float*)o)[0]; }
template<class T,class F> static void qs1(F f){ const T*a=gQA; float s=0; for(size_t i=0;i<gN;i++) s+=(float)f(a[i]); gSink+=s; }
template<class T,class F> static void qs2(F f){ const T*a=gQA; const T*b=gQB; float s=0; for(size_t i=0;i<gN;i++) s+=(float)f(a[i],b[i]); gSink+=s; }
// GLM 侧 quat 用 GQA/GQB
template<class F> static void gq1(F f){ const glm::quat*a=GQA; glm::quat*o=(glm::quat*)gO; for(size_t i=0;i<gN;i++) o[i]=f(a[i]); gSink+=((const float*)o)[0]; }
template<class F> static void gq2(F f){ const glm::quat*a=GQA; const glm::quat*b=GQB; glm::quat*o=(glm::quat*)gO; for(size_t i=0;i<gN;i++) o[i]=f(a[i],b[i]); gSink+=((const float*)o)[0]; }
template<class F> static void gq2s(F f){ const glm::quat*a=GQA; const glm::quat*b=GQB; glm::quat*o=(glm::quat*)gO; for(size_t i=0;i<gN;i++) o[i]=f(a[i],b[i],0.25f); gSink+=((const float*)o)[0]; }
template<class F> static void gqs1(F f){ const glm::quat*a=GQA; float s=0; for(size_t i=0;i<gN;i++) s+=(float)f(a[i]); gSink+=s; }
template<class F> static void gqs2(F f){ const glm::quat*a=GQA; const glm::quat*b=GQB; float s=0; for(size_t i=0;i<gN;i++) s+=(float)f(a[i],b[i]); gSink+=s; }

struct GABB { glm::vec2 lo, hi; };   // 与 fv_aabb2 同布局

// ------------------------- 计时 -------------------------
static double med(std::vector<double> v){ if(v.empty())return 0./0.; std::sort(v.begin(),v.end()); size_t n=v.size(); return n&1?v[n/2]:0.5*(v[n/2-1]+v[n/2]); }
static double qt(std::vector<double> v,double q){ if(v.empty())return 0./0.; std::sort(v.begin(),v.end()); return v[(size_t)(q*(v.size()-1)+0.5)]; }
struct Case { const char *grp, *name; Fn a, b; };
static void run(const Case &c, int reps){
    c.a(); c.b(); c.a(); c.b();
    std::vector<double> ta,tb,ra;
    for(int r=0;r<reps;r++){ double t0,t1,t2,t3;
        if(r&1){ t0=now_ms(); c.b(); t1=now_ms(); t2=now_ms(); c.a(); t3=now_ms(); tb.push_back(t1-t0); ta.push_back(t3-t2); }
        else   { t0=now_ms(); c.a(); t1=now_ms(); t2=now_ms(); c.b(); t3=now_ms(); ta.push_back(t1-t0); tb.push_back(t3-t2); }
        if(tb.back()>0&&ta.back()>0) ra.push_back(tb.back()/ta.back()); }
    double ma=med(ta), mb=med(tb);
    printf("N=%zu\t%-3s\t%-26s\tA_ms=%.4f\tB_ms=%.4f\tratio=%.4f\tr10=%.4f\tr90=%.4f\n",
           gN,c.grp,c.name,ma,mb,mb/ma,qt(ra,0.10),qt(ra,0.90)); fflush(stdout);
}

// ------------------------- 校验 -------------------------
static void verify(){
    float d4=0,d3=0,dq=0;
    {   fv_vec4 iv={0.7f,0.3f,1.1f,1.0f},of; fv_mat4_transform_points(&gFvM4,&iv,&of,1);
        glm::vec4 og=gGlM4*glm::vec4(0.7f,0.3f,1.1f,1.0f);
        for(int k=0;k<3;k++) d4=fmaxf(d4,fabsf((&of.x)[k]-glm::value_ptr(og)[k])); }
    {   fv_vec3 v={0.4f,0.9f,0.2f}; fv_vec3 o=fv_mat3_transform_vec3(&gFvM3,v);
        glm::vec3 g=gGlM3*glm::vec3(0.4f,0.9f,0.2f);
        for(int k=0;k<3;k++) d3=fmaxf(d3,fabsf((&o.x)[k]-glm::value_ptr(g)[k])); }
    {   fv_vec3 v={0.4f,0.9f,0.2f}; fv_vec3 o=fv_quat_rotate_vec3(gFvQ,v);
        glm::vec3 g=glm::rotate(gGlQ,glm::vec3(0.4f,0.9f,0.2f));
        for(int k=0;k<3;k++) dq=fmaxf(dq,fabsf((&o.x)[k]-glm::value_ptr(g)[k])); }
    fprintf(stderr,"[verify] mat4=%.1e mat3=%.1e quat_rot=%.1e  %s\n", d4,d3,dq,
            (d4<1e-4f&&d3<1e-4f)?"矩阵侧 MATCH":"矩阵侧 MISMATCH");
    if(dq>=1e-4f) fprintf(stderr,"[verify] 四元数乘法/旋转约定两侧不同(fv 与 GLM 的乘法序或共轭方向不一致)，quat 组的比值只反映等量运算开销，不是同一旋转的复现\n");
}

// ------------------------- 注册表 -------------------------
static std::vector<Case> build_cases(){
    std::vector<Case> v;
    auto add=[&](const char*g,const char*n,Fn a,Fn b){ v.push_back({g,n,a,b}); };
    // ---- vec2 标量 ----
    add("P","vec2_add", +[]{ p2<fv_vec2>([](fv_vec2 a, fv_vec2 b){ return fv_vec2_add(a,b); }); },
                        +[]{ p2<glm::vec2>([](glm::vec2 a, glm::vec2 b){ return a+b; }); });
    add("P","vec2_sub", +[]{ p2<fv_vec2>([](fv_vec2 a, fv_vec2 b){ return fv_vec2_sub(a,b); }); },
                        +[]{ p2<glm::vec2>([](glm::vec2 a, glm::vec2 b){ return a-b; }); });
    add("P","vec2_neg", +[]{ p1<fv_vec2>([](fv_vec2 a){ return fv_vec2_neg(a); }); },
                        +[]{ p1<glm::vec2>([](glm::vec2 a){ return -a; }); });
    add("P","vec2_scale", +[]{ p1s<fv_vec2>([](fv_vec2 a, float s){ return fv_vec2_scale(a,s); }); },
                          +[]{ p1s<glm::vec2>([](glm::vec2 a, float s){ return a*s; }); });
    add("P","vec2_lerp", +[]{ p2s<fv_vec2>([](fv_vec2 a, fv_vec2 b, float t){ return fv_vec2_lerp(a,b,t); }); },
                         +[]{ p2s<glm::vec2>([](glm::vec2 a, glm::vec2 b, float t){ return glm::mix(a,b,t); }); });
    add("P","vec2_min", +[]{ p2<fv_vec2>([](fv_vec2 a, fv_vec2 b){ return fv_vec2_min(a,b); }); },
                        +[]{ p2<glm::vec2>([](glm::vec2 a, glm::vec2 b){ return glm::min(a,b); }); });
    add("P","vec2_max", +[]{ p2<fv_vec2>([](fv_vec2 a, fv_vec2 b){ return fv_vec2_max(a,b); }); },
                        +[]{ p2<glm::vec2>([](glm::vec2 a, glm::vec2 b){ return glm::max(a,b); }); });
    add("P","vec2_clamp", +[]{ p3<fv_vec2>([](fv_vec2 a, fv_vec2 b, fv_vec2 c){ return fv_vec2_clamp(a,b,c); }); },
                          +[]{ p3<glm::vec2>([](glm::vec2 a, glm::vec2 b, glm::vec2 c){ return glm::clamp(a,b,c); }); });
    add("P","vec2_madd", +[]{ p2s<fv_vec2>([](fv_vec2 a, fv_vec2 b, float s){ return fv_vec2_madd(a,b,s); }); },
                         +[]{ p2s<glm::vec2>([](glm::vec2 a, glm::vec2 b, float s){ return a*b+s; }); });
    add("P","vec2_addmul", +[]{ p2s<fv_vec2>([](fv_vec2 a, fv_vec2 b, float s){ return fv_vec2_addmul(a,b,s); }); },
                           +[]{ p2s<glm::vec2>([](glm::vec2 a, glm::vec2 b, float s){ return a+b*s; }); });
    add("P","vec2_maddv", +[]{ p3<fv_vec2>([](fv_vec2 a, fv_vec2 b, fv_vec2 c){ return fv_vec2_maddv(a,b,c); }); },
                          +[]{ p3<glm::vec2>([](glm::vec2 a, glm::vec2 b, glm::vec2 c){ return a*b+c; }); });
    add("P","vec2_dot", +[]{ s2<fv_vec2>([](fv_vec2 a, fv_vec2 b){ return fv_vec2_dot(a,b); }); },
                        +[]{ s2<glm::vec2>([](glm::vec2 a, glm::vec2 b){ return glm::dot(a,b); }); });
    add("P","vec2_len_sqr", +[]{ s1<fv_vec2>([](fv_vec2 a){ return fv_vec2_len_sqr(a); }); },
                           +[]{ s1<glm::vec2>([](glm::vec2 a){ return glm::length2(a); }); });
    add("P","vec2_len", +[]{ s1<fv_vec2>([](fv_vec2 a){ return fv_vec2_len(a); }); },
                        +[]{ s1<glm::vec2>([](glm::vec2 a){ return glm::length(a); }); });
    add("P","vec2_dist_sqr", +[]{ s2<fv_vec2>([](fv_vec2 a, fv_vec2 b){ return fv_vec2_dist_sqr(a,b); }); },
                           +[]{ s2<glm::vec2>([](glm::vec2 a, glm::vec2 b){ return glm::distance2(a,b); }); });
    add("P","vec2_dist", +[]{ s2<fv_vec2>([](fv_vec2 a, fv_vec2 b){ return fv_vec2_dist(a,b); }); },
                         +[]{ s2<glm::vec2>([](glm::vec2 a, glm::vec2 b){ return glm::distance(a,b); }); });
    add("P","vec2_normalize", +[]{ p1<fv_vec2>([](fv_vec2 a){ return fv_vec2_normalize(a); }); },
                             +[]{ p1<glm::vec2>([](glm::vec2 a){ return glm::normalize(a); }); });
    add("P","vec2_project", +[]{ p2<fv_vec2>([](fv_vec2 a, fv_vec2 b){ return fv_vec2_project(a,b); }); },
                            +[]{ p2<glm::vec2>([](glm::vec2 a, glm::vec2 b){ return glm::dot(a,b)/glm::length2(b)*b; }); });
    add("P","vec2_near", +[]{ si2<fv_vec2>([](fv_vec2 a, fv_vec2 b){ return fv_vec2_near(a,b,0.05f); }); },
                          +[]{ si2<glm::vec2>([](glm::vec2 a, glm::vec2 b){ return (int)glm::all(glm::epsilonEqual(a,b,0.05f)); }); });
    // ---- vec3 标量 ----
    add("P","vec3_add", +[]{ p2<fv_vec3>([](fv_vec3 a, fv_vec3 b){ return fv_vec3_add(a,b); }); },
                        +[]{ p2<glm::vec3>([](glm::vec3 a, glm::vec3 b){ return a+b; }); });
    add("P","vec3_sub", +[]{ p2<fv_vec3>([](fv_vec3 a, fv_vec3 b){ return fv_vec3_sub(a,b); }); },
                        +[]{ p2<glm::vec3>([](glm::vec3 a, glm::vec3 b){ return a-b; }); });
    add("P","vec3_neg", +[]{ p1<fv_vec3>([](fv_vec3 a){ return fv_vec3_neg(a); }); },
                        +[]{ p1<glm::vec3>([](glm::vec3 a){ return -a; }); });
    add("P","vec3_scale", +[]{ p1s<fv_vec3>([](fv_vec3 a, float s){ return fv_vec3_scale(a,s); }); },
                          +[]{ p1s<glm::vec3>([](glm::vec3 a, float s){ return a*s; }); });
    add("P","vec3_lerp", +[]{ p2s<fv_vec3>([](fv_vec3 a, fv_vec3 b, float t){ return fv_vec3_lerp(a,b,t); }); },
                         +[]{ p2s<glm::vec3>([](glm::vec3 a, glm::vec3 b, float t){ return glm::mix(a,b,t); }); });
    add("P","vec3_min", +[]{ p2<fv_vec3>([](fv_vec3 a, fv_vec3 b){ return fv_vec3_min(a,b); }); },
                        +[]{ p2<glm::vec3>([](glm::vec3 a, glm::vec3 b){ return glm::min(a,b); }); });
    add("P","vec3_max", +[]{ p2<fv_vec3>([](fv_vec3 a, fv_vec3 b){ return fv_vec3_max(a,b); }); },
                        +[]{ p2<glm::vec3>([](glm::vec3 a, glm::vec3 b){ return glm::max(a,b); }); });
    add("P","vec3_clamp", +[]{ p3<fv_vec3>([](fv_vec3 a, fv_vec3 b, fv_vec3 c){ return fv_vec3_clamp(a,b,c); }); },
                          +[]{ p3<glm::vec3>([](glm::vec3 a, glm::vec3 b, glm::vec3 c){ return glm::clamp(a,b,c); }); });
    add("P","vec3_madd", +[]{ p2s<fv_vec3>([](fv_vec3 a, fv_vec3 b, float s){ return fv_vec3_madd(a,b,s); }); },
                         +[]{ p2s<glm::vec3>([](glm::vec3 a, glm::vec3 b, float s){ return a*b+s; }); });
    add("P","vec3_addmul", +[]{ p2s<fv_vec3>([](fv_vec3 a, fv_vec3 b, float s){ return fv_vec3_addmul(a,b,s); }); },
                           +[]{ p2s<glm::vec3>([](glm::vec3 a, glm::vec3 b, float s){ return a+b*s; }); });
    add("P","vec3_cross", +[]{ p2<fv_vec3>([](fv_vec3 a, fv_vec3 b){ return fv_vec3_cross(a,b); }); },
                          +[]{ p2<glm::vec3>([](glm::vec3 a, glm::vec3 b){ return glm::cross(a,b); }); });
    add("P","vec3_dot", +[]{ s2<fv_vec3>([](fv_vec3 a, fv_vec3 b){ return fv_vec3_dot(a,b); }); },
                        +[]{ s2<glm::vec3>([](glm::vec3 a, glm::vec3 b){ return glm::dot(a,b); }); });
    add("P","vec3_len_sqr", +[]{ s1<fv_vec3>([](fv_vec3 a){ return fv_vec3_len_sqr(a); }); },
                           +[]{ s1<glm::vec3>([](glm::vec3 a){ return glm::length2(a); }); });
    add("P","vec3_len", +[]{ s1<fv_vec3>([](fv_vec3 a){ return fv_vec3_len(a); }); },
                        +[]{ s1<glm::vec3>([](glm::vec3 a){ return glm::length(a); }); });
    add("P","vec3_dist_sqr", +[]{ s2<fv_vec3>([](fv_vec3 a, fv_vec3 b){ return fv_vec3_dist_sqr(a,b); }); },
                           +[]{ s2<glm::vec3>([](glm::vec3 a, glm::vec3 b){ return glm::distance2(a,b); }); });
    add("P","vec3_dist", +[]{ s2<fv_vec3>([](fv_vec3 a, fv_vec3 b){ return fv_vec3_dist(a,b); }); },
                         +[]{ s2<glm::vec3>([](glm::vec3 a, glm::vec3 b){ return glm::distance(a,b); }); });
    add("P","vec3_normalize", +[]{ p1<fv_vec3>([](fv_vec3 a){ return fv_vec3_normalize(a); }); },
                             +[]{ p1<glm::vec3>([](glm::vec3 a){ return glm::normalize(a); }); });
    add("P","vec3_project", +[]{ p2<fv_vec3>([](fv_vec3 a, fv_vec3 b){ return fv_vec3_project(a,b); }); },
                            +[]{ p2<glm::vec3>([](glm::vec3 a, glm::vec3 b){ return glm::dot(a,b)/glm::length2(b)*b; }); });
    add("P","vec3_near", +[]{ si2<fv_vec3>([](fv_vec3 a, fv_vec3 b){ return fv_vec3_near(a,b,0.05f); }); },
                          +[]{ si2<glm::vec3>([](glm::vec3 a, glm::vec3 b){ return (int)glm::all(glm::epsilonEqual(a,b,0.05f)); }); });
    // ---- vec4 标量 ----
    add("P","vec4_add", +[]{ p2<fv_vec4>([](fv_vec4 a, fv_vec4 b){ return fv_vec4_add(a,b); }); },
                        +[]{ p2<glm::vec4>([](glm::vec4 a, glm::vec4 b){ return a+b; }); });
    add("P","vec4_sub", +[]{ p2<fv_vec4>([](fv_vec4 a, fv_vec4 b){ return fv_vec4_sub(a,b); }); },
                        +[]{ p2<glm::vec4>([](glm::vec4 a, glm::vec4 b){ return a-b; }); });
    add("P","vec4_neg", +[]{ p1<fv_vec4>([](fv_vec4 a){ return fv_vec4_neg(a); }); },
                        +[]{ p1<glm::vec4>([](glm::vec4 a){ return -a; }); });
    add("P","vec4_scale", +[]{ p1s<fv_vec4>([](fv_vec4 a, float s){ return fv_vec4_scale(a,s); }); },
                          +[]{ p1s<glm::vec4>([](glm::vec4 a, float s){ return a*s; }); });
    add("P","vec4_lerp", +[]{ p2s<fv_vec4>([](fv_vec4 a, fv_vec4 b, float t){ return fv_vec4_lerp(a,b,t); }); },
                         +[]{ p2s<glm::vec4>([](glm::vec4 a, glm::vec4 b, float t){ return glm::mix(a,b,t); }); });
    add("P","vec4_min", +[]{ p2<fv_vec4>([](fv_vec4 a, fv_vec4 b){ return fv_vec4_min(a,b); }); },
                        +[]{ p2<glm::vec4>([](glm::vec4 a, glm::vec4 b){ return glm::min(a,b); }); });
    add("P","vec4_max", +[]{ p2<fv_vec4>([](fv_vec4 a, fv_vec4 b){ return fv_vec4_max(a,b); }); },
                        +[]{ p2<glm::vec4>([](glm::vec4 a, glm::vec4 b){ return glm::max(a,b); }); });
    add("P","vec4_clamp", +[]{ p3<fv_vec4>([](fv_vec4 a, fv_vec4 b, fv_vec4 c){ return fv_vec4_clamp(a,b,c); }); },
                          +[]{ p3<glm::vec4>([](glm::vec4 a, glm::vec4 b, glm::vec4 c){ return glm::clamp(a,b,c); }); });
    add("P","vec4_madd", +[]{ p2s<fv_vec4>([](fv_vec4 a, fv_vec4 b, float s){ return fv_vec4_madd(a,b,s); }); },
                         +[]{ p2s<glm::vec4>([](glm::vec4 a, glm::vec4 b, float s){ return a*b+s; }); });
    add("P","vec4_addmul", +[]{ p2s<fv_vec4>([](fv_vec4 a, fv_vec4 b, float s){ return fv_vec4_addmul(a,b,s); }); },
                           +[]{ p2s<glm::vec4>([](glm::vec4 a, glm::vec4 b, float s){ return a+b*s; }); });
    add("P","vec4_dot", +[]{ s2<fv_vec4>([](fv_vec4 a, fv_vec4 b){ return fv_vec4_dot(a,b); }); },
                        +[]{ s2<glm::vec4>([](glm::vec4 a, glm::vec4 b){ return glm::dot(a,b); }); });
    add("P","vec4_len_sqr", +[]{ s1<fv_vec4>([](fv_vec4 a){ return fv_vec4_len_sqr(a); }); },
                           +[]{ s1<glm::vec4>([](glm::vec4 a){ return glm::length2(a); }); });
    add("P","vec4_len", +[]{ s1<fv_vec4>([](fv_vec4 a){ return fv_vec4_len(a); }); },
                        +[]{ s1<glm::vec4>([](glm::vec4 a){ return glm::length(a); }); });
    add("P","vec4_dist", +[]{ s2<fv_vec4>([](fv_vec4 a, fv_vec4 b){ return fv_vec4_dist(a,b); }); },
                         +[]{ s2<glm::vec4>([](glm::vec4 a, glm::vec4 b){ return glm::distance(a,b); }); });
    add("P","vec4_normalize", +[]{ p1<fv_vec4>([](fv_vec4 a){ return fv_vec4_normalize(a); }); },
                             +[]{ p1<glm::vec4>([](glm::vec4 a){ return glm::normalize(a); }); });
    // ---- aabb2 ----
    add("P","aabb2_merge", +[]{ p2<fv_aabb2>([](fv_aabb2 a, fv_aabb2 b){ fv_aabb2 r; r.lo=fv_vec2_min(a.lo,b.lo); r.hi=fv_vec2_max(a.hi,b.hi); return r; }); },
                         +[]{ p2<GABB>([](GABB a, GABB b){ GABB r; r.lo=glm::min(a.lo,b.lo); r.hi=glm::max(a.hi,b.hi); return r; }); });
    add("P","aabb_merge_point", +[]{ p2<fv_aabb2>([](fv_aabb2 box, fv_aabb2 p){ return fv_aabb2_merge(box, p.lo); }); },
                              +[]{ p2<GABB>([](GABB box, GABB p){ GABB r; r.lo=glm::min(box.lo,p.lo); r.hi=glm::max(box.hi,p.lo); return r; }); });
    add("P","aabb2_overlaps", +[]{ si2<fv_aabb2>([](fv_aabb2 a, fv_aabb2 b){ return fv_aabb2_overlaps(a,b); }); },
                             +[]{ si2<GABB>([](GABB a, GABB b){ return (int)(a.lo.x<=b.hi.x&&b.lo.x<=a.hi.x&&a.lo.y<=b.hi.y&&b.lo.y<=a.hi.y); }); });
    // ---- mat3 / mat4 标量 ----
    add("P","mat3_mul", +[]{ pm2<fv_mat3>([](fv_mat3 a, fv_mat3 b){ return fv_mat3_mul(&a,&b); }); },
                        +[]{ pm2<glm::mat3>([](glm::mat3 a, glm::mat3 b){ return a*b; }); });
    add("P","mat3_inverse", +[]{ pm1<fv_mat3>([](fv_mat3 a){ return fv_mat3_inverse(&a); }); },
                           +[]{ pm1<glm::mat3>([](glm::mat3 a){ return glm::inverse(a); }); });
    add("P","mat3_transpose", +[]{ pm1<fv_mat3>([](fv_mat3 a){ return fv_mat3_transpose(&a); }); },
                             +[]{ pm1<glm::mat3>([](glm::mat3 a){ return glm::transpose(a); }); });
    add("P","mat3_determinant", +[]{ sm1<fv_mat3>([](fv_mat3 a){ return fv_mat3_determinant(&a); }); },
                              +[]{ sm1<glm::mat3>([](glm::mat3 a){ return glm::determinant(a); }); });
    add("P","mat4_mul", +[]{ pm2<fv_mat4>([](fv_mat4 a, fv_mat4 b){ return fv_mat4_mul(&a,&b); }); },
                        +[]{ pm2<glm::mat4>([](glm::mat4 a, glm::mat4 b){ return a*b; }); });
    add("P","mat4_inverse", +[]{ pm1<fv_mat4>([](fv_mat4 a){ return fv_mat4_inverse(&a); }); },
                           +[]{ pm1<glm::mat4>([](glm::mat4 a){ return glm::inverse(a); }); });
    // ---- 常量算子 + vec 数组 ----
    add("P","mat3_transform_vec3", +[]{ c1<fv_vec3>([](fv_vec3 v){ return fv_mat3_transform_vec3(&gFvM3,v); }); },
                                 +[]{ c1<glm::vec3>([](glm::vec3 v){ return gGlM3*v; }); });
    add("P","quat_rotate_vec3", +[]{ c1<fv_vec3>([](fv_vec3 v){ return fv_quat_rotate_vec3(gFvQ,v); }); },
                              +[]{ c1<glm::vec3>([](glm::vec3 v){ return glm::rotate(gGlQ,v); }); });
    // ---- quat 标量 ----
    add("P","quat_mul", +[]{ q2<fv_quat>([](fv_quat a, fv_quat b){ return fv_quat_mul(a,b); }); },
                        +[]{ gq2([](glm::quat a, glm::quat b){ return a*b; }); });
    add("P","quat_add", +[]{ q2<fv_quat>([](fv_quat a, fv_quat b){ return fv_quat_add(a,b); }); },
                        +[]{ gq2([](glm::quat a, glm::quat b){ return a+b; }); });
    add("P","quat_sub", +[]{ q2<fv_quat>([](fv_quat a, fv_quat b){ return fv_quat_sub(a,b); }); },
                        +[]{ gq2([](glm::quat a, glm::quat b){ return a-b; }); });
    add("P","quat_neg", +[]{ q1<fv_quat>([](fv_quat a){ return fv_quat_neg(a); }); },
                        +[]{ gq1([](glm::quat a){ return -a; }); });
    add("P","quat_scale", +[]{ q1<fv_quat>([](fv_quat a){ return fv_quat_scale(a,1.5f); }); },
                          +[]{ gq1([](glm::quat a){ return a*1.5f; }); });
    add("P","quat_conjugate", +[]{ q1<fv_quat>([](fv_quat a){ return fv_quat_conjugate(a); }); },
                            +[]{ gq1([](glm::quat a){ return glm::conjugate(a); }); });
    add("P","quat_inverse", +[]{ q1<fv_quat>([](fv_quat a){ return fv_quat_inverse(a); }); },
                          +[]{ gq1([](glm::quat a){ return glm::inverse(a); }); });
    add("P","quat_normalize", +[]{ q1<fv_quat>([](fv_quat a){ return fv_quat_normalize(a); }); },
                            +[]{ gq1([](glm::quat a){ return glm::normalize(a); }); });
    add("P","quat_len", +[]{ qs1<fv_quat>([](fv_quat a){ return fv_quat_len(a); }); },
                        +[]{ gqs1([](glm::quat a){ return glm::length(a); }); });
    add("P","quat_dot", +[]{ qs2<fv_quat>([](fv_quat a, fv_quat b){ return fv_quat_dot(a,b); }); },
                        +[]{ gqs2([](glm::quat a, glm::quat b){ return glm::dot(a,b); }); });
    add("P","quat_nlerp_vs_fastMix", +[]{ q2s<fv_quat>([](fv_quat a, fv_quat b, float t){ return fv_quat_nlerp(a,b,t); }); },
                          +[]{ gq2s([](glm::quat a, glm::quat b, float t){ return glm::fastMix(a,b,t); }); });
    add("P","quat_slerp", +[]{ q2s<fv_quat>([](fv_quat a, fv_quat b, float t){ return fv_quat_slerp(a,b,t); }); },
                          +[]{ gq2s([](glm::quat a, glm::quat b, float t){ return glm::slerp(a,b,t); }); });
    add("P","mat3_from_quat", +[]{ xq1<fv_quat,fv_mat3>([](fv_quat q){ return fv_mat3_from_quat(q); }); },
                            +[]{ xq1<glm::quat,glm::mat3>([](glm::quat q){ return glm::mat3_cast(q); }); });
    add("P","quat_from_mat3", +[]{ xq1<fv_mat3,fv_quat>([](fv_mat3 m){ return fv_quat_from_mat3(&m); }); },
                            +[]{ xq1<glm::mat3,glm::quat>([](glm::mat3 m){ return glm::quat_cast(m); }); });
    return v;
}

// ------------------------- B 组：batch / SoA / stream -------------------------
#define BA(nm,FVCALL,...) static void b_##nm(){ FVCALL; gSink+=gO[0]; } static void bg_##nm(){ __VA_ARGS__; gSink+=gO[0]; }
BA(v2_add_x4,   fv_vec2_add_x4((fv_vec2*)gA,(fv_vec2*)gB,(fv_vec2*)gO,gN),      { const glm::vec2*a=(glm::vec2*)gA,*b=(glm::vec2*)gB; glm::vec2*o=(glm::vec2*)gO; for(size_t i=0;i<gN;i++) o[i]=a[i]+b[i]; })
BA(v2_sub_x4,   fv_vec2_sub_x4((fv_vec2*)gA,(fv_vec2*)gB,(fv_vec2*)gO,gN),      { const glm::vec2*a=(glm::vec2*)gA,*b=(glm::vec2*)gB; glm::vec2*o=(glm::vec2*)gO; for(size_t i=0;i<gN;i++) o[i]=a[i]-b[i]; })
BA(v2_scale_x4, fv_vec2_scale_x4((fv_vec2*)gA,1.5f,(fv_vec2*)gO,gN),            { const glm::vec2*a=(glm::vec2*)gA; glm::vec2*o=(glm::vec2*)gO; for(size_t i=0;i<gN;i++) o[i]=a[i]*1.5f; })
BA(v2_lerp_x4,  fv_vec2_lerp_x4((fv_vec2*)gA,(fv_vec2*)gB,0.25f,(fv_vec2*)gO,gN), { const glm::vec2*a=(glm::vec2*)gA,*b=(glm::vec2*)gB; glm::vec2*o=(glm::vec2*)gO; for(size_t i=0;i<gN;i++) o[i]=glm::mix(a[i],b[i],0.25f); })
BA(v2_madd_x4,  fv_vec2_madd_x4((fv_vec2*)gA,(fv_vec2*)gB,1.5f,(fv_vec2*)gO,gN), { const glm::vec2*a=(glm::vec2*)gA,*b=(glm::vec2*)gB; glm::vec2*o=(glm::vec2*)gO; for(size_t i=0;i<gN;i++) o[i]=a[i]*b[i]+1.5f; })
BA(v2_maddv_x4, fv_vec2_maddv_x4((fv_vec2*)gA,(fv_vec2*)gB,(fv_vec2*)gC,(fv_vec2*)gO,gN), { const glm::vec2*a=(glm::vec2*)gA,*b=(glm::vec2*)gB,*c=(glm::vec2*)gC; glm::vec2*o=(glm::vec2*)gO; for(size_t i=0;i<gN;i++) o[i]=a[i]*b[i]+c[i]; })
BA(v2_addmul_x4,fv_vec2_addmul_x4((fv_vec2*)gA,(fv_vec2*)gB,1.5f,gN),           { glm::vec2*x=(glm::vec2*)gA; const glm::vec2*v=(glm::vec2*)gB; for(size_t i=0;i<gN;i++) x[i]=x[i]+v[i]*1.5f; })
BA(v2_add_x4_mask, fv_vec2_add_x4_mask((fv_vec2*)gA,(fv_vec2*)gB,(fv_vec2*)gO,gMask,gN), { const glm::vec2*a=(glm::vec2*)gA,*b=(glm::vec2*)gB; glm::vec2*o=(glm::vec2*)gO; for(size_t i=0;i<gN;i++) o[i]=gMask[i]?a[i]+b[i]:a[i]; })
BA(v2_add_x4_idx,  fv_vec2_add_x4_idx((fv_vec2*)gA,(fv_vec2*)gB,(fv_vec2*)gO,gIdx,gN),   { const glm::vec2*a=(glm::vec2*)gA,*b=(glm::vec2*)gB; glm::vec2*o=(glm::vec2*)gO; for(size_t k=0;k<gN;k++){ uint32_t i=gIdx[k]; o[i]=a[i]+b[i]; } })
BA(v4_add_x2,   fv_vec4_add_x2((fv_vec4*)gA,(fv_vec4*)gB,(fv_vec4*)gO,gN),      { const glm::vec4*a=(glm::vec4*)gA,*b=(glm::vec4*)gB; glm::vec4*o=(glm::vec4*)gO; for(size_t i=0;i<gN;i++) o[i]=a[i]+b[i]; })
BA(v4_sub_x2,   fv_vec4_sub_x2((fv_vec4*)gA,(fv_vec4*)gB,(fv_vec4*)gO,gN),      { const glm::vec4*a=(glm::vec4*)gA,*b=(glm::vec4*)gB; glm::vec4*o=(glm::vec4*)gO; for(size_t i=0;i<gN;i++) o[i]=a[i]-b[i]; })
BA(v4_scale_x2, fv_vec4_scale_x2((fv_vec4*)gA,1.5f,(fv_vec4*)gO,gN),            { const glm::vec4*a=(glm::vec4*)gA; glm::vec4*o=(glm::vec4*)gO; for(size_t i=0;i<gN;i++) o[i]=a[i]*1.5f; })
BA(v4_lerp_x2,  fv_vec4_lerp_x2((fv_vec4*)gA,(fv_vec4*)gB,0.25f,(fv_vec4*)gO,gN), { const glm::vec4*a=(glm::vec4*)gA,*b=(glm::vec4*)gB; glm::vec4*o=(glm::vec4*)gO; for(size_t i=0;i<gN;i++) o[i]=glm::mix(a[i],b[i],0.25f); })
BA(v4_madd_x2,  fv_vec4_madd_x2((fv_vec4*)gA,(fv_vec4*)gB,1.5f,(fv_vec4*)gO,gN), { const glm::vec4*a=(glm::vec4*)gA,*b=(glm::vec4*)gB; glm::vec4*o=(glm::vec4*)gO; for(size_t i=0;i<gN;i++) o[i]=a[i]*b[i]+1.5f; })
BA(v4_addmul_x2,fv_vec4_addmul_x2((fv_vec4*)gA,(fv_vec4*)gB,1.5f,gN),           { glm::vec4*x=(glm::vec4*)gA; const glm::vec4*v=(glm::vec4*)gB; for(size_t i=0;i<gN;i++) x[i]=x[i]+v[i]*1.5f; })
BA(m4_transform_points, fv_mat4_transform_points(&gFvM4,(fv_vec4*)gA,(fv_vec4*)gO,gN), { const glm::vec4*a=(glm::vec4*)gA; glm::vec4*o=(glm::vec4*)gO; glm::mat4 m=gGlM4; for(size_t i=0;i<gN;i++) o[i]=m*a[i]; })
BA(m4_transform_dirs,   fv_mat4_transform_dirs(&gFvM4,(fv_vec4*)gA,(fv_vec4*)gO,gN),   { const glm::vec4*a=(glm::vec4*)gA; glm::vec4*o=(glm::vec4*)gO; glm::mat4 m=gGlM4; for(size_t i=0;i<gN;i++){ glm::vec4 r=m*a[i]; o[i]=glm::vec4(r); o[i].w=0.f; } })
// _stream（non-temporal store）：对手是同函数族的普通版，GLM 侧无对应物
BA(v2_add_x4_stream,  fv_vec2_add_x4_stream((fv_vec2*)gA,(fv_vec2*)gB,(fv_vec2*)gO,gN),      fv_vec2_add_x4((fv_vec2*)gA,(fv_vec2*)gB,(fv_vec2*)gO,gN))
BA(v4_add_x2_stream,  fv_vec4_add_x2_stream((fv_vec4*)gA,(fv_vec4*)gB,(fv_vec4*)gO,gN),      fv_vec4_add_x2((fv_vec4*)gA,(fv_vec4*)gB,(fv_vec4*)gO,gN))
BA(m4_points_stream,  fv_mat4_transform_points_stream(&gFvM4,(fv_vec4*)gA,(fv_vec4*)gO,gN),  fv_mat4_transform_points(&gFvM4,(fv_vec4*)gA,(fv_vec4*)gO,gN))
// SoA
#define BS(nm,FVCALL,...) static void s_##nm(){ FVCALL; gSink+=SA[0][0]; } static void sg_##nm(){ __VA_ARGS__; gSink+=SA[0][0]; }
BS(v2_add_soa, fv_vec2_add_soa(SA[0],SA[1],SB[0],SB[1],SO[0],SO[1],gN),
   { for(size_t i=0;i<gN;i++){ glm::vec2 r = glm::vec2(SA[0][i],SA[1][i])+glm::vec2(SB[0][i],SB[1][i]); SO[0][i]=r.x; SO[1][i]=r.y; } })
BS(v2_scale_soa, fv_vec2_scale_soa(SA[0],SA[1],1.5f,SO[0],SO[1],gN),
   { for(size_t i=0;i<gN;i++){ glm::vec2 r = glm::vec2(SA[0][i],SA[1][i])*1.5f; SO[0][i]=r.x; SO[1][i]=r.y; } })
BS(v2_len_soa, fv_vec2_len_soa(SA[0],SA[1],SO[0],gN),
   { for(size_t i=0;i<gN;i++) SO[0][i]=glm::length(glm::vec2(SA[0][i],SA[1][i])); })
BS(v2_normalize_soa, fv_vec2_normalize_soa(SA[0],SA[1],SO[0],SO[1],gN),
   { for(size_t i=0;i<gN;i++){ glm::vec2 r=glm::normalize(glm::vec2(SA[0][i],SA[1][i])); SO[0][i]=r.x; SO[1][i]=r.y; } })
BS(v2_dot_soa, gSink+=fv_vec2_dot_soa(SA[0],SA[1],SB[0],SB[1],gN),
   { float s=0; for(size_t i=0;i<gN;i++) s+=glm::dot(glm::vec2(SA[0][i],SA[1][i]),glm::vec2(SB[0][i],SB[1][i])); gSink+=s; })
BS(v2_madd_normalize_soa, fv_vec2_madd_normalize_soa(SA[0],SA[1],SB[0],SB[1],1.5f,SO[0],SO[1],gN),
   { for(size_t i=0;i<gN;i++){ glm::vec2 a(SA[0][i],SA[1][i]),b(SB[0][i],SB[1][i]); glm::vec2 r=glm::normalize(a*b+glm::vec2(1.5f)); SO[0][i]=r.x; SO[1][i]=r.y; } })
BS(v3_add_soa, fv_vec3_add_soa(SA[0],SA[1],SA[2],SB[0],SB[1],SB[2],SO[0],SO[1],SO[2],gN),
   { for(size_t i=0;i<gN;i++){ glm::vec3 r=glm::vec3(SA[0][i],SA[1][i],SA[2][i])+glm::vec3(SB[0][i],SB[1][i],SB[2][i]); SO[0][i]=r.x; SO[1][i]=r.y; SO[2][i]=r.z; } })
BS(v3_scale_soa, fv_vec3_scale_soa(SA[0],SA[1],SA[2],1.5f,SO[0],SO[1],SO[2],gN),
   { for(size_t i=0;i<gN;i++){ glm::vec3 r=glm::vec3(SA[0][i],SA[1][i],SA[2][i])*1.5f; SO[0][i]=r.x; SO[1][i]=r.y; SO[2][i]=r.z; } })
BS(v3_len_soa, fv_vec3_len_soa(SA[0],SA[1],SA[2],SO[0],gN),
   { for(size_t i=0;i<gN;i++) SO[0][i]=glm::length(glm::vec3(SA[0][i],SA[1][i],SA[2][i])); })
BS(v3_normalize_soa, fv_vec3_normalize_soa(SA[0],SA[1],SA[2],SO[0],SO[1],SO[2],gN),
   { for(size_t i=0;i<gN;i++){ glm::vec3 r=glm::normalize(glm::vec3(SA[0][i],SA[1][i],SA[2][i])); SO[0][i]=r.x; SO[1][i]=r.y; SO[2][i]=r.z; } })
BS(v3_dot_soa, gSink+=fv_vec3_dot_soa(SA[0],SA[1],SA[2],SB[0],SB[1],SB[2],gN),
   { float s=0; for(size_t i=0;i<gN;i++) s+=glm::dot(glm::vec3(SA[0][i],SA[1][i],SA[2][i]),glm::vec3(SB[0][i],SB[1][i],SB[2][i])); gSink+=s; })
BS(v4_len_soa, fv_vec4_len_soa(SA[0],SA[1],SA[2],SA[3],SO[0],gN),
   { for(size_t i=0;i<gN;i++) SO[0][i]=glm::length(glm::vec4(SA[0][i],SA[1][i],SA[2][i],SA[3][i])); })
BS(v4_normalize_soa, fv_vec4_normalize_soa(SA[0],SA[1],SA[2],SA[3],SO[0],SO[1],SO[2],SO[3],gN),
   { for(size_t i=0;i<gN;i++){ glm::vec4 r=glm::normalize(glm::vec4(SA[0][i],SA[1][i],SA[2][i],SA[3][i])); SO[0][i]=r.x; SO[1][i]=r.y; SO[2][i]=r.z; SO[3][i]=r.w; } })
BS(v4_dot_soa, gSink+=fv_vec4_dot_soa(SA[0],SA[1],SA[2],SA[3],SB[0],SB[1],SB[2],SB[3],gN),
   { float s=0; for(size_t i=0;i<gN;i++) s+=glm::dot(glm::vec4(SA[0][i],SA[1][i],SA[2][i],SA[3][i]),glm::vec4(SB[0][i],SB[1][i],SB[2][i],SB[3][i])); gSink+=s; })
BS(v2_normalize_soa_fast, fv_vec2_normalize_soa_fast(SA[0],SA[1],SO[0],SO[1],gN),
   fv_vec2_normalize_soa(SA[0],SA[1],SO[0],SO[1],gN))
BS(v3_normalize_soa_fast, fv_vec3_normalize_soa_fast(SA[0],SA[1],SA[2],SO[0],SO[1],SO[2],gN),
   fv_vec3_normalize_soa(SA[0],SA[1],SA[2],SO[0],SO[1],SO[2],gN))
BS(v4_normalize_soa_fast, fv_vec4_normalize_soa_fast(SA[0],SA[1],SA[2],SA[3],SO[0],SO[1],SO[2],SO[3],gN),
   fv_vec4_normalize_soa(SA[0],SA[1],SA[2],SA[3],SO[0],SO[1],SO[2],SO[3],gN))
BS(q_mul_soa, fv_quat_mul_soa(SA[0],SA[1],SA[2],SA[3],SB[0],SB[1],SB[2],SB[3],SO[0],SO[1],SO[2],SO[3],gN),
   { for(size_t i=0;i<gN;i++){ glm::quat a(SA[3][i],SA[0][i],SA[1][i],SA[2][i]),b(SB[3][i],SB[0][i],SB[1][i],SB[2][i]),r=a*b; SO[0][i]=r.x; SO[1][i]=r.y; SO[2][i]=r.z; SO[3][i]=r.w; } })
BS(q_normalize_soa, fv_quat_normalize_soa(SA[0],SA[1],SA[2],SA[3],SO[0],SO[1],SO[2],SO[3],gN),
   { for(size_t i=0;i<gN;i++){ glm::quat r=glm::normalize(glm::quat(SA[3][i],SA[0][i],SA[1][i],SA[2][i])); SO[0][i]=r.x; SO[1][i]=r.y; SO[2][i]=r.z; SO[3][i]=r.w; } })
BS(q_rotate_soa, fv_quat_rotate_vec3_soa(SA[0],SA[1],SA[2],SA[3],SB[0],SB[1],SB[2],SO[0],SO[1],SO[2],gN),
   { for(size_t i=0;i<gN;i++){ glm::vec3 r=glm::rotate(glm::quat(SA[3][i],SA[0][i],SA[1][i],SA[2][i]),glm::vec3(SB[0][i],SB[1][i],SB[2][i])); SO[0][i]=r.x; SO[1][i]=r.y; SO[2][i]=r.z; } })
BS(m3_mul_soa, fv_mat3_mul_soa(MA[0],MA[1],MA[2],MA[3],MA[4],MA[5],MA[6],MA[7],MA[8],MB[0],MB[1],MB[2],MB[3],MB[4],MB[5],MB[6],MB[7],MB[8],MO[0],MO[1],MO[2],MO[3],MO[4],MO[5],MO[6],MO[7],MO[8],gN),
   { for(size_t i=0;i<gN;i++){ glm::mat3 A(MA[0][i],MA[1][i],MA[2][i],MA[3][i],MA[4][i],MA[5][i],MA[6][i],MA[7][i],MA[8][i]);
       glm::mat3 B(MB[0][i],MB[1][i],MB[2][i],MB[3][i],MB[4][i],MB[5][i],MB[6][i],MB[7][i],MB[8][i]); glm::mat3 R=A*B;
       MO[0][i]=R[0][0]; MO[1][i]=R[0][1]; MO[2][i]=R[0][2]; MO[3][i]=R[1][0]; MO[4][i]=R[1][1]; MO[5][i]=R[1][2]; MO[6][i]=R[2][0]; MO[7][i]=R[2][1]; MO[8][i]=R[2][2]; } })
BS(m3_inverse_soa, fv_mat3_inverse_soa(MA[0],MA[1],MA[2],MA[3],MA[4],MA[5],MA[6],MA[7],MA[8],MO[0],MO[1],MO[2],MO[3],MO[4],MO[5],MO[6],MO[7],MO[8],gN),
   { for(size_t i=0;i<gN;i++){ glm::mat3 A(MA[0][i],MA[1][i],MA[2][i],MA[3][i],MA[4][i],MA[5][i],MA[6][i],MA[7][i],MA[8][i]); glm::mat3 R=glm::inverse(A);
       MO[0][i]=R[0][0]; MO[1][i]=R[0][1]; MO[2][i]=R[0][2]; MO[3][i]=R[1][0]; MO[4][i]=R[1][1]; MO[5][i]=R[1][2]; MO[6][i]=R[2][0]; MO[7][i]=R[2][1]; MO[8][i]=R[2][2]; } })
BS(m3_transform_soa, fv_mat3_transform_vec3_soa(MA[0],MA[1],MA[2],MA[3],MA[4],MA[5],MA[6],MA[7],MA[8],SB[0],SB[1],SB[2],SO[0],SO[1],SO[2],gN),
   { for(size_t i=0;i<gN;i++){ glm::mat3 A(MA[0][i],MA[1][i],MA[2][i],MA[3][i],MA[4][i],MA[5][i],MA[6][i],MA[7][i],MA[8][i]);
       glm::vec3 r=A*glm::vec3(SB[0][i],SB[1][i],SB[2][i]); SO[0][i]=r.x; SO[1][i]=r.y; SO[2][i]=r.z; } })
// S 组：batch 入口在小 n 下的调用开销（对手 = GLM 同数据内联循环）
static void small1_fv(){ for(size_t i=0;i<gN;i++) fv_vec4_add_x2((fv_vec4*)gA+i,(fv_vec4*)gB+i,(fv_vec4*)gO+i,1); gSink+=gO[0]; }
static void small2_fv(){ for(size_t i=0;i+1<gN;i+=2) fv_vec4_add_x2((fv_vec4*)gA+i,(fv_vec4*)gB+i,(fv_vec4*)gO+i,2); gSink+=gO[0]; }
static void small4_fv(){ for(size_t i=0;i+3<gN;i+=4) fv_vec2_add_x4((fv_vec2*)gA+i,(fv_vec2*)gB+i,(fv_vec2*)gO+i,4); gSink+=gO[0]; }
static void small_g(){ const glm::vec4*a=(glm::vec4*)gA,*b=(glm::vec4*)gB; glm::vec4*o=(glm::vec4*)gO; for(size_t i=0;i<gN;i++) o[i]=a[i]+b[i]; gSink+=gO[0]; }

static std::vector<Case> build_b(){
    std::vector<Case> v;
    auto add=[&](const char*n,Fn a,Fn b){ v.push_back({"B",n,a,b}); };
    add("v2_add_x4",b_v2_add_x4,bg_v2_add_x4);            add("v2_sub_x4",b_v2_sub_x4,bg_v2_sub_x4);
    add("v2_scale_x4",b_v2_scale_x4,bg_v2_scale_x4);      add("v2_lerp_x4",b_v2_lerp_x4,bg_v2_lerp_x4);
    add("v2_madd_x4",b_v2_madd_x4,bg_v2_madd_x4);         add("v2_maddv_x4",b_v2_maddv_x4,bg_v2_maddv_x4);
    add("v2_addmul_x4",b_v2_addmul_x4,bg_v2_addmul_x4);   add("v2_add_x4_mask",b_v2_add_x4_mask,bg_v2_add_x4_mask);
    add("v2_add_x4_idx",b_v2_add_x4_idx,bg_v2_add_x4_idx);
    add("v4_add_x2",b_v4_add_x2,bg_v4_add_x2);            add("v4_sub_x2",b_v4_sub_x2,bg_v4_sub_x2);
    add("v4_scale_x2",b_v4_scale_x2,bg_v4_scale_x2);      add("v4_lerp_x2",b_v4_lerp_x2,bg_v4_lerp_x2);
    add("v4_madd_x2",b_v4_madd_x2,bg_v4_madd_x2);         add("v4_addmul_x2",b_v4_addmul_x2,bg_v4_addmul_x2);
    add("m4_transform_points",b_m4_transform_points,bg_m4_transform_points);
    add("m4_transform_dirs",b_m4_transform_dirs,bg_m4_transform_dirs);
    add("v2_add_x4_stream_NTvs普通",b_v2_add_x4_stream,bg_v2_add_x4_stream);
    add("v4_add_x2_stream_NTvs普通",b_v4_add_x2_stream,bg_v4_add_x2_stream);
    add("m4_points_stream_NTvs普通",b_m4_points_stream,bg_m4_points_stream);
    add("v2_add_soa",s_v2_add_soa,sg_v2_add_soa);                 add("v2_scale_soa",s_v2_scale_soa,sg_v2_scale_soa);
    add("v2_len_soa",s_v2_len_soa,sg_v2_len_soa);                 add("v2_normalize_soa",s_v2_normalize_soa,sg_v2_normalize_soa);
    add("v2_dot_soa",s_v2_dot_soa,sg_v2_dot_soa);                 add("v2_madd_normalize_soa",s_v2_madd_normalize_soa,sg_v2_madd_normalize_soa);
    add("v3_add_soa",s_v3_add_soa,sg_v3_add_soa);                 add("v3_scale_soa",s_v3_scale_soa,sg_v3_scale_soa);
    add("v3_len_soa",s_v3_len_soa,sg_v3_len_soa);                 add("v3_normalize_soa",s_v3_normalize_soa,sg_v3_normalize_soa);
    add("v3_dot_soa",s_v3_dot_soa,sg_v3_dot_soa);
    add("v4_len_soa",s_v4_len_soa,sg_v4_len_soa);                 add("v4_normalize_soa",s_v4_normalize_soa,sg_v4_normalize_soa);
    add("v4_dot_soa",s_v4_dot_soa,sg_v4_dot_soa);
    add("v2_normalize_fast_vs精确",s_v2_normalize_soa_fast,sg_v2_normalize_soa_fast);
    add("v3_normalize_fast_vs精确",s_v3_normalize_soa_fast,sg_v3_normalize_soa_fast);
    add("v4_normalize_fast_vs精确",s_v4_normalize_soa_fast,sg_v4_normalize_soa_fast);
    add("q_mul_soa",s_q_mul_soa,sg_q_mul_soa);                   add("q_normalize_soa",s_q_normalize_soa,sg_q_normalize_soa);
    add("q_rotate_soa",s_q_rotate_soa,sg_q_rotate_soa);
    add("m3_mul_soa",s_m3_mul_soa,sg_m3_mul_soa);                add("m3_inverse_soa",s_m3_inverse_soa,sg_m3_inverse_soa);
    add("m3_transform_soa",s_m3_transform_soa,sg_m3_transform_soa);
    v.push_back({"S","batch_n1_vs_glm_loop",small1_fv,small_g});
    v.push_back({"S","batch_n2_vs_glm_loop",small2_fv,small_g});
    v.push_back({"S","batch_n4_vs_glm_loop",small4_fv,small_g});
    v.push_back({"AA","add4_same_impl",+[]{ aa_add4_a((fv_vec4*)gA,(fv_vec4*)gB,(fv_vec4*)gO,gN); },
                                    +[]{ aa_add4_b((fv_vec4*)gA,(fv_vec4*)gB,(fv_vec4*)gO,gN); }});
    return v;
}

// ------------------------- main -------------------------
int main(int argc, char **argv){
    size_t N = argc>1 ? (size_t)strtoul(argv[1],nullptr,10) : 8192u;
    int reps = argc>2 ? atoi(argv[2]) : 21;
    const char *only = argc>3 ? argv[3] : "";
    gN = N & ~(size_t)1; gM = gN/4;
    size_t fN = gN*4;
    gA=(float*)a64(fN*4); gB=(float*)a64(fN*4); gC=(float*)a64(fN*4); gO=(float*)a64(fN*4);
    for(int c=0;c<4;c++){ SA[c]=(float*)a64(gN*4); SB[c]=(float*)a64(gN*4); SO[c]=(float*)a64(gN*4); }
    for(int c=0;c<9;c++){ MA[c]=(float*)a64(gN*4); MB[c]=(float*)a64(gN*4); MO[c]=(float*)a64(gN*4); }
    gMask=(uint8_t*)a64(gN); gIdx=(uint32_t*)a64(gN*4); gAng=(float*)a64(gN*4);
    gQA=(fv_quat*)a64(gN*16); gQB=(fv_quat*)a64(gN*16);
    GQA=(glm::quat*)a64(gN*sizeof(glm::quat)); GQB=(glm::quat*)a64(gN*sizeof(glm::quat));
    fill(gA,fN); fill(gB,fN); fill(gC,fN);
    for(int c=0;c<4;c++){ fill(SA[c],gN); fill(SB[c],gN); }
    for(int c=0;c<9;c++){ fill(MA[c],gN); fill(MB[c],gN); }
    unsigned s=24680u;
    for(size_t i=0;i<gN;i++){ gMask[i]=(uint8_t)((i*7+3)%5!=0); gIdx[i]=(uint32_t)((i*2654435761u)%(gN?gN:1)); gAng[i]=frand(&s); }
    for(size_t i=0;i<gN;i++){                          // quat 数据：同一逻辑值，两种布局各自填
        float x=frand(&s),y=frand(&s),z=frand(&s),w=0.3f+frand(&s);
        float n=1.0f/sqrtf(x*x+y*y+z*z+w*w); x*=n;y*=n;z*=n;w*=n;
        gQA[i]=(fv_quat){x,y,z,w}; gQB[i]=(fv_quat){z,x,w,y};
        GQA[i]=glm::quat(w,x,y,z); GQB[i]=glm::quat(y,z,x,w);
    }
    // 常量算子：fv 行主序，GLM 列主序（转置映射），mat4 取仿射形式
    float Am4[16],Am3[9];
    for(int r=0;r<4;r++) for(int c=0;c<4;c++) Am4[r*4+c] = (r==3)?(float)(c==3):0.1f+0.05f*(float)((r*4+c)%17);
    for(int r=0;r<3;r++) for(int c=0;c<3;c++) Am3[r*3+c] = 0.1f+0.07f*(float)((r*3+c)%11);
    memcpy(&gFvM4,Am4,64); memcpy(&gFvM3,Am3,36);
    for(int r=0;r<4;r++) for(int c=0;c<4;c++) glm::value_ptr(gGlM4)[c*4+r]=Am4[r*4+c];
    for(int r=0;r<3;r++) for(int c=0;c<3;c++) glm::value_ptr(gGlM3)[c*3+r]=Am3[r*3+c];
    gFvQ=(fv_quat){0.2f,0.3f,0.4f,0.84f}; gGlQ=glm::quat(0.84f,0.2f,0.3f,0.4f);
    verify();
    printf("# flashvec=A GLM=B; ratio=t_glm/t_fv, >1 = flashvec 快。stream/fast 组对手是 flashvec 自己的精确版。\n");
    printf("# build=%s N=%zu reps=%d\n", BUILD_TAG, N, reps);
    auto cs=build_cases(); for(auto &c:cs) if(!*only||!strcmp(c.grp,only)) run(c,reps);
    auto bs=build_b();     for(auto &c:bs) if(!*only||!strcmp(c.grp,only)) run(c,reps);
    fprintf(stderr,"[sink] %f\n",(double)gSink);
    return 0;
}
