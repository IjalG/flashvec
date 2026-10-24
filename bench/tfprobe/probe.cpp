
#include <immintrin.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <vector>

struct V4 { float x, y, z, w; };
struct Mat { float m[12]; };   /* 3 rows x 4 cols, affine, row-major */

static double now_ns() { timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec*1e9 + t.tv_nsec; }

/* ---------------- A: plain loop, compiler auto-vectorizes ---------------- */
__attribute__((noinline)) void tf_aos_plain(const Mat& M, const V4* __restrict in, V4* __restrict out, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        float x = in[i].x, y = in[i].y, z = in[i].z;
        out[i].x = M.m[0]*x + M.m[1]*y + M.m[2]*z + M.m[3];
        out[i].y = M.m[4]*x + M.m[5]*y + M.m[6]*z + M.m[7];
        out[i].z = M.m[8]*x + M.m[9]*y + M.m[10]*z + M.m[11];
        out[i].w = 1.0f;
    }
}
__attribute__((noinline)) void tf_aos_plain2(const Mat& M, const V4* __restrict in, V4* __restrict out, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        float x = in[i].x, y = in[i].y, z = in[i].z;
        out[i].x = M.m[0]*x + M.m[1]*y + M.m[2]*z + M.m[3];
        out[i].y = M.m[4]*x + M.m[5]*y + M.m[6]*z + M.m[7];
        out[i].z = M.m[8]*x + M.m[9]*y + M.m[10]*z + M.m[11];
        out[i].w = 1.0f;
    }
}

/* ---------------- B: hand-written AVX2, AoS, 2 points per ymm ---------------- */
__attribute__((target("avx2"))) __attribute__((noinline))
void tf_aos_avx2(const Mat& M, const V4* __restrict in, V4* __restrict out, size_t n) {
    const __m256 c0 = _mm256_setr_ps(M.m[0],M.m[4],M.m[8],0.f, M.m[0],M.m[4],M.m[8],0.f);
    const __m256 c1 = _mm256_setr_ps(M.m[1],M.m[5],M.m[9],0.f, M.m[1],M.m[5],M.m[9],0.f);
    const __m256 c2 = _mm256_setr_ps(M.m[2],M.m[6],M.m[10],0.f,M.m[2],M.m[6],M.m[10],0.f);
    const __m256 c3 = _mm256_setr_ps(M.m[3],M.m[7],M.m[11],0.f,M.m[3],M.m[7],M.m[11],0.f);
    const __m256 wone = _mm256_setr_ps(0,0,0,1.f, 0,0,0,1.f);
    const __m256i ix = _mm256_setr_epi32(0,0,0,0,4,4,4,4);
    const __m256i iy = _mm256_setr_epi32(1,1,1,1,5,5,5,5);
    const __m256i iz = _mm256_setr_epi32(2,2,2,2,6,6,6,6);
    size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        const __m256 p = _mm256_loadu_ps(&in[i].x);
        const __m256 x = _mm256_permutevar8x32_ps(p, ix);
        const __m256 y = _mm256_permutevar8x32_ps(p, iy);
        const __m256 z = _mm256_permutevar8x32_ps(p, iz);
        __m256 r = _mm256_mul_ps(c0, x);
        r = _mm256_add_ps(r, _mm256_mul_ps(c1, y));
        r = _mm256_add_ps(r, _mm256_mul_ps(c2, z));
        r = _mm256_add_ps(r, c3);
        r = _mm256_blend_ps(r, wone, 0x88);
        _mm256_storeu_ps(&out[i].x, r);
    }
    for (; i < n; ++i) {
        float x = in[i].x, y = in[i].y, z = in[i].z;
        out[i].x = M.m[0]*x + M.m[1]*y + M.m[2]*z + M.m[3];
        out[i].y = M.m[4]*x + M.m[5]*y + M.m[6]*z + M.m[7];
        out[i].z = M.m[8]*x + M.m[9]*y + M.m[10]*z + M.m[11];
        out[i].w = 1.0f;
    }
}

/* ---------------- SoA ---------------- */
__attribute__((noinline)) void tf_soa_plain(const Mat& M, const float* __restrict X, const float* __restrict Y, const float* __restrict Z,
                                            float* __restrict OX, float* __restrict OY, float* __restrict OZ, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        const float x = X[i], y = Y[i], z = Z[i];
        OX[i] = M.m[0]*x + M.m[1]*y + M.m[2]*z + M.m[3];
        OY[i] = M.m[4]*x + M.m[5]*y + M.m[6]*z + M.m[7];
        OZ[i] = M.m[8]*x + M.m[9]*y + M.m[10]*z + M.m[11];
    }
}
__attribute__((noinline)) void tf_soa_plain2(const Mat& M, const float* __restrict X, const float* __restrict Y, const float* __restrict Z,
                                             float* __restrict OX, float* __restrict OY, float* __restrict OZ, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        const float x = X[i], y = Y[i], z = Z[i];
        OX[i] = M.m[0]*x + M.m[1]*y + M.m[2]*z + M.m[3];
        OY[i] = M.m[4]*x + M.m[5]*y + M.m[6]*z + M.m[7];
        OZ[i] = M.m[8]*x + M.m[9]*y + M.m[10]*z + M.m[11];
    }
}
__attribute__((target("avx2"))) __attribute__((noinline))
void tf_soa_avx2(const Mat& M, const float* __restrict X, const float* __restrict Y, const float* __restrict Z,
                 float* __restrict OX, float* __restrict OY, float* __restrict OZ, size_t n) {
    const __m256 a0=_mm256_set1_ps(M.m[0]), a1=_mm256_set1_ps(M.m[1]), a2=_mm256_set1_ps(M.m[2]), a3=_mm256_set1_ps(M.m[3]);
    const __m256 b0=_mm256_set1_ps(M.m[4]), b1=_mm256_set1_ps(M.m[5]), b2=_mm256_set1_ps(M.m[6]), b3=_mm256_set1_ps(M.m[7]);
    const __m256 c0=_mm256_set1_ps(M.m[8]), c1=_mm256_set1_ps(M.m[9]), c2=_mm256_set1_ps(M.m[10]), c3=_mm256_set1_ps(M.m[11]);
    size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        const __m256 x = _mm256_loadu_ps(X+i), y = _mm256_loadu_ps(Y+i), z = _mm256_loadu_ps(Z+i);
        __m256 r = _mm256_mul_ps(a0,x);
        r = _mm256_add_ps(r, _mm256_mul_ps(a1,y));
        r = _mm256_add_ps(r, _mm256_mul_ps(a2,z));
        r = _mm256_add_ps(r, a3);
        _mm256_storeu_ps(OX+i, r);
        r = _mm256_mul_ps(b0,x);
        r = _mm256_add_ps(r, _mm256_mul_ps(b1,y));
        r = _mm256_add_ps(r, _mm256_mul_ps(b2,z));
        r = _mm256_add_ps(r, b3);
        _mm256_storeu_ps(OY+i, r);
        r = _mm256_mul_ps(c0,x);
        r = _mm256_add_ps(r, _mm256_mul_ps(c1,y));
        r = _mm256_add_ps(r, _mm256_mul_ps(c2,z));
        r = _mm256_add_ps(r, c3);
        _mm256_storeu_ps(OZ+i, r);
    }
    for (; i < n; ++i) {
        const float x=X[i], y=Y[i], z=Z[i];
        OX[i] = M.m[0]*x + M.m[1]*y + M.m[2]*z + M.m[3];
        OY[i] = M.m[4]*x + M.m[5]*y + M.m[6]*z + M.m[7];
        OZ[i] = M.m[8]*x + M.m[9]*y + M.m[10]*z + M.m[11];
    }
}

static double med(std::vector<double> v) { std::sort(v.begin(), v.end()); return v[v.size()/2]; }

int main() {
    Mat M; const float md[12] = {0.87f,-0.31f,0.38f,1.5f, 0.42f,0.90f,-0.11f,-2.5f, -0.30f,0.28f,0.91f,0.75f};
    memcpy(M.m, md, sizeof md);

    const size_t sizes[] = {8192, 32768, 131072, 1048576, 8388608};
    std::printf("%-6s %-10s %12s %12s %8s %8s %s\n", "layout","N","plain","avx2","ratio","A/A","bytes/pt");
    std::printf("%s\n", "------------------------------------------------------------------------------");

    for (size_t n : sizes) {
        V4 *in, *o1, *o2, *o3;
        float *X,*Y,*Z,*OX,*OY,*OZ,*QX,*QY,*QZ;
        posix_memalign((void**)&in,64,n*16); posix_memalign((void**)&o1,64,n*16);
        posix_memalign((void**)&o2,64,n*16); posix_memalign((void**)&o3,64,n*16);
        posix_memalign((void**)&X,64,n*4); posix_memalign((void**)&Y,64,n*4); posix_memalign((void**)&Z,64,n*4);
        posix_memalign((void**)&OX,64,n*4); posix_memalign((void**)&OY,64,n*4); posix_memalign((void**)&OZ,64,n*4);
        posix_memalign((void**)&QX,64,n*4); posix_memalign((void**)&QY,64,n*4); posix_memalign((void**)&QZ,64,n*4);
        for (size_t i=0;i<n;i++){ in[i].x=(float)(i%17)*0.37f; in[i].y=(float)(i%13)*0.11f; in[i].z=(float)(i%7)*0.53f; in[i].w=1.f;
            X[i]=in[i].x; Y[i]=in[i].y; Z[i]=in[i].z; }

        /* correctness + bit-exactness check */
        tf_aos_plain(M,in,o1,n); tf_aos_avx2(M,in,o2,n);
        int bitdiff=0; float maxrel=0;
        for(size_t i=0;i<n;i++){ if(memcmp(&o1[i],&o2[i],16)) {bitdiff++; for(int k=0;k<3;k++){float a=((float*)&o1[i])[k],b=((float*)&o2[i])[k]; float d=fabsf(a-b)/(fabsf(a)+1e-30f); if(d>maxrel)maxrel=d;} } }

        std::vector<double> ta,tb,taa;
        const int R = 15;
        for (int r=0;r<R;r++){
            double t0=now_ns(); tf_aos_plain (M,in,o1,n); double t1=now_ns();
            double t2=now_ns(); tf_aos_avx2  (M,in,o2,n); double t3=now_ns();
            double t4=now_ns(); tf_aos_plain2(M,in,o3,n); double t5=now_ns();
            ta.push_back((t1-t0)/n); tb.push_back((t3-t2)/n); taa.push_back((t5-t4)/n);
        }
        double a=med(ta), b=med(tb), aa=med(taa);
        std::printf("%-6s %-10zu %12.4f %12.4f %8.3f %8.3f %zu%s\n","AoS",n,a,b,a/b, a/aa, (size_t)32, bitdiff? "  [bits differ]":"  [bit-exact]");
        if (bitdiff) std::printf("       max rel diff %.3e, differing elements %d/%zu\n", maxrel, bitdiff, n);

        std::vector<double> sa,sb,saa;
        for (int r=0;r<R;r++){
            double t0=now_ns(); tf_soa_plain (M,X,Y,Z,OX,OY,OZ,n); double t1=now_ns();
            double t2=now_ns(); tf_soa_avx2  (M,X,Y,Z,QX,QY,QZ,n); double t3=now_ns();
            double t4=now_ns(); tf_soa_plain2(M,X,Y,Z,OX,OY,OZ,n); double t5=now_ns();
            sa.push_back((t1-t0)/n); sb.push_back((t3-t2)/n); saa.push_back((t5-t4)/n);
        }
        double x=med(sa), y=med(sb), xx=med(saa);
        std::printf("%-6s %-10zu %12.4f %12.4f %8.3f %8.3f %zu\n","SoA",n,x,y,x/y, x/xx, (size_t)24);

        free(in);free(o1);free(o2);free(o3);
        free(X);free(Y);free(Z);free(OX);free(OY);free(OZ);free(QX);free(QY);free(QZ);
    }
    return 0;
}
