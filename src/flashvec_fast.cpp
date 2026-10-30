/* Approximate math entries -- currently fv_atan_approx.
 *
 * The polynomial is the standard 5-term atan approximation on [-1, 1]:
 *   atan(x) ~= x*(0.9998660 + x^2*(-0.3302995 + x^2*(0.1801410 + x^2*(-0.0851330 + x^2*0.0208351))))
 * evaluated with Horner's rule. The measured max absolute error against atanf is about 1.2e-5.
 *
 * Outside [-1, 1] (and for NaN) it returns NaN instead of clamping: a clamped angle would be silently
 * wrong, while the caller can reduce |x| > 1 with atan(x) = sign(x)*pi/2 - atan(1/x). This is an
 * approximate entry point; it does not promise libm's result. The unit carries fp-contract=off so the
 * polynomial is evaluated as written, but no bit-reproducibility claim is made for it.
 */
#include "flashvec.h"

#include <math.h>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("fp-contract=off")
#endif

fv_real fv_atan_approx(fv_real x)
{
    if (!(x >= (fv_real)-1 && x <= (fv_real)1)) {
        return (fv_real)NAN;
    }
    {
        fv_real x2 = x * x;
        fv_real p = (fv_real)0.0208351;
        p = p * x2 + (fv_real)-0.0851330;
        p = p * x2 + (fv_real)0.1801410;
        p = p * x2 + (fv_real)-0.3302995;
        p = p * x2 + (fv_real)0.9998660;
        return x * p;
    }
}
