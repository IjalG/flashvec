/* bench/ab/ab_scalar_lib.cpp -- the **library side** of the "before the refactor"
 * anchor: the HEAD snapshot's src/flashvec.cpp compiled as-is into its own TU, so the
 * element-wise scalar API has real external symbols here.
 *
 * Why two TUs (this file + ab_scalar_loop.cpp): only by calling it from **another** TU can
 * GCC be kept from seeing the body -- no inlining, no vectorization -- which is what a
 * consumer linking a prebuilt static library really gets (tier 1 in the README). In a
 * single TU it inlines into tier 2 (scalar_flat), and that number is not an anchor.
 */
#include "snapshot/include/flashvec.h"
#include "snapshot/src/flashvec.cpp"
