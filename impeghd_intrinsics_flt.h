#ifndef IMPEGHD_INTRINSICS_FLT_H
#define IMPEGHD_INTRINSICS_FLT_H
#define ITT_INTRINSICS
#include <math.h>
#include <stdlib.h>
#define ia_sqrt_flt sqrt
#define ia_floor_flt floor
#define ia_fabs_flt(a) ((FLOAT32)fabs(a))
#ifdef WIN32
#define ia_max_flt(a, b) (max(a, b))
#define ia_min_flt(a, b) (min(a, b))
#define ia_max_int(a, b) (max(a, b))
#define ia_min_int(a, b) (min(a, b))
#else
#define ia_max_flt(a, b) (fmax(a, b))
#define ia_min_flt(a, b) (fmin(a, b))
#define ia_max_int(a, b) (((a) > (b)) ? (a) : (b))
#define ia_min_int(a, b) (((a) < (b)) ? (a) : (b))
#define ia_eq_int(a, b) ((a) == (b) ? (1) : (0))
#endif
#define ia_abs_int abs
#define ia_negate_flt(a) (-a)
#define ia_sub_flt(a, b) ((a) - (b))
#define ia_ceil_flt(a) ceil(a)
#define ia_div_q15_flt(a) (a) / 32768
#define ia_add_flt(a, b) ((a) + (b))
#define ia_add_double(a, b) ((a) + (b))
#define ia_mul_flt(a, b) ((a) * (b))
#define ia_mul_double_flt(a, b) ((a) * (b))
#define ia_lt_flt(a, b) ((a) < (b) ? (1) : (0))
#define ia_lteq_flt(a, b) ((a) <= (b) ? (1) : (0))
#define ia_eq_flt(a, b) ((a) == (b) ? (1) : (0))
#define ia_mac_flt(x, a, b) ((x) + (a) * (b))
#define ia_msu_flt(x, a, b) ((x) - (a) * (b))
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif
#endif /*IMPEGHD_INTRINSICS_FLT_H*/