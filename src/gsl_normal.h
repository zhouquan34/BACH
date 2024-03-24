/* The function gsl_cdf_ugaussian_P is copied from GSL library */

#ifndef GSL_NORMAL
#define GSL_NORMAL

#include <cmath>

/* From GSL: gsl_machine.h */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
// end

/* From GSL: gsl_math.h */ 
#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif
//end

double gsl_cdf_ugaussian_P (const double x);

#endif 

