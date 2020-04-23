/*
 * rtGetNaN.h
 *
 * Code generation for model "f14_acc".
 *
 * Model version              : 1.13
 * Simulink Coder version : 8.11 (R2016b Prerelease) 21-Feb-2016
 * C source code generated on : Thu Mar 10 11:15:53 2016
 *
 * Target selection: accel.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Emulation hardware selection:
 *    Differs from embedded hardware (MATLAB Host)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef dsp_math_utils_h_
#define dsp_math_utils_h_
#include <stddef.h>
#include "rtwtypes.h"
#include "rt_nonfinite.h"

extern float rtGetNaN(void);
extern float rtGetNaNF(void);
extern float rtGetInf(void);
extern float rtGetInfF(void);
extern float rtGetMinusInf(void);
extern float rtGetMinusInfF(void);

#endif                                 /* RTW_HEADER_rtGetNaN_h_ */
