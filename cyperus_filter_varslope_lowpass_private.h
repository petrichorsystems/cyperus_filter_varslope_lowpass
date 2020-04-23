/*
 * cyperus_lowpass_module_private.h
 *
 * Code generation for model "cyperus_lowpass_module".
 *
 * Model version              : 1.2
 * Simulink Coder version : 9.1 (R2019a) 23-Nov-2018
 * C source code generated on : Mon Apr 13 12:16:16 2020
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objective: Execution efficiency
 * Validation result: Not run
 */

#ifndef RTW_HEADER_cyperus_filter_varslope_lowpass_private_h_
#define RTW_HEADER_cyperus_filter_varslope_lowpass_private_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmSetTFinal
# define rtmSetTFinal(rtm, val)        ((rtm)->Timing.tFinal = (val))
#endif

#ifndef UCHAR_MAX
#include <limits.h>
#endif

extern float rt_roundd_snf(float u);
extern float rt_remd_snf(float u0, float u1);
extern float rt_hypotd_snf(float u0, float u1);
extern float rt_atan2d_snf(float u0, float u1);

#endif

