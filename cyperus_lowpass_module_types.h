/*
 * cyperus_lowpass_module_types.h
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

#ifndef RTW_HEADER_cyperus_lowpass_module_types_h_
#define RTW_HEADER_cyperus_lowpass_module_types_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#ifndef struct_tag_spGKsvEVm7uA89hv31XX4LH
#define struct_tag_spGKsvEVm7uA89hv31XX4LH

struct tag_spGKsvEVm7uA89hv31XX4LH
{
  uint32_T MissingPlacement;
  uint32_T ComparisonMethod;
};

#endif                                 /*struct_tag_spGKsvEVm7uA89hv31XX4LH*/

#ifndef typedef_spGKsvEVm7uA89hv31XX4LH_cyper_T
#define typedef_spGKsvEVm7uA89hv31XX4LH_cyper_T

typedef struct tag_spGKsvEVm7uA89hv31XX4LH spGKsvEVm7uA89hv31XX4LH_cyper_T;

#endif                               /*typedef_spGKsvEVm7uA89hv31XX4LH_cyper_T*/

#ifndef struct_tag_skA4KFEZ4HPkJJBOYCrevdH
#define struct_tag_skA4KFEZ4HPkJJBOYCrevdH

struct tag_skA4KFEZ4HPkJJBOYCrevdH
{
  uint32_T SafeEq;
  uint32_T Absolute;
  uint32_T NaNBias;
  uint32_T NaNWithFinite;
  uint32_T FiniteWithNaN;
  uint32_T NaNWithNaN;
};

#endif                                 /*struct_tag_skA4KFEZ4HPkJJBOYCrevdH*/

#ifndef typedef_skA4KFEZ4HPkJJBOYCrevdH_cyper_T
#define typedef_skA4KFEZ4HPkJJBOYCrevdH_cyper_T

typedef struct tag_skA4KFEZ4HPkJJBOYCrevdH skA4KFEZ4HPkJJBOYCrevdH_cyper_T;

#endif                               /*typedef_skA4KFEZ4HPkJJBOYCrevdH_cyper_T*/

#ifndef struct_md4c3b2b4d3766deb09f5f407fd721e562
#define struct_md4c3b2b4d3766deb09f5f407fd721e562

struct md4c3b2b4d3766deb09f5f407fd721e562
{
  int32_T S0_isInitialized;
  real_T W0_FILT_STATES[16];
  int32_T W1_PreviousNumChannels;
  real_T P0_ICRTP;
};

#endif                             /*struct_md4c3b2b4d3766deb09f5f407fd721e562*/

#ifndef typedef_dsp_BiquadFilter_1_cyperus_lo_T
#define typedef_dsp_BiquadFilter_1_cyperus_lo_T

typedef struct md4c3b2b4d3766deb09f5f407fd721e562
  dsp_BiquadFilter_1_cyperus_lo_T;

#endif                               /*typedef_dsp_BiquadFilter_1_cyperus_lo_T*/

#ifndef struct_mdWUmD80HR4BS85gk2pz9XaC
#define struct_mdWUmD80HR4BS85gk2pz9XaC

struct mdWUmD80HR4BS85gk2pz9XaC
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  dsp_BiquadFilter_1_cyperus_lo_T cSFunObject;
};

#endif                                 /*struct_mdWUmD80HR4BS85gk2pz9XaC*/

#ifndef typedef_dspcodegen_BiquadFilter_cyper_T
#define typedef_dspcodegen_BiquadFilter_cyper_T

typedef struct mdWUmD80HR4BS85gk2pz9XaC dspcodegen_BiquadFilter_cyper_T;

#endif                               /*typedef_dspcodegen_BiquadFilter_cyper_T*/

#ifndef struct_tag_sJCxfmxS8gBOONUZjbjUd9E
#define struct_tag_sJCxfmxS8gBOONUZjbjUd9E

struct tag_sJCxfmxS8gBOONUZjbjUd9E
{
  boolean_T CaseSensitivity;
  boolean_T StructExpand;
  char_T PartialMatching[6];
  boolean_T IgnoreNulls;
};

#endif                                 /*struct_tag_sJCxfmxS8gBOONUZjbjUd9E*/

#ifndef typedef_sJCxfmxS8gBOONUZjbjUd9E_cyper_T
#define typedef_sJCxfmxS8gBOONUZjbjUd9E_cyper_T

typedef struct tag_sJCxfmxS8gBOONUZjbjUd9E sJCxfmxS8gBOONUZjbjUd9E_cyper_T;

#endif                               /*typedef_sJCxfmxS8gBOONUZjbjUd9E_cyper_T*/

#ifndef struct_mdkrVvm4QB5TjYRI09eGPp6C
#define struct_mdkrVvm4QB5TjYRI09eGPp6C

struct mdkrVvm4QB5TjYRI09eGPp6C
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T TunablePropsChanged;
  real_T PrivateSampleRate;
  real_T PassbandFrequency;
  real_T Amplitude;
  real_T Slope;
  real_T Num[12];
  real_T Den[8];
  dspcodegen_BiquadFilter_cyper_T Biquad;
};

#endif                                 /*struct_mdkrVvm4QB5TjYRI09eGPp6C*/

#ifndef typedef_CyperusLowpassFilter_cyperus__T
#define typedef_CyperusLowpassFilter_cyperus__T

typedef struct mdkrVvm4QB5TjYRI09eGPp6C CyperusLowpassFilter_cyperus__T;

#endif                               /*typedef_CyperusLowpassFilter_cyperus__T*/

/* Parameters (default storage) */
typedef struct P_cyperus_lowpass_module_T_ P_cyperus_lowpass_module_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_cyperus_lowpass_modul_T RT_MODEL_cyperus_lowpass_modu_T;

#endif                          /* RTW_HEADER_cyperus_lowpass_module_types_h_ */
