/*
 * cyperus_lowpass_module.h
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

#ifndef RTW_HEADER_cyperus_lowpass_module_h_
#define RTW_HEADER_cyperus_lowpass_module_h_
#include <float.h>
#include <math.h>
#include <string.h>
#include <stddef.h>
#ifndef cyperus_lowpass_module_COMMON_INCLUDES_
# define cyperus_lowpass_module_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_logging.h"
#include "HostLib_MMFile.h"
#include "HostLib_Multimedia.h"
#include "HostLib_Audio.h"
#endif                             /* cyperus_lowpass_module_COMMON_INCLUDES_ */

#include "cyperus_lowpass_module_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rtGetNaN.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"
#include "rtGetInf.h"

struct cyperus_parameters {
  /* block processing */
  /* rtqueue_t block_fifo; */

  /* generic */
  float freq;
  float amp;
  float amt;
  float delta;
  float mix;
  float res; /* resonance */
  float shift;
  float fb; /* feedback */
  float delay_amt; /* delay amount, 0-1 */
  float delay_time; /* init this with
		       = seconds * sample_rate */

  float attack;
  float decay;
  float scale;

  int pos;
  int delay_pos;
  float avg;
  float max;
  int vocoder_pos;

  float q;
  
  /* sine */
  int exp;
  int skip;
  int skip_count;
  float last_freq;
  float hold_freq;
  float phase_delta;
  float phase;
  
  /* phase vocoders */
  float *fft_buffer;
  float *last_phase;
  float *sum_phase;
  float *output_accumulator;
  
  float *signal_buffer;
  float *signal_buffer_out;
  
  float in;

  float state0;
  float state1;
  float state2;
  float state3;
  float state4;
  float state5;
  float state6;
  float state7;
  float state8;
  
  float tempval;

  float lastinval;
  float lastinval1;
  float lastinval2;
  float lastinval3;
  float lastoutval;
  float lastoutval1;
  float lastoutval2;
  float lastoutval3;

  float x0;
  float x1;
  float x2;
  float x3;
  float x4;

  float y0;
  float y1;
  float y2;
  float y3;
  float y4;

  int slope;
  int fs;
  float B[12];
  float A[8];

  float W0_FILT_STATES[16];
};

/* Macros for accessing real-time model data structure */
#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetRTWLogInfo
# define rtmGetRTWLogInfo(rtm)         ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  ((rtm)->Timing.taskTime0)
#endif

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               (&(rtm)->Timing.taskTime0)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T FromMultimediaFile[2048];     /* '<Root>/From Multimedia File' */
  real_T MATLABSystem[2048];           /* '<Root>/MATLAB System' */
  real_T dv0[2048];
} B_cyperus_lowpass_module_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  CyperusLowpassFilter_cyperus__T obj; /* '<Root>/MATLAB System' */
  real_T FromMultimediaFile_HostLib[137];/* '<Root>/From Multimedia File' */
  real_T FromMultimediaFile_AudioInfo[5];/* '<Root>/From Multimedia File' */
  real_T FromMultimediaFile_VideoInfo[11];/* '<Root>/From Multimedia File' */
  uint64m_T thisPtr;                   /* '<Root>/MATLAB System' */
  uint8_T AudioDeviceWriter_AudioDeviceLi[1096];/* '<Root>/Audio Device Writer' */
  boolean_T objisempty;                /* '<Root>/MATLAB System' */
  boolean_T thisPtr_not_empty;         /* '<Root>/MATLAB System' */
  boolean_T isInitialized;             /* '<Root>/MATLAB System' */
} DW_cyperus_lowpass_module_T;

/* Parameters (default storage) */
struct P_cyperus_lowpass_module_T_ {
  real_T MATLABSystem_PassbandFrequency;/* Expression: 672.9902424112955
                                         * Referenced by: '<Root>/MATLAB System'
                                         */
  real_T MATLABSystem_Amplitude;       /* Expression: 0.2027734502156575
                                        * Referenced by: '<Root>/MATLAB System'
                                        */
  real_T MATLABSystem_Slope;           /* Expression: 10
                                        * Referenced by: '<Root>/MATLAB System'
                                        */
  real_T AudioDeviceWriter_P1;         /* Expression: 0
                                        * Referenced by: '<Root>/Audio Device Writer'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_cyperus_lowpass_modul_T {
  const char_T *errorStatus;
  RTWLogInfo *rtwLogInfo;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T taskTime0;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    time_T tFinal;
    boolean_T stopRequestedFlag;
  } Timing;
};

/* Block parameters (default storage) */
extern P_cyperus_lowpass_module_T cyperus_lowpass_module_P;

/* Block signals (default storage) */
extern B_cyperus_lowpass_module_T cyperus_lowpass_module_B;

/* Block states (default storage) */
extern DW_cyperus_lowpass_module_T cyperus_lowpass_module_DW;

/* Model entry point functions */
extern void cyperus_lowpass_module_initialize(void);
extern void cyperus_lowpass_module_step(void);
extern void cyperus_lowpass_module_terminate(void);

/* Real-time Model object */
extern RT_MODEL_cyperus_lowpass_modu_T *const cyperus_lowpass_module_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'cyperus_lowpass_module'
 */
#endif                                /* RTW_HEADER_cyperus_lowpass_module_h_ */
