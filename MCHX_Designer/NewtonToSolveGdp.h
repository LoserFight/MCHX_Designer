//
// MATLAB Compiler: 7.0.1 (R2019a)
// Date: Thu Jun  2 09:34:56 2022
// Arguments:
// "-B""macro_default""-W""cpplib:NewtonToSolveGdp""-T""link:lib""NewtonToSolveG
// dp.m""-C"
//

#ifndef __NewtonToSolveGdp_h
#define __NewtonToSolveGdp_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#include "mclcppclass.h"
#ifdef __cplusplus
extern "C" {
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_NewtonToSolveGdp_C_API 
#define LIB_NewtonToSolveGdp_C_API /* No special import/export declaration */
#endif

/* GENERAL LIBRARY FUNCTIONS -- START */

extern LIB_NewtonToSolveGdp_C_API 
bool MW_CALL_CONV NewtonToSolveGdpInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_NewtonToSolveGdp_C_API 
bool MW_CALL_CONV NewtonToSolveGdpInitialize(void);

extern LIB_NewtonToSolveGdp_C_API 
void MW_CALL_CONV NewtonToSolveGdpTerminate(void);

extern LIB_NewtonToSolveGdp_C_API 
void MW_CALL_CONV NewtonToSolveGdpPrintStackTrace(void);

/* GENERAL LIBRARY FUNCTIONS -- END */

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_NewtonToSolveGdp_C_API 
bool MW_CALL_CONV mlxNewtonToSolveGdp(int nlhs, mxArray *plhs[], int nrhs, mxArray 
                                      *prhs[]);

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

#ifdef __cplusplus
}
#endif


/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__MINGW64__)

#ifdef EXPORTING_NewtonToSolveGdp
#define PUBLIC_NewtonToSolveGdp_CPP_API __declspec(dllexport)
#else
#define PUBLIC_NewtonToSolveGdp_CPP_API __declspec(dllimport)
#endif

#define LIB_NewtonToSolveGdp_CPP_API PUBLIC_NewtonToSolveGdp_CPP_API

#else

#if !defined(LIB_NewtonToSolveGdp_CPP_API)
#if defined(LIB_NewtonToSolveGdp_C_API)
#define LIB_NewtonToSolveGdp_CPP_API LIB_NewtonToSolveGdp_C_API
#else
#define LIB_NewtonToSolveGdp_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_NewtonToSolveGdp_CPP_API void MW_CALL_CONV NewtonToSolveGdp(int nargout, mwArray& r, mwArray& F, mwArray& R, const mwArray& x0, const mwArray& RL, const mwArray& RT, const mwArray& Gsum, const mwArray& eps);

/* C++ INTERFACE -- WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */
#endif

#endif
