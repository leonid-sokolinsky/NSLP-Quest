/*==============================================================================
Project: NSLP (Non-Stationary Linear Programming)
Theme: Quest Phase
Module: Problem-bsfParameters.h (BSF Skeleton Parameters)
Prefix: PP_BSF
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/

//=========================== BSF Skeleton Parameters =========================
#define PP_BSF_PRECISION 5		// Decimal precision on output
#define PP_BSF_ITER_OUTPUT		// If PP_BSF_ITER_OUTPUT is defined then Iteration Output is performed
#define PP_BSF_TRACE_COUNT 1	// Each PP_BSF_TRACE_COUNT-th iteration to be outputted
//--------------------------- OpenMP Parameters ---------------------------
#define PP_BSF_OMP				// If PP_BSF_OMP is defined then OpenMP is turned on for Map Step
//#define PP_BSF_NUM_THREADS 12	// If PP_BSF_NUM_THREADS is udefined then all accessable threads are used