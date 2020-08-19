/*==============================================================================
Project: NSLP (Non-Stationary Linear Programming)
Theme: Quest Phase
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
//=========================== Problem Parameters =========================
#define PP_N 10			// Dimension of Space
#define PP_M (2*PP_N+2)		// Number of inequations
#define PP_SF	200			// Scale factor

#define PP_MAX_ITER_COUNT	1E5	// Maximal count of iterations
#define PP_LAMBDA		1		// Relaxation factor
#define PP_EPS_RELAX	1E-6	// Precision
#define PP_DIST_TO_APEX 1E5	// Distance to Apex Point 

//-------------------------- Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	11	// Number of Elements to output
#define PP_MATRIX_OUTPUT	// Output Matrix
