/*==============================================================================
Project: NSLP (Non-Stationary Linear Programming)
Theme: Quest Phase
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
#pragma once
//=========================== Problem Parameters =========================
#define PP_N 2			// Dimension of Space
#define PP_M (2*PP_N+2)		// Number of inequations
#define PP_SF	200			// Scale factor

#define PP_INIT_APPROX (-PP_SF)

#define PP_MAX_ITER_COUNT	1E5	// Maximal count of iterations
#define PP_LAMBDA		140		// Relaxation factor
#define PP_EPS_RELAX	1E-6	// Precision
#define PP_EPS_IN		1E-6
#define PP_EPS_LAG		100

//_____________________________________________
#define PP_VELOCITY 0 // Units per Second
//_____________________________________________

//-------------------------- Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	11	// Number of Elements to output
#define PP_MATRIX_OUTPUT	// Output Matrix
