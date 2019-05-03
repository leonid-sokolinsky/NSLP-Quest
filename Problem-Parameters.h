/*==============================================================================
Project: NSLP (Non-Stationary Linear Programming)
Theme: Quest Phase
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
#pragma once
//=========================== Skeleton Parameters =========================
#define PP_BSF_PRECISION 4		// Decimal precision on output
#define PP_BSF_ITER_OUTPUT		// If PP_BSF_ITER_OUTPUT is defined then Iteration Output is performed
#define PP_BSF_TRACE_COUNT	100
//--------------------------- OpenMP Parameters ---------------------------
#define PP_BSF_OMP			// If PP_BSF_OMP is defined then OpenMP is turned on for Map Step
#define PP_BSF_NUM_THREADS 12	// If PP_BSF_NUM_THREADS is udefined then all accessable threads are used

//=========================== Problem Parameters =========================
//#define PP_CIMMINO
#define PP_N 54000			// Dimension of Space
#define PP_M (2*PP_N+2)		// Number of inequations
#define PP_SF	200			// Scale factor

#define PP_INIT_APPROX (-PP_SF * 0.01)

#define PP_MAX_ITER_COUNT	1E8	// Maximal count of iterations
#define PP_LAMBDA		1.99	// Relaxation factor
#define PP_EPS_RELAX	1E-6			// Precision
#define PP_EPS_IN		1E-6
#define PP_EPS_LAG		0.001

//_____________________________________________
#define PP_VELOCITY 0.0025 // Units per Second	
//_____________________________________________

//-------------------------- Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	11	// Number of Elements to output
//#define PP_MATRIX_OUTPUT	// Output Matrix
