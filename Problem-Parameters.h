/*==============================================================================
Project: CoFePro
Theme: Projection Algorithm for Solving Convex Feasibility Problems
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
//=========================== Problem Parameters =========================
#define PP_N 5											// Dimension of Space
#define PP_NUM_OF_RND_INEQUALITIES (2 * PP_N)			// Number of random inequalities
#define PP_M (2*PP_N + PP_NUM_OF_RND_INEQUALITIES + 1)	// Total number of inequalities of given system
#define PP_SF	200										// Scale factor

#define PP_MAX_ITER_COUNT	1E5		// Maximal count of iterations
#define PP_LAMBDA			1		// Relaxation factor
#define PP_EPS_RELAX		1E-7	// Precision
#define PP_EPS_ZERO			1E-9	// Precision of comparison with zero

//-------------------------- Input/Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	11	// Number of Elements to output
#define PP_MATRIX_OUTPUT	// To output Matrix
#define PP_SETW 14
#define PP_PATH "C:/TEMP/"
#define PP_LPP_FILE "inequalities.txt"
#define PP_EXTERIOR_POINT_FILE "exteriorPoint.txt"
#define PP_INTERIOR_POINT_FILE "interiorPoint.txt"