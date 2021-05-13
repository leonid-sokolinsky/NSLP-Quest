/*==============================================================================
Project: CoFePro
Theme: Projection Algorithm for Solving Convex Feasibility Problems
Module: Problem-bsfTypes.h (Predefined Problem Types)
Prefix: PT_bsf
Author: Leonid B. Sokolinsky
This source code is a part of BSF Skeleton
==============================================================================*/
#pragma once
#include "Problem-Types.h"		// Problem Types 
//=========================== BSF Types =========================
struct PT_bsf_parameter_T {				// Parameter for workers (current approximation)
	PT_vector_T x;		// Current approximation
};

struct PT_bsf_mapElem_T {	// Element of the map list
	PT_vector_T a;		// Inequality: a[0]*x_1+...+a[n-1]*x_n <= b
	PT_float_T b;
	PT_float_T normSquare; // a[0]*a[0]+...+a[n-1]*a[n-1]
};

struct PT_bsf_reduceElem_T {	// Element of reduce list	
	PT_vector_T projection;	// ProjectionVector of vector onto hyperplane
};

struct PT_bsf_reduceElem_T_1 {
	// optional filling
};

struct PT_bsf_reduceElem_T_2 {
	// optional filling
};

struct PT_bsf_reduceElem_T_3 {
	// optional filling
};