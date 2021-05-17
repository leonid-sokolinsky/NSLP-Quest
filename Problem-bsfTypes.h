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
	PT_float_T* a;			// Pointer to row: (a_0,...,a_{n-1})
	PT_float_T* b;			// Pointer to constant term b: a_0x_0+...+a_{n-1}x_{n-1} \leq b
	PT_float_T normSquare;	// a_0^2+...+a_{n-1}^2
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