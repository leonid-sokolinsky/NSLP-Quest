/*==============================================================================
Project: NSLP (Non-Stationary Linear Programming)
Theme: Quest Phase
Module: Problem-Types.h (BSF Types)
Prefix: PT
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
#pragma once						
#include "Problem-Parameters.h"		// Problem Parameters 

//=========================== Problem Types =========================
typedef float PT_point_T[PP_N];			// Point in n-Dimensional Space 
typedef float PT_shift_T[PP_N];	// Shift (non-stationarity)
typedef float PT_inequality_T[PP_N + 1];		// Inequality: inequality[0]*x_1+...+inequality[n-1]*x_n <= inequality[n]
//---------------------- Way_b Types (Fourier series) -------------------------
struct PT_amplitude_T {	// Amplitudes
	float _1;
	float _2;
	float _3;
};
struct PT_theta_T {	// Phases
	float _1;
	float _2;
	float _3;
};
//=========================== BSF Types =========================
struct PT_bsf_data_T {				// Data for workers (current approximation)
	float approximation[PP_N];		// Current approximation
	float shift[PP_N];	// Way_b of Polytope
};

struct PT_bsf_mapElem_T {			// Element of the map list
	float* inequality;	// Inequality: inequality[0]*x_1+...+inequality[n-1]*x_n <= inequality[n] 
};

struct PT_bsf_reduceElem_T {		// Element of reduce list	
	float point[PP_N];			// Current Approximation
	bool pointIn;				// Point is in Polytope
};