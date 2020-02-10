/*==============================================================================
Project: NSLP (Non-Stationary Linear Programming)
Theme: Quest Phase
Module: Problem-Data.h (Problem Data)
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Types.h"		// Problem Parameters 
//========================== Problem variables ====================================
static float* PD_shift;
static float* PD_approximation;

//========================== Problem structures ====================================

/*------------ Way_b pararters (Fourier series) ---------------*/
static float* PD_direction;			// Amplitudes

static float PD_previousShift[PP_N];
