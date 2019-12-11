/*==============================================================================
Project: NSLP (Non-Stationary Linear Programming)
Theme: Quest Phase
Module: Problem-Data.h (Problem Data)
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Types.h"		// Problem Parameters 
//========================== Problem variables ====================================

//========================== Problem structures ====================================
static float PD_normSquare_a[PP_M];

/*------------ Way_b pararters (Fourier series) ---------------*/
static float* PD_direction;			// Amplitudes

static float PD_previousShift[PP_N];
