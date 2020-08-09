/*==============================================================================
Project: NSLP (Non-Stationary Linear Programming)
Theme: Quest Phase
Module: Problem-Types.h (BSF Types)
Prefix: PT
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/			
#pragma once
#include "Problem-Include.h"		// Problem "Include" Files
#include "Problem-Parameters.h"		// Problem Parameters 
//=========================== Problem Types =========================
typedef double PT_point_T[PP_N];			// Point in n-Dimensional Space 
typedef double PT_inequality_T[PP_N + 1];		// Inequality: inequality[0]*x_1+...+inequality[n-1]*x_n <= inequality[n]