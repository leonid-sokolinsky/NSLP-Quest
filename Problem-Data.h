/*==============================================================================
Project: NSLP (Non-Stationary Linear Programming)
Theme: Quest Phase
Module: Problem-Data.h (Problem Data)
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Types.h"		// Problem Parameters 

//=========================== Variables for BSF-skeleton Parameters =========================
static int PP_BSF_addressOffset;
static int PP_BSF_iterCounter;
static int PP_BSF_jobCase;
static int PP_BSF_mpiRank;
static int PP_BSF_numberInSublist;
static int PP_BSF_numOfWorkers;
static int PP_BSF_sublistLength;

//========================== Problem variables ====================================
static double* PD_x_P; // Pointer to Current Approximation

//========================== Problem structures ====================================
static double PD_A[PP_M][PP_N];
static double PD_b[PP_M];
static double PD_c[PP_N];		// Objective Function Coefficients
static double PD_normSquare_a[PP_M]; // a_i1*a_i1+...+a_in*a_in
static PT_point_T PD_apex;		// Apex Point