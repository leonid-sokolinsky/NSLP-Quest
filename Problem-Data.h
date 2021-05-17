/*==============================================================================
Project: CoFePro
Theme: Projection Algorithm for Solving Convex Feasibility Problems
Module: Problem-Data.h (Problem Data)
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Types.h"		// Problem Parameters 
using namespace std;
//========================== INput/Output ====================================
static string PD_inequalitiesFile; /* LPP file in the following format:
------------ begin of file -------------
m n
A_11 A_12 ... A_1n b_1
A_21 A_22 ... A_2n b_2
...
A_m1 A_m2 ... A_mn b_m
------------ end of file----------------*/

static string PD_exteriorPointFile; /* Exterior point file in the following format:
------------ begin of file -------------
n
x_1 x_2 ... x_n
------------ end of file----------------*/

static string PD_interiorPointFile; /* Interior point file in the following format:
------------ begin of file -------------
n
x_1 x_2 ... x_n
------------ end of file----------------*/

//========================== Problem structures ====================================
static PT_matrix_T PD_A;					// Matrix of coefficients of inequalities 
static PT_column_T PD_b;					// Column of the constant terms of the system Ax <= PD_b
static PT_vector_T PD_exteriorPoint;		