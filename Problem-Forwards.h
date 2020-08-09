/*==============================================================================
Project: NSLP (Non-Stationary Linear Programming)
Theme: Quest Phase
Module: Problem-Forwards.h (Problem Function Forwards)
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-bsfTypes.h"
#include "Problem-Types.h"
//====================== Problem Functions ===========================
static double	DotProduct(PT_point_T x, PT_point_T y);
static bool		ExitCondition(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter);
static double	NormSquare(PT_point_T x);

//====================== Macros ================================
#define PF_MIN(x,y) (x<y?x:y)
#define PF_MAX(x,y) (x>y?x:y)