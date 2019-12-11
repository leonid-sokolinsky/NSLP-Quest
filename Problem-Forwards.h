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
static float	DotProduct(PT_point_T x, PT_point_T y);
static float	Way_b(PT_inequality_T inequality, PT_shift_T shift);
static bool		ExitCondition(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_data_T* data, int iterCounter);
static void		GetShift(PT_shift_T shift);
static float	NormSquare(PT_point_T x);
static bool		PointIn(PT_point_T point, PT_inequality_T inequality, float wayed_b);

//====================== Macros ================================
#define PF_MIN(x,y) (x<y?x:y)
#define PF_MAX(x,y) (x>y?x:y)