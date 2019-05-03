/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF-MRA Skeleton
Module: Problem-bsf-Forwards.h (Predefinite Problem Function Forwards)
Author(s): Leonid B. Sokolinsky
This source code belonges to the BSF-skeleton
==============================================================================*/
#pragma once
void PI_bsf_AssignListSize(int* listSize);
void PI_bsf_CopyData(PT_bsf_data_T* dataIn, PT_bsf_data_T* dataOut);
void PI_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int count, PT_bsf_data_T data,
	int iterCount, double elapsedTime);
void PI_bsf_Init(bool* success);
void PI_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int index, PT_bsf_data_T* data, int* success);
void PI_bsf_ParametersOutput(int numOfWorkers, PT_bsf_data_T data);
void PI_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int count, PT_bsf_data_T data,
	int iterCount, double t, double t_L, double t_s_L, double t_S, double t_r_L,
	double BD_t_w, double BD_t_A_w, double BD_t_A_m, double BD_t_p);
void PI_bsf_ProcessResults(bool* exit, PT_bsf_reduceElem_T* reduceResult, int count, PT_bsf_data_T* data);
void PI_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z);
void PI_bsf_SetInitApproximation(PT_bsf_data_T* data);
void PI_bsf_SetMapSubList(PT_bsf_mapElem_T* subList, int count, int offset, bool* success);