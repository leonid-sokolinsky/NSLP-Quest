/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF Skeleton
Module: BSF-ProblemFunctions.h (Predefined Problem Function Forwards)
Author: Leonid B. Sokolinsky

This source code is a part of BSF Skeleton
==============================================================================*/
#include "Problem-bsfTypes.h"	// Predefined Problem Types

void PC_bsf_AssignListSize(int* listSize);
void PC_bsf_CopyData(PT_bsf_data_T* dataIn, PT_bsf_data_T* dataOut);
void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_data_T data,
	int iterCounter, double elapsedTime);
void PC_bsf_Init(bool* success);
void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, PT_bsf_data_T* data, int numberInSublist, 
	int sublistLength, int offset, int* success);
void PC_bsf_ParametersOutput(int numOfWorkers, PT_bsf_data_T data);
void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_data_T data,
	int iterCounter, double t);
void PC_bsf_ProcessResults(bool* exit, PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_data_T* data, int iterCounter);
void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z);
void PC_bsf_SetInitData(PT_bsf_data_T* data);
void PC_bsf_SetMapSubList(PT_bsf_mapElem_T* sublist, int sublistLength, int offset, bool* success);
