/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF-MRA Skeleton
Module: BSF-Forwards.h (BSF Function Forwards)
Author(s): Leonid B. Sokolinsky
This source code belonges to the BSF-skeleton
==============================================================================*/
#pragma once
static void	BI_Init(bool* success);		// BSF Initialization
static void BI_Master();	// Master Process
static void BI_MasterMap(bool exit);
static void BI_MasterReduce();
static void BI_MeasureTimeParameters();
static void BI_ProcessExtendedReduceList(BT_extendedReduceElem_T* reduceList, int index, int length, 
	BT_extendedReduceElem_T** extendedReduceResult_P);
static void BI_Worker();	// Worker Process
static bool BI_WorkerMap();
static void BI_WorkerReduce();

static void MpiRun();
