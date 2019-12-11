/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF Skeleton
Module: BSF-Forwards.h (Problem Independent Function Forwards)
Author: Leonid B. Sokolinsky 

This source code is a part of BSF Skeleton
==============================================================================*/
#include "BSF-Types.h"					// Problem Independent Types 

static void	BC_Init(bool* success);		// BSF Initialization
static void BC_Master();				// Master Process
static void BC_MasterMap(bool exit);
static void BC_MasterReduce();
static void BC_MpiRun();
static void BC_ProcessExtendedReduceList(BT_extendedReduceElem_T* reduceList, int index, int length, 
	BT_extendedReduceElem_T** extendedReduceResult_P);
static void BC_Worker();	// Worker Process
static bool BC_WorkerMap();
static void BC_WorkerReduce();

