/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF Skeleton
Module: BSF-Types.h (Problem Independent Types)
Prefix: BT
Author: Leonid B. Sokolinsky 

This source code is a part of BSF Skeleton
==============================================================================*/
#pragma once
#include "Problem-bsfTypes.h"	// Predefined BSF Problem Types
//=========================== BSF Types ===============================
struct BT_order_T {  
	char exit;		// true, if worker must stop
	PT_bsf_data_T data;
};

struct BT_extendedReduceElem_T {// Extended element type of reduce list
	PT_bsf_reduceElem_T elem;	// Element of reduce list
	int reduceCounter;			// Reduce Counter
};
