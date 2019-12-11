/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF Skeleton
Module: BSF-Data.h (Problem Independent Variables & Data Structures)
Prefix: BD
Author: Leonid B. Sokolinsky 

This source code is a part of BSF Skeleton
==============================================================================*/
#include "BSF-Include.h"		// Problem Independent "Include" Files
#include "BSF-Types.h"			// Problem Independent Types
// ========================= Constants ========================
#define BD_EXIT	true			// End of Work
// ========================= Variables ========================
static int	BD_listSize;		// Size of the list
static int	BD_rank;			// Rank of the current process MPI 
static int	BD_size;			// Number of MPI processes
static int	BD_masterRank;		// Rank of Master = Number of MPI processes - 1
static int	BD_numOfWorkers;	// Number of Workers = Number of MPI processes - 1
static int	BD_elemsPerWorker;	// Number of elements per worker
static int	BD_tailLength;		// Length of the remainder of the list after dividing by the number of workers
static bool BD_exit;			// Indicator of ending the calculations
static bool BD_success;			// Indicator of Successfulness of Initialization
// ========================= Time Variables ========================
static double BD_t;				// Total time
static int BD_iterCounter;
// ========================= Data stuctures ========================
static PT_bsf_data_T BD_data;
// ========================= Lists ========================
PT_bsf_mapElem_T* BD_mapSubList;
BT_extendedReduceElem_T* BD_extendedReduceList; 
BT_extendedReduceElem_T* BD_extendedReduceResult_P;
BT_order_T* BD_order;		// Orders for Workers
MPI_Status*  BD_status;		// MPI ststus
MPI_Request* BD_request;	// MPI request
int* BD_sublistSize;		// SubList Sizes for Workers
int* BD_offset;				// List offsets for Workers