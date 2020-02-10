/*==============================================================================
Project: NSLP (Non-Stationary Linear Programming)
Theme: Quest Phase
Module: Problem-bsfParameters.h (BSF Skeleton Parameters)
Prefix: PP_BSF
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/

//=========================== BSF Skeleton Parameters =========================
#define PP_BSF_PRECISION 5		// Decimal precision on output
#define PP_BSF_ITER_OUTPUT		// If PP_BSF_ITER_OUTPUT is defined then Iteration Output is performed
#define PP_BSF_TRACE_COUNT 1	// Each PP_BSF_TRACE_COUNT-th iteration to be outputted
#define PP_BSF_MAX_JOB_CASE 0
#define PP_BSF_FRAGMENTED_MAP_LIST	// If it is defined, the worker-node stores only its part of the map-list. 
										// Otherwise, the map-list is replicated on each worker-node. 
//--------------------------- OpenMP Parameters ---------------------------
//#define PP_BSF_OMP				// If PP_BSF_OMP is defined then OpenMP is turned on for Map Step
//#define PP_BSF_NUM_THREADS 12	// If PP_BSF_NUM_THREADS is udefined then all accessable threads are used
//=========================== Skeleton Variables =========================
static int PP_BSF_addressOffset;
static int PP_BSF_iterCounter;
static int PP_BSF_jobCase;
static void* PP_BSF_mapSubList;
static int PP_BSF_mpiRank;
static int PP_BSF_numberInSublist;
static int PP_BSF_numOfWorkers;
static int PP_BSF_sublistLength;