/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF Skeleton
Module: BSF-Code.cpp (Problem Independent Code)
Prefix: BI
Author: Leonid B. Sokolinsky 

This source code is a part of BSF Skeleton
==============================================================================*/
#include "BSF-Data.h"				// Problem Independent Variables & Data Structures 
#include "BSF-Forwards.h"			// Problem Independent Function Forwards
#include "BSF-ProblemFunctions.h"	// Predefined Problem Function Forwards
#include "Problem-bsfParameters.h"	// Predefined Problem Parameters
using namespace std;

int main(int argc, char *argv[]) {
	char emptystring[] = "";
	char* message = emptystring;
	unsigned success;

	BC_MpiRun();
	
	BD_success = true;
	PC_bsf_Init(&BD_success);
	MPI_Allreduce(&BD_success, &success, 1, MPI_UNSIGNED, MPI_LAND, MPI_COMM_WORLD);
	if (!success) {
		if (BD_rank == BD_masterRank) 
			cout << "Error: PC_bsf_Init failed!" << endl;
		MPI_Finalize();
		exit(1);
	};

	BD_success = true; 
	BC_Init(&BD_success);
	MPI_Allreduce(&BD_success, &success, 1, MPI_UNSIGNED, MPI_LAND, MPI_COMM_WORLD);
	if (!success) {
		if (BD_rank == BD_masterRank) cout << "Error: BC_Init failed (not enough memory)!" << endl;
		MPI_Finalize();
		exit(1);
	};

	if (BD_rank == BD_masterRank)
		BC_Master();
	else
		BC_Worker();
	
	MPI_Finalize();
	return 0;
};
static void BC_Master() {// Master Process
	PC_bsf_ParametersOutput(BD_numOfWorkers, BD_parameter);
	BD_iterCounter = 0;

	BD_t -= MPI_Wtime();
	do {
		BD_jobCase = BD_newJobCase;
		if (BD_jobCase > PP_BSF_MAX_JOB_CASE) {
			cout << "BC_Master:Error: Job Case = " << BD_jobCase << " > PP_BSF_MAX_JOB_CASE = " << PP_BSF_MAX_JOB_CASE << endl;
			BD_exit = BD_EXIT;
			break;
		};
		BC_MasterMap(!BD_EXIT);
		BC_MasterReduce();
		switch (BD_jobCase) {
		case 0:
			PC_bsf_ProcessResults(
				&BD_exit, &BD_extendedReduceResult_P->elem, 
				BD_extendedReduceResult_P->reduceCounter,
				&BD_parameter, 
				&BD_newJobCase);
#ifdef PP_BSF_ITER_OUTPUT
			if (BD_iterCounter % PP_BSF_TRACE_COUNT == 0)
				PC_bsf_IterOutput(
					&BD_extendedReduceResult_P->elem, 
					BD_extendedReduceResult_P->reduceCounter, 
					BD_parameter,
					BD_t + MPI_Wtime(),
					BD_newJobCase);
#endif // PP_BSF_ITER_OUTPUT
			break;
		case 1:
			PC_bsf_ProcessResults_1(&BD_exit, 
				&BD_extendedReduceResult_P_1->elem, 
				BD_extendedReduceResult_P_1->reduceCounter,
				&BD_parameter, 
				&BD_newJobCase);
#ifdef PP_BSF_ITER_OUTPUT
			if (BD_iterCounter % PP_BSF_TRACE_COUNT == 0)
				PC_bsf_IterOutput_1(
					&BD_extendedReduceResult_P_1->elem, 
					BD_extendedReduceResult_P_1->reduceCounter, 
					BD_parameter,
					BD_t + MPI_Wtime(),
					BD_newJobCase);
#endif // PP_BSF_ITER_OUTPUT
			break;
		case 2:
			PC_bsf_ProcessResults_2(
				&BD_exit, 
				&BD_extendedReduceResult_P_2->elem, 
				BD_extendedReduceResult_P_2->reduceCounter,
				&BD_parameter, 
				&BD_newJobCase);
#ifdef PP_BSF_ITER_OUTPUT
			if (BD_iterCounter % PP_BSF_TRACE_COUNT == 0)
				PC_bsf_IterOutput_2(
					&BD_extendedReduceResult_P_2->elem, 
					BD_extendedReduceResult_P_2->reduceCounter, 
					BD_parameter,
					BD_t + MPI_Wtime(),
					BD_newJobCase);
#endif // PP_BSF_ITER_OUTPUT
			break;
		default:
			cout << "BC_Master: Undefined job type!" << endl;
			break;
		};
		BD_iterCounter++;
	} while (!BD_exit);
	BD_t += MPI_Wtime();

	BC_MasterMap(BD_EXIT);

	switch (BD_jobCase) {
	case 0:
		PC_bsf_ProblemOutput(&BD_extendedReduceResult_P->elem, BD_extendedReduceResult_P->reduceCounter, BD_parameter, BD_t);
		break;
	case 1:
		PC_bsf_ProblemOutput_1(&BD_extendedReduceResult_P_1->elem, BD_extendedReduceResult_P_1->reduceCounter, BD_parameter, BD_t);
		break;
	case 2:
		PC_bsf_ProblemOutput_2(&BD_extendedReduceResult_P_2->elem, BD_extendedReduceResult_P_2->reduceCounter, BD_parameter, BD_t);
		break;
	default:
		cout << "BC_Master: Undefined job type!" << endl;
		break;
	};
};
static void BC_Worker() {// Worker Process
	bool exit;
	while (true) {
		exit = BC_WorkerMap();
		/*debug*//*
#define RANK 0 
		if (BD_rank == RANK) {
			cout << RANK << ":BC_Worke: exit = " << exit << endl;
		};/*end debug*/
		if (exit) break;
		BC_WorkerReduce();
	};
};

static void BC_MasterMap(bool exit) {
	BV_AssignJobCase(BD_jobCase);
	BV_AssignIterCounter(BD_iterCounter);

	//*debug*/ cout << "BC_MasterMap: Job Case = " << BD_jobCase << endl;
	for (int rank = 0; rank < BD_numOfWorkers; rank++) {
		PC_bsf_CopyParameter(&BD_parameter, &(BD_order[rank].parameter));
		BD_order[rank].exit = exit;
		BD_order[rank].jobCase = BD_jobCase;
		BD_order[rank].iterCounter = BD_iterCounter;
		MPI_Isend(
			&BD_order[rank],
			sizeof(BT_order_T),
			MPI_BYTE,
			rank,
			0,
			MPI_COMM_WORLD,
			&BD_request[rank]);
	};
	MPI_Waitall(BD_numOfWorkers, BD_request, BD_status);
};

static void BC_MasterReduce() {
	//*debug*/ cout << "BC_MasterReduce: Job Case = " << BD_jobCase << endl;

	for (int rank = 0; rank < BD_numOfWorkers; rank++)
		switch (BD_jobCase) {
		case 0:
			MPI_Irecv(&BD_extendedReduceList[rank], sizeof(BT_extendedReduceElem_T), MPI_BYTE, rank, 0, MPI_COMM_WORLD, &BD_request[rank]);
			break;
		case 1:
			MPI_Irecv(&BD_extendedReduceList_1[rank], sizeof(BT_extendedReduceElem_T_1), MPI_BYTE, rank, 0, MPI_COMM_WORLD, &BD_request[rank]);
			break;
		case 2:
			MPI_Irecv(&BD_extendedReduceList_2[rank], sizeof(BT_extendedReduceElem_T_2), MPI_BYTE, rank, 0, MPI_COMM_WORLD, &BD_request[rank]);
			break;
		default:
			cout << "BC_MasterReduce Error: Undefined job type!" << endl;
			break;
		};

	MPI_Waitall(BD_numOfWorkers, BD_request, BD_status);

	switch (BD_jobCase) {
	case 0:
		BC_ProcessExtendedReduceList(BD_extendedReduceList, 0, BD_numOfWorkers, &BD_extendedReduceResult_P);
		break;
	case 1:
		BC_ProcessExtendedReduceList_1(BD_extendedReduceList_1, 0, BD_numOfWorkers, &BD_extendedReduceResult_P_1);
		break;
	case 2:
		BC_ProcessExtendedReduceList_2(BD_extendedReduceList_2, 0, BD_numOfWorkers, &BD_extendedReduceResult_P_2);
		break;
	default:
		cout << "BC_MasterReduce Error: Undefined job type!" << endl;
		break;
	};
};

static bool BC_WorkerMap() {
	/*debug*//*
#define RANK 0 
	if (BD_rank == RANK) {
		cout << RANK << ":BC_WorkerMap:";
	};/*end debug*/

	MPI_Recv(
		&BD_order[BD_rank],
		sizeof(BT_order_T),
		MPI_BYTE,
		BD_masterRank,
		0,
		MPI_COMM_WORLD,
		&BD_status[BD_rank]);

	if (BD_order[BD_rank].exit)
		return BD_EXIT;
	/*debug*//* if (BD_rank == RANK) {
		cout << " Job Case = " << BD_order[BD_rank].jobCase << endl;
	};/*end debug*/

	BV_AssignJobCase(BD_order[BD_rank].jobCase);
	BV_AssignIterCounter(BD_order[BD_rank].iterCounter);
	BV_AssignSublistLength(BD_sublistSize[BD_rank]);
	BV_AssignAddressOffset(BD_offset[BD_rank]);
	PC_bsf_SetParameter(&BD_order[BD_rank].parameter);

	/*debug*//*
	if (RANK == 0) cout << RANK << ":BC_WorkerMap: sublistSize = " 
		<< BD_sublistSize[BD_rank] << "\toffset = " << BD_offset[BD_rank] << endl;/*end debug*/

#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
#pragma omp parallel for num_threads(PP_BSF_NUM_THREADS)
#else
#pragma omp parallel for
#endif // PP_BSF_NUM_THREADS
#endif // PP_BSF_OMP
	for (int index = BD_offset[BD_rank]; index < BD_offset[BD_rank] + BD_sublistSize[BD_rank]; index++) {
		int mapIndex;
#ifdef PP_BSF_FRAGMENTED_MAP_LIST
		mapIndex = index - BD_offset[BD_rank];
#else
		mapIndex = index;
#endif
		BV_AssignNumberInSublist(mapIndex);
		switch (BD_order[BD_rank].jobCase) {
		case 0:
			BD_extendedReduceList[index].reduceCounter = 1;
			PC_bsf_MapF(&BD_mapSubList[mapIndex], &BD_extendedReduceList[index].elem,
				&BD_extendedReduceList[index].reduceCounter);
			break;
		case 1:
			BD_extendedReduceList_1[index].reduceCounter = 1;
			PC_bsf_MapF_1(&BD_mapSubList[mapIndex], &BD_extendedReduceList_1[index].elem,
				&BD_extendedReduceList_1[index].reduceCounter);
			break;
		case 2:
			BD_extendedReduceList_2[index].reduceCounter = 1;
			PC_bsf_MapF_2(&BD_mapSubList[mapIndex], &BD_extendedReduceList_2[index].elem,
				&BD_extendedReduceList_2[index].reduceCounter);
			break;
		default:
			cout << "BC_WorkerReduce Error: Undefined job type!" << endl;
			break;
		};
	};
	return !BD_EXIT;
};

static void BC_WorkerReduce() {
	/*debug*//*
#define RANK 0 
	if (BD_rank == RANK) {
		cout << RANK << ":BC_WorkerReduce: Job Case = " << BD_order[BD_rank].jobCase << endl;
	};/*end debug*/

	switch (BD_order[BD_rank].jobCase) {
	case 0:
		BC_ProcessExtendedReduceList(BD_extendedReduceList, BD_offset[BD_rank], BD_sublistSize[BD_rank],
			&BD_extendedReduceResult_P);
		MPI_Send(BD_extendedReduceResult_P, sizeof(BT_extendedReduceElem_T), MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD);
		break;
	case 1:
		BC_ProcessExtendedReduceList_1(BD_extendedReduceList_1, BD_offset[BD_rank], BD_sublistSize[BD_rank],
			&BD_extendedReduceResult_P_1);
		MPI_Send(BD_extendedReduceResult_P_1, sizeof(BT_extendedReduceElem_T_1), MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD);
		break;
	case 2:
		BC_ProcessExtendedReduceList_2(BD_extendedReduceList_2, BD_offset[BD_rank], BD_sublistSize[BD_rank],
			&BD_extendedReduceResult_P_2);
		MPI_Send(BD_extendedReduceResult_P_2, sizeof(BT_extendedReduceElem_T_2), MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD);
		break;
	default:
		cout << "BC_WorkerReduce Error: Undefined job type!" << endl;
		break;
	};
};

static void BC_ProcessExtendedReduceList(BT_extendedReduceElem_T* reduceList, int index, int length,
	BT_extendedReduceElem_T** extendedReduceResult_P) {
	int firstSuccessIndex = -1;

	*extendedReduceResult_P = &reduceList[index];

	for (int i = index; i < index + length; i++) {
		if (reduceList[i].reduceCounter > 0) {
			*extendedReduceResult_P = &reduceList[i];
			firstSuccessIndex = i;
			break;
		};
	};

	if (firstSuccessIndex >= 0) {
		for (int i = firstSuccessIndex + 1; i < index + length; i++)
			if (BD_extendedReduceList[i].reduceCounter > 0) {
				PC_bsf_ReduceF(&(*extendedReduceResult_P)->elem, &BD_extendedReduceList[i].elem,
					&(*extendedReduceResult_P)->elem);
				(*extendedReduceResult_P)->reduceCounter += BD_extendedReduceList[i].reduceCounter;
			};
	};
};

static void BC_ProcessExtendedReduceList_1(BT_extendedReduceElem_T_1* reduceList, int index, int length,
	BT_extendedReduceElem_T_1** extendedReduceResult_P) {
	int firstSuccessIndex = -1;

	*extendedReduceResult_P = &reduceList[index];

	for (int i = index; i < index + length; i++) {
		if (reduceList[i].reduceCounter > 0) {
			*extendedReduceResult_P = &reduceList[i];
			firstSuccessIndex = i;
			break;
		};
	};

	if (firstSuccessIndex >= 0) {
		for (int i = firstSuccessIndex + 1; i < index + length; i++)
			if (BD_extendedReduceList_1[i].reduceCounter > 0) {
				PC_bsf_ReduceF_1(&(*extendedReduceResult_P)->elem, &BD_extendedReduceList_1[i].elem,
					&(*extendedReduceResult_P)->elem);
				(*extendedReduceResult_P)->reduceCounter += BD_extendedReduceList_1[i].reduceCounter;
			};
	};
};

static void BC_ProcessExtendedReduceList_2(BT_extendedReduceElem_T_2* reduceList, int index, int length,
	BT_extendedReduceElem_T_2** extendedReduceResult_P) {
	int firstSuccessIndex = -1;

	*extendedReduceResult_P = &reduceList[index];

	for (int i = index; i < index + length; i++) {
		if (reduceList[i].reduceCounter > 0) {
			*extendedReduceResult_P = &reduceList[i];
			firstSuccessIndex = i;
			break;
		};
	};

	if (firstSuccessIndex >= 0) {
		for (int i = firstSuccessIndex + 1; i < index + length; i++)
			if (BD_extendedReduceList_2[i].reduceCounter > 0) {
				PC_bsf_ReduceF_2(&(*extendedReduceResult_P)->elem, &BD_extendedReduceList_2[i].elem,
					&(*extendedReduceResult_P)->elem);
				(*extendedReduceResult_P)->reduceCounter += BD_extendedReduceList_2[i].reduceCounter;
			};
	};
};

static void BC_Init(bool* success) {// Initialization
	//* debug */int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); if (rank == 0) *success = false;
	cout << setprecision(PP_BSF_PRECISION);
	PC_bsf_AssignListSize(&BD_listSize);

	BD_extendedReduceList = (BT_extendedReduceElem_T*)malloc(BD_listSize * sizeof(BT_extendedReduceElem_T));
	if (BD_extendedReduceList == NULL) {*success = false;return;};
	
	if (PP_BSF_MAX_JOB_CASE > 0) {
		BD_extendedReduceList_1 = (BT_extendedReduceElem_T_1*)malloc(BD_listSize * sizeof(BT_extendedReduceElem_T_1));
		if (BD_extendedReduceList_1 == NULL) { *success = false; return; };
	};

	if (PP_BSF_MAX_JOB_CASE > 1) {
		BD_extendedReduceList_2 = (BT_extendedReduceElem_T_2*)malloc(BD_listSize * sizeof(BT_extendedReduceElem_T_2));
		if (BD_extendedReduceList_2 == NULL) { *success = false; return; };
	};

	if (BD_size > BD_listSize + 1) {
		if (BD_rank == 0) cout << "Error: MPI_SIZE must be < Map List Size + 2 =" << BD_listSize + 2 << endl;
		MPI_Finalize();
		exit(1);
	};
	
	BD_masterRank = BD_size - 1;
	BD_numOfWorkers = BD_size - 1;
	BD_elemsPerWorker = BD_listSize / BD_numOfWorkers;
	BD_tailLength = BD_listSize - BD_elemsPerWorker * BD_numOfWorkers;
	BD_status = (MPI_Status*)malloc(BD_numOfWorkers * sizeof(MPI_Status));
	BD_request = (MPI_Request*)malloc(BD_numOfWorkers * sizeof(MPI_Request));
	BD_order = (BT_order_T*)malloc(BD_numOfWorkers * sizeof(BT_order_T));
	BD_sublistSize = (int*)malloc(BD_numOfWorkers * sizeof(int));
	BD_offset = (int*)malloc(BD_numOfWorkers * sizeof(int));

	PC_bsf_SetInitParameter(&BD_parameter);

	int offset = 0;
	for (int rank = 0; rank < BD_numOfWorkers; rank++) {
		BD_sublistSize[rank] = BD_elemsPerWorker + (rank < BD_tailLength ? 1 : 0);
		BD_offset[rank] = offset;
		offset += BD_sublistSize[rank];
	};
	if (BD_rank != BD_masterRank) {
#ifdef PP_BSF_FRAGMENTED_MAP_LIST
		BD_mapSubList = (PT_bsf_mapElem_T*)malloc(BD_sublistSize[BD_rank] * sizeof(PT_bsf_mapElem_T*));
#else
		BD_mapSubList = (PT_bsf_mapElem_T*)malloc(BD_listSize * sizeof(PT_bsf_mapElem_T*));
#endif
		if (BD_mapSubList == NULL) {
			*success = false;
			return;
		};

		BV_AssignMapSubList(BD_mapSubList);

#ifdef PP_BSF_FRAGMENTED_MAP_LIST
		PC_bsf_SetMapSubList(BD_mapSubList, BD_sublistSize[BD_rank], BD_offset[BD_rank], success);
#else
		PC_bsf_SetMapSubList(BD_mapSubList, BD_listSize, 0, success);
#endif
	};
};

static void BC_MpiRun() {
	int rc;
	rc = MPI_Init(NULL, NULL);	// Starting MPI
	if (rc != MPI_SUCCESS) {
		printf("Error starting MPI program. Terminating! \n");
		MPI_Abort(MPI_COMM_WORLD, rc);
	};
	MPI_Comm_rank(MPI_COMM_WORLD, &BD_rank);
	BV_AssignMpiRank(BD_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &BD_size);
	BV_AssignNumOfWorkers(BD_size - 1);
	if (BD_size < 2) {
		if (BD_rank == 0) cout << "Error: MPI_SIZE must be > 1" << endl;
		MPI_Finalize();
		exit(1);
	};
};