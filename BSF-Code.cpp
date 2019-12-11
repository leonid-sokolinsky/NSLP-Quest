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
			cout << "Error: PC_bsf_Init failed (not enough memory)!" << endl;
		MPI_Finalize();
		exit(1);
	};

	BD_success = true; 
	BC_Init(&BD_success);
	MPI_Allreduce(&BD_success, &success, 1, MPI_UNSIGNED, MPI_LAND, MPI_COMM_WORLD);
	if (!success) {
		if (BD_rank == BD_masterRank) cout << "Error: BC_Init failed (not enough memory). N = " << endl;
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
	PC_bsf_ParametersOutput(BD_numOfWorkers, BD_data);
	BD_iterCounter = 0;

	BD_t -= MPI_Wtime();
	do {
		BC_MasterMap(!BD_EXIT);
		BC_MasterReduce();
		PC_bsf_ProcessResults(&BD_exit, &BD_extendedReduceResult_P->elem, BD_extendedReduceResult_P->reduceCounter, &BD_data, BD_iterCounter);
		BD_iterCounter++;
#ifdef PP_BSF_ITER_OUTPUT
		if (BD_iterCounter % PP_BSF_TRACE_COUNT == 0)
			PC_bsf_IterOutput(&BD_extendedReduceResult_P->elem, BD_extendedReduceResult_P->reduceCounter, BD_data,
				BD_iterCounter, BD_t + MPI_Wtime());
#endif // PP_BSF_ITER_OUTPUT
	} while (!BD_exit);
	BD_t += MPI_Wtime();

	BC_MasterMap(BD_EXIT);

	PC_bsf_ProblemOutput(&BD_extendedReduceResult_P->elem, BD_extendedReduceResult_P->reduceCounter, BD_data, 
		BD_iterCounter, BD_t);
};
static void BC_Worker() {// Worker Process
	bool exit;
	while (true) {
		exit = BC_WorkerMap();
		if (exit) break;
		BC_WorkerReduce();
	};
};

static void BC_MasterMap(bool exit) {
	for (int rank = 0; rank < BD_numOfWorkers; rank++) {
		PC_bsf_CopyData(&BD_data, &(BD_order[rank].data));
		BD_order[rank].exit = exit;
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
	for (int rank = 0; rank < BD_numOfWorkers; rank++)
		MPI_Irecv(
			&BD_extendedReduceList[rank],
			sizeof(BT_extendedReduceElem_T),
			MPI_BYTE,
			rank,
			0,
			MPI_COMM_WORLD,
			&BD_request[rank]);
	MPI_Waitall(BD_numOfWorkers, BD_request, BD_status);
	BC_ProcessExtendedReduceList(BD_extendedReduceList, 0, BD_numOfWorkers, &BD_extendedReduceResult_P);
};

static bool BC_WorkerMap() {
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

#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
#pragma omp parallel for num_threads(PP_BSF_NUM_THREADS)
#else
#pragma omp parallel for
#endif // PP_BSF_NUM_THREADS
#endif // PP_BSF_OMP
	//*debug*/int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); if (rank == 3) cout << "BC_WorkerMap: sublistSize = " << BD_sublistSize[BD_rank] << "\toffset = " << BD_offset[BD_rank] << endl;
	for (int index = BD_offset[BD_rank]; index < BD_offset[BD_rank] + BD_sublistSize[BD_rank]; index++) {
		BD_extendedReduceList[index].reduceCounter = 1;
		PC_bsf_MapF(&BD_mapSubList[index - BD_offset[BD_rank]], &BD_extendedReduceList[index].elem, &BD_order[BD_rank].data,
			index - BD_offset[BD_rank], BD_sublistSize[BD_rank], BD_offset[BD_rank], &BD_extendedReduceList[index].reduceCounter);
	};
	return !BD_EXIT;
};

static void BC_WorkerReduce() {
	BC_ProcessExtendedReduceList(BD_extendedReduceList, BD_offset[BD_rank], BD_sublistSize[BD_rank], 
		&BD_extendedReduceResult_P);

	MPI_Send(
		BD_extendedReduceResult_P,
		sizeof(BT_extendedReduceElem_T),
		MPI_BYTE,
		BD_masterRank,
		0,
		MPI_COMM_WORLD);
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

static void BC_Init(bool* success) {// Initialization
	//* debug */int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); if (rank == 0) *success = false;
	cout << setprecision(PP_BSF_PRECISION);
	PC_bsf_AssignListSize(&BD_listSize);

	BD_extendedReduceList = (BT_extendedReduceElem_T*)malloc(BD_listSize * sizeof(BT_extendedReduceElem_T));
	
	if (BD_extendedReduceList == NULL) {
		*success = false;
		return;
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

	PC_bsf_SetInitData(&BD_data);

	int offset = 0;
	for (int rank = 0; rank < BD_numOfWorkers; rank++) {
		BD_sublistSize[rank] = BD_elemsPerWorker + (rank < BD_tailLength ? 1 : 0);
		BD_offset[rank] = offset;
		offset += BD_sublistSize[rank];
	};
	if (BD_rank != BD_masterRank) {
		BD_mapSubList = (PT_bsf_mapElem_T*)malloc(BD_sublistSize[BD_rank] * sizeof(PT_bsf_mapElem_T*));
		if (BD_mapSubList == NULL) {
			*success = false;
			return;
		};
	PC_bsf_SetMapSubList(BD_mapSubList, BD_sublistSize[BD_rank], BD_offset[BD_rank], success);
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
	MPI_Comm_size(MPI_COMM_WORLD, &BD_size);
	if (BD_size < 2) {
		if (BD_rank == 0) cout << "Error: MPI_SIZE must be > 1" << endl;
		MPI_Finalize();
		exit(1);
	};
};
