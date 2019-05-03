/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: BSF-MRA Skeleton
Module: BSF-Implementation.cpp (Implementation of BSF-MRA Skeleton)
Prefix: BI
Author(s): Leonid B. Sokolinsky
==============================================================================*/
#pragma once
#include "BSF-Include.h"
#include "Problem-Types.h"			// Problem Types 
#include "BSF-Types.h"				// BSF Types 
#include "BSF-Forwards.h"			// BSF Forwards
#include "BSF-Data.h"				// BSF Data 
#include "Problem-bsf-Forwards.h"		// Problem Forwards 

using namespace std;

int main(int argc, char *argv[]) {
	char emptystring[] = "";
	char* message = emptystring;
	unsigned success;

	MpiRun();
	
	BD_success = true;
	PI_bsf_Init(&BD_success);

	MPI_Allreduce(&BD_success, &success, 1, MPI_UNSIGNED, MPI_LAND, MPI_COMM_WORLD);

	if (!success) {
		if (BD_rank == BD_masterRank) 
			cout << "Error: PI_bsf_Init failed (not enough memory)!" << endl;
		MPI_Finalize();
		exit(1);
	};

	BD_success = true; 
	BI_Init(&BD_success);

	MPI_Allreduce(&BD_success, &success, 1, MPI_UNSIGNED, MPI_LAND, MPI_COMM_WORLD);

	if (!success) {
		if (BD_rank == BD_masterRank) cout << "Error: BI_Init failed (not enough memory). N = " << endl;
		MPI_Finalize();
		exit(1);
	};

	//BI_MeasureTimeParameters();

	if (BD_rank == BD_masterRank)
		BI_Master();
	else
		BI_Worker();
	
	MPI_Finalize();
	return 0;
};
static void BI_Master() {// Master Process
	PI_bsf_ParametersOutput(BD_numOfWorkers, BD_data);
	BD_iterCount = 0;

	BD_t -= MPI_Wtime();
	do {
		BI_MasterMap(!BD_EXIT);
		BI_MasterReduce();
		BD_t_p -= MPI_Wtime();
		PI_bsf_ProcessResults(&BD_exit, &BD_extendedReduceResult_P->elem, BD_extendedReduceResult_P->counter, &BD_data);
		BD_t_p += MPI_Wtime();
		BD_iterCount++;
#ifdef PP_BSF_ITER_OUTPUT
		if (BD_iterCount % PP_BSF_TRACE_COUNT == 0)
			PI_bsf_IterOutput(&BD_extendedReduceResult_P->elem, BD_extendedReduceResult_P->counter, BD_data,
				BD_iterCount, BD_t + MPI_Wtime());
#endif // PP_BSF_ITER_OUTPUT
	} while (!BD_exit);
	BD_t += MPI_Wtime();

	BI_MasterMap(BD_EXIT);

	BD_t_A_w /= BD_iterCount;
	BD_t_W /= BD_iterCount;
	BD_t_p /= BD_iterCount;
	BD_t_S /= BD_iterCount;

	PI_bsf_ProblemOutput(&BD_extendedReduceResult_P->elem, BD_extendedReduceResult_P->counter, BD_data, 
		BD_iterCount, BD_t, BD_t_L, BD_t_s_L, BD_t_S, BD_t_r_L, BD_t_W, BD_t_A_w, BD_t_A_m, BD_t_p);
};
static void BI_Worker() {// Worker Process
	bool exit;
	while (true) {
		exit = BI_WorkerMap();
		if (exit) break;
		BI_WorkerReduce();
	};
};

static void BI_MasterMap(bool exit) {
	BD_t_S -= MPI_Wtime();
	for (int rank = 0; rank < BD_numOfWorkers; rank++) {
		PI_bsf_CopyData(&BD_data, &(BD_order[rank].data));
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
	BD_t_S += MPI_Wtime();
};

static void BI_MasterReduce() {
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

	for (int rank = 0; rank < BD_numOfWorkers; rank++)
		BD_t_W += BD_extendedReduceList[rank].t_W;

	for (int rank = 0; rank < BD_numOfWorkers; rank++)
		BD_t_A_w += BD_extendedReduceList[rank].t_A_w;

	BD_t_A_m -= MPI_Wtime();
	BI_ProcessExtendedReduceList(BD_extendedReduceList, 0, BD_numOfWorkers, &BD_extendedReduceResult_P);
	BD_t_A_m += MPI_Wtime();
};

static bool BI_WorkerMap() {
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

	BD_t_W -= MPI_Wtime();

#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
#pragma omp parallel for num_threads(PP_BSF_NUM_THREADS)
#else
#pragma omp parallel for
#endif // PP_BSF_NUM_THREADS
#endif // PP_BSF_OMP/**/
	for (int index = BD_offset[BD_rank]; index < BD_offset[BD_rank] + BD_subListSize[BD_rank]; index++) {
		BD_extendedReduceList[index].counter = 1;
		PI_bsf_MapF(&BD_mapSubList[index - BD_offset[BD_rank]], &BD_extendedReduceList[index].elem, index, &BD_order[BD_rank].data,
			&BD_extendedReduceList[index].counter);
	};

	BD_t_W += MPI_Wtime();

	return !BD_EXIT;
};

static void BI_WorkerReduce() {
	BD_t_A_w -= MPI_Wtime();
	BI_ProcessExtendedReduceList(BD_extendedReduceList, BD_offset[BD_rank], BD_subListSize[BD_rank], 
		&BD_extendedReduceResult_P);

	BD_t_A_w += MPI_Wtime();

	BD_extendedReduceResult_P->t_W = BD_t_W;
	BD_extendedReduceResult_P->t_A_w = BD_t_A_w;

	MPI_Send(
		BD_extendedReduceResult_P,
		sizeof(BT_extendedReduceElem_T),
		MPI_BYTE,
		BD_masterRank,
		0,
		MPI_COMM_WORLD);
};

static void BI_ProcessExtendedReduceList(BT_extendedReduceElem_T* reduceList, int index, int length,
	BT_extendedReduceElem_T** extendedReduceResult_P) {
	int firstSuccessIndex = -1;

	*extendedReduceResult_P = &reduceList[index];

	for (int i = index; i < index + length; i++) {
		if (reduceList[i].counter > 0) {
			*extendedReduceResult_P = &reduceList[i];
			firstSuccessIndex = i;
			break;
		};
	};

	if (firstSuccessIndex >= 0) {
		for (int i = firstSuccessIndex + 1; i < index + length; i++)
			if (BD_extendedReduceList[i].counter > 0) {
				PI_bsf_ReduceF(&(*extendedReduceResult_P)->elem, &BD_extendedReduceList[i].elem, 
					&(*extendedReduceResult_P)->elem);
				(*extendedReduceResult_P)->counter += BD_extendedReduceList[i].counter;
			};
	};	
};

static void BI_Init(bool* success) {// Initialization
	//* debug */int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); if (rank == 0) *success = false;

	PI_bsf_AssignListSize(&BD_listSize);

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
	BD_subListSize = (int*)malloc(BD_numOfWorkers * sizeof(int));
	BD_offset = (int*)malloc(BD_numOfWorkers * sizeof(int));
	PI_bsf_SetInitApproximation(&BD_data);
	int offset = 0;
	for (int rank = 0; rank < BD_numOfWorkers; rank++) {
		BD_subListSize[rank] = BD_elemsPerWorker + (rank < BD_tailLength ? 1 : 0);
		BD_offset[rank] = offset;
		offset += BD_subListSize[rank];
	};
	if (BD_rank != BD_masterRank) {
		BD_mapSubList = (PT_bsf_mapElem_T*)malloc(BD_subListSize[BD_rank] * sizeof(PT_bsf_mapElem_T*));
		if (BD_mapSubList == NULL) {
			*success = false;
			return;
		};
	PI_bsf_SetMapSubList(BD_mapSubList, BD_subListSize[BD_rank], BD_offset[BD_rank], success);
	};
};

static void BI_MeasureTimeParameters() {
	bool dummyByte1;
	bool dummyByte2;
	BT_order_T dummyOrded;
//	BT_extendedReduceElem_T dummyExtendedReduceElem;
	
	if (BD_rank == BD_masterRank) {

		MPI_Barrier(MPI_COMM_WORLD); // Measuring latency time 
		BD_t_L = -MPI_Wtime();
		MPI_Send(&dummyByte1, sizeof(bool), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
		BD_t_L += MPI_Wtime();

		MPI_Barrier(MPI_COMM_WORLD); // Measuring send time 
		BD_t_s_L = -MPI_Wtime();
		MPI_Send(&dummyOrded, sizeof(BT_order_T), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
		BD_t_s_L += MPI_Wtime();

		/*MPI_Barrier(MPI_COMM_WORLD); // Measuring receive time
		BD_t_r_L = -MPI_Wtime();
		MPI_Recv(&dummyExtendedReduceElem, sizeof(BT_extendedReduceElem_T), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &BD_status[BD_rank]);
		BD_t_r_L += MPI_Wtime();/**/

		MPI_Barrier(MPI_COMM_WORLD); // Measuring latency time 
		BD_t_L = -MPI_Wtime();
		MPI_Send(&dummyByte2, sizeof(bool), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
		BD_t_L += MPI_Wtime();/**/

		MPI_Barrier(MPI_COMM_WORLD);
	} else {
		MPI_Barrier(MPI_COMM_WORLD);// Measuring latency time
		if (BD_rank == 0)
			MPI_Recv(&dummyByte1, sizeof(bool), MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD, &BD_status[BD_rank]);

		MPI_Barrier(MPI_COMM_WORLD);// Measurement of send time 
		if (BD_rank == 0)
			MPI_Recv(&dummyOrded, sizeof(BT_order_T), MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD, &BD_status[BD_rank]);

		/*MPI_Barrier(MPI_COMM_WORLD);// Measurement of receive time 
		if (BD_rank == 0)
			MPI_Send(&dummyExtendedReduceElem, sizeof(BT_extendedReduceElem_T), MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD);/**/
		
		MPI_Barrier(MPI_COMM_WORLD);// Measuring latency time
		if (BD_rank == 0)
			MPI_Recv(&dummyByte2, sizeof(bool), MPI_BYTE, BD_masterRank, 0, MPI_COMM_WORLD, &BD_status[BD_rank]);/**/

		MPI_Barrier(MPI_COMM_WORLD);
	};
};

static void MpiRun() {
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
