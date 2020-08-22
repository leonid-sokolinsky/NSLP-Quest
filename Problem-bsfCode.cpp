/*==============================================================================
Project: NSLP (Non-Stationary Linear Programming)
Theme: Quest Phase
Module: Problem-bsfCode.cpp (Implementation of the Problem)
Prefix: PI
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	for (int j = 0; j < PP_N; j++) // Generating initial approximation
		parameter->x[j] = PD_apex[j];
};

void PC_bsf_Init(bool* success) {

	// Generating Objective Function Coefficients
	for (int j = 0; j < PP_N; j++)
		PD_c[j] = PP_N - j;

	// Generating Coordinates of Apex Point
	double c_normSquare = NormSquare(PD_c);
	for (int j = 0; j < PP_N; j++) 
		PD_apex[j] = PP_DIST_TO_APEX * PD_c[j] / sqrt(c_normSquare);
	
	// Generating A & b
	for (int i = 0; i < PP_N; i++) {
		for (int j = 0; j < PP_N; j++)
			PD_A[i][j] = 0;
		PD_A[i][i] = 1;
		PD_b[i] = PP_SF;
	};
	for (int j = 0; j < PP_N; j++)
		PD_A[PP_N][j] = 1;
	PD_b[PP_N] = PP_SF * (PP_N - 1) + PP_SF / 2;
	for (int j = 0; j < PP_N; j++)
		PD_A[PP_N + 1][j] = -1;
	PD_b[PP_N + 1] = -PP_SF / 2;
	for (int i = PP_N + 2; i < PP_M; i++) {
		for (int j = 0; j < PP_N; j++)
			PD_A[i][j] = 0;
		PD_A[i][i - PP_N - 2] = -1;
		PD_b[i] = 0;
	};

	// Calculating a_i norm squares
	for (int i = 0; i < PP_M; i++) {
		PD_normSquare_a[i] = 0;
		for (int j = 0; j < PP_N; j++)
			PD_normSquare_a[i] += PD_A[i][j] * PD_A[i][j];
	};
	
	/* debug *//* if (PP_BSF_mpiRank == 0) {
		cout << "----------PC_bsf_Init-------------" << endl;
		for (int i = 0; i < PP_M; i++)
			cout << "InequalityNo = " << i << "\tNorm Square = " << PD_normSquare_a[i] << endl;
		//system("pause");
	};/* end debug */

};

void PC_bsf_SetListSize(int* listSize) {
	*listSize = PP_M;
};

void PC_bsf_SetMapSubList(PT_bsf_mapElem_T* sublist, int sublistLength, int offset) {
	//*debug*/int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); cout << rank << ":=========>PC_bsf_SetMapSubList: sublistLength = " << sublistLength << "\toffset = " << offset << endl; cout.flush();
	for (int i = 0; i < sublistLength; i++)
			sublist[i].inequalityNo = offset + i;
};

void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success // 1 - reduceElem was produced successfully; 0 - otherwise
){
	double factor;

	factor = (PD_b[mapElem->inequalityNo] - DotProduct(BSF_sv_parameter.x, PD_A[mapElem->inequalityNo])) / PD_normSquare_a[mapElem->inequalityNo];

	if (factor > 0)
		*success = false;
	else
		for (int j = 0; j < PP_N; j++)
			reduceElem->projection[j] = factor * PD_A[mapElem->inequalityNo][j];

	/* debug *//* if (PP_BSF_mpiRank == 0) {
		cout << "Hyperplane No = " << mapElem->inequalityNo << "\tProjection: ";
		if (factor >= 0)
			cout << "\tsuccess = false" << endl;
		else {
			cout << "\tProjection Vector: ";
			for (int j = 0; j < PP_N; j++)
				cout << setw(12) << reduceElem->projection[j];
		};
		cout << endl;
		system("pause");
	};/* end debug */
};

void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem,
	int* success // 1 - reduceElem was produced successfully (default); 0 - otherwise
) {
	/* not used */
};

void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem,
	int* success // 1 - reduceElem was produced successfully (default); 0 - otherwise
) {
	/* not used */
};

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem,
	int* success // 1 - reduceElem was produced successfully (default); 0 - otherwise
	) {
	// optional filling
};

void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	for (int j = 0; j < PP_N; j++)
		z->projection[j] = x->projection[j] + y->projection[j];
};
void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {/* not used */};
void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {/* not used */};
void PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {/* not used */ }

void PC_bsf_ProcessResults(
	PT_bsf_reduceElem_T* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	if (ExitCondition(reduceResult, reduceCounter, parameter))
		*exit = true;
	else {
		*exit = false;

		/*// New method
		double norm = 0;
		for (int j = 0; j < PP_N; j++)
			norm += reduceResult->projection[j] * reduceResult->projection[j];
		norm = sqrt(norm);/**/

		for (int j = 0; j < PP_N; j++)
			//parameter->x[j] += (PP_LAMBDA * reduceResult->projection[j] / norm);
			parameter->x[j] += PP_LAMBDA * reduceResult->projection[j] / reduceCounter;
	};
	//cout << "Number of successfully produced Elrments of Reduce List: " << reduceCounter << endl;
};

void PC_bsf_ProcessResults_1(
	PT_bsf_reduceElem_T_1* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// optional filling
};

void PC_bsf_ProcessResults_2(
	PT_bsf_reduceElem_T_2* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// optional filling
};

void PC_bsf_ProcessResults_3(
	PT_bsf_reduceElem_T_3* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// optional filling
};

void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int jobCase) {
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	// optional filling

};

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	cout << "=================================================== Quest Modified Fejer ====================================================" << endl;
	cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;
#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
	cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
#else
	cout << "Number of Threads: " << omp_get_num_procs() << endl;
#endif // PP_BSF_NUM_THREADS
#else
	cout << "OpenMP is turned off!" << endl;
#endif // PP_BSF_OMP
	cout << "Dimension: N = " << PP_N << endl;
	cout << "Number of Constraints: M = " << PP_M << endl;
	cout << "Scale Factor: SF = " << PP_SF << endl;
	cout << "Relaxation Factor: LAMBDA = " << PP_LAMBDA << endl;
	cout << "Eps_Relax = " << PP_EPS_RELAX << endl;
	cout << "Distance to Apex = " << PP_DIST_TO_APEX << endl;

#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix A & Column b -------" << endl;
	for (int i = 0; i < PP_M; i++) {
		for (int j = 0; j < PP_N; j++)
			cout << setw(5) << PD_A[i][j];
		cout << "\t<=" << setw(5) << PD_b[i] << endl;
	};
#endif // PP_MATRIX_OUTPUT
	cout << "Objective Function: "; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(2) << PD_c[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;
	cout << "Apex Point: "; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(12) << PD_apex[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;
	cout << "-------------------------------------------" << endl;
	//* debug */ system("pause");/* end debug */
};

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	for (int i = 0; i < PP_N; i++)
		parameterOutP->x[i] = parameterIn.x[i];
};

void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob) {
	cout << "-------------------- " << BSF_sv_iterCounter << " -------------------" << endl;
	//cout << "Reduce Result:"; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(10) << reduceResult->point [j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;/**/
	cout << "Approximation:\t\t"; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(12) << parameter.x[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;/**/
};

void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	/* not used */
};

void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	/* not used */
};

void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {// Output Function
	cout << "=============================================" << endl;
	cout << "Time: " << t << endl;
	cout << "Iterations: " << BSF_sv_iterCounter << endl;
	cout << "Solution: "; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(12) << parameter.x[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;/**/
};

void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {// Output Function
	// optional filling
};

void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {// Output Function
	// optional filling
};

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {// Output Function
	// optional filling
};

//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; };
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; };
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; };
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; };
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; };
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; };
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; };

//----------------------- Problem functions ---------------------------
static double DotProduct(PT_point_T x, PT_point_T y) {
	double sum = 0;
	for (int j = 0; j < PP_N; j++)
		sum += x[j] * y[j];
	return sum;
};

static double NormSquare(PT_point_T x) {
	double sum = 0;
	for (int j = 0; j < PP_N; j++)
		sum += x[j] * x[j];
	return sum;
};

static bool ExitCondition(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter) {
	static double shift_approx_prev = FLT_MAX, shift_approx_next;

#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = " << PP_MAX_ITER_COUNT << endl;
		return true;
	};
#endif // PP_MAX_ITER_COUNT

	double sum = 0;
	for (int j = 0; j < PP_N; j++) {
		double s;
		s = (PP_LAMBDA * reduceResult->projection[j] / reduceCounter);
		sum += s * s;
	};

	if (sum < PP_EPS_RELAX) 
		return true;
	
	return false;
};