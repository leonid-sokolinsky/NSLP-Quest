/*==============================================================================
Project: NSLP (Non-Stationary Linear Programming)
Theme: Quest Phase
Module: Problem-bsfCode.cpp (Implementation of the Problem)
Prefix: PI
Author(s): Leonid B. Sokolinsky, Irina M. Sokolinskaya
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-bsfParameters.h"	// Predefined Problem Parameters
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "BSF-VariableAccess.h"
using namespace std;

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	for (int j = 0; j < PP_N; j++) // Generating coordinates of initial appriximation
		parameter->approximation[j] = PP_INIT_APPROX;
//	parameter->approximation[0] = PP_INIT_APPROX;
};

void PC_bsf_Init(bool* success) {

	// Parameters for Waying
	PD_direction = (float*)malloc(PP_N * sizeof(float));
	if (PD_direction == NULL) { 
		*success = false;
		return;
	};

	for (int j = 0; j < PP_N; j++)  // Direction of moving
		PD_direction[j] = 1;
};

void PC_bsf_AssignListSize(int* listSize) {
	*listSize = PP_M;
};

void PC_bsf_CopyParameter(PT_bsf_parameter_T* parameterIn, PT_bsf_parameter_T* parameterOut) {
	for (int j = 0; j < PP_N; j++) {
		parameterOut->approximation[j] = parameterIn->approximation[j];
		parameterOut->shift[j] = parameterIn->shift[j];
	};
};

void PC_bsf_SetMapSubList(PT_bsf_mapElem_T* sublist, int sublistLength, int offset, bool* success) {
	//*debug*/int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); cout << rank << ":=========>PC_bsf_SetMapSubList 1: sublistLength = " << sublistLength << "\toffset = " << offset << endl;
	for (int i = 0; i < sublistLength; i++) {
		sublist[i].inequality = (float*)calloc(PP_N + 1,sizeof(float));
		if (sublist[i].inequality == NULL) {
			*success = false;
			return;
		};
		for (int j = 0; j < PP_N + 1; j++)
			sublist[i].inequality[j] = 0;
	};
	for (int i = offset; i < PF_MIN(PP_N, offset + sublistLength); i++) {
		sublist[i - offset].inequality[i] = 1;
		sublist[i - offset].inequality[PP_N] = PP_SF;
	};
	if ((offset <= PP_N) && (offset + sublistLength > PP_N)) {
		for (int j = 0; j < PP_N; j++)
			sublist[PP_N - offset].inequality[j] = 1;
		sublist[PP_N - offset].inequality[PP_N] = PP_SF * (PP_N - 1) + PP_SF / 2;
	};
	if ((offset <= PP_N + 1) && (offset + sublistLength > PP_N + 1)) {
		for (int j = 0; j < PP_N; j++)
			sublist[PP_N + 1 - offset].inequality[j] = -1;
		sublist[PP_N + 1 - offset].inequality[PP_N] = -PP_SF / 2;
	};
	for (int i = PF_MAX(PP_N + 2, offset); i < offset + sublistLength; i++) 
		sublist[i - offset].inequality[i - PP_N - 2] = -1;
};

void PC_bsf_SetParameter(PT_bsf_parameter_T* parameter) {
	PD_shift = parameter->shift;
	PD_approximation = parameter->approximation;
};

void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success // 1 - reduceElem was produced successfully; 0 - otherwise
){
	float factor;
	float wayed_b;

	wayed_b = Way_b(mapElem->inequality, PD_shift);

	factor = PF_MIN(wayed_b - DotProduct(PD_approximation, mapElem->inequality), 0) / NormSquare(mapElem->inequality);

	if (factor == 0)
		*success = 0;

	for (int j = 0; j < PP_N; j++)
		reduceElem->point[j] = factor * mapElem->inequality[j];

	reduceElem->pointIn = PointIn(PD_approximation, mapElem->inequality, wayed_b);

	//*debug*/cout << "PC_bsf_MapF: pointIn=" << reduceElem->pointIn << endl;
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

void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	for (int j = 0; j < PP_N; j++)
		z->point[j] = x->point[j] + y->point[j];
	z->pointIn = x->pointIn && y->pointIn;
};
void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {/* not used */};
void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {/* not used */};

void PC_bsf_ProcessResults(
	bool* exit, // "true" if Stopping Criterion is satisfied, and "false" otherwise
	PT_bsf_reduceElem_T* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* jobCase
) {
	for (int j = 0; j < PP_N; j++) PD_previousShift[j] = parameter->shift[j];
	GetShift(parameter->shift);
	//*debug*/cout << "PC_bsf_ProcessResults: pointIn=" << reduceResult->pointIn << endl;
	if (ExitCondition(reduceResult, reduceCounter, parameter))
		*exit = true;
	else {
		*exit = false;

		// New method
		float norm = 0;
		for (int j = 0; j < PP_N; j++)
			norm += reduceResult->point[j] * reduceResult->point[j];
		norm = sqrtf(norm);

		for (int j = 0; j < PP_N; j++)
			parameter->approximation[j] += (float)(PP_LAMBDA * reduceResult->point[j] / norm);
			//parameter->approximation[j] += (float)(PP_LAMBDA * reduceResult->point[j] / reduceCounter);
	};
	//cout << "Number of successfully produced Elrments of Reduce List: " << reduceCounter << endl;
};

void PC_bsf_ProcessResults_1(
	bool* exit, // "true" if Stopping Criterion is satisfied, and "false" otherwise
	PT_bsf_reduceElem_T_1* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* jobCase
) {
	// optional filling
};

void PC_bsf_ProcessResults_2(
	bool* exit, // "true" if Stopping Criterion is satisfied, and "false" otherwise
	PT_bsf_reduceElem_T_2* reduceResult,
	int reduceCounter, // Number of successfully produced Elrments of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* jobCase
) {
	// optional filling
};

void PC_bsf_ParametersOutput(int numOfWorkers, PT_bsf_parameter_T parameter) {
	cout << "=================================================== Quest Modified Fejer ====================================================" << endl;
	cout << "Number of Workers: " << numOfWorkers << endl;
	cout << "Velosity: " << PP_VELOCITY << endl;
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
	cout << "Eps_In = " << PP_EPS_IN << endl;
	cout << "Eps_Lag = " << PP_EPS_LAG << endl;

#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix A & Column b -------" << endl;
	for (int i = 0; i < PP_N; i++) {
		for (int j = 0; j < i; j++)
			cout << setw(5) << 0;
		cout << setw(5) << 1;
		for (int j = i + 1; j < PP_N; j++)
			cout << setw(5) << 0;
		cout << setw(5) << PP_SF << endl;
	};
	for (int j = 0; j < PP_N; j++)
		cout << setw(5) << 1;
	cout << setw(5) << PP_SF * (PP_N - 1) + PP_SF / 2 << endl;
	for (int j = 0; j < PP_N; j++)
		cout << setw(5) << -1;
	cout << setw(5) << -PP_SF / 2 << endl;
	for (int i = 0; i < PP_N; i++) {
		for (int j = 0; j < i; j++)
			cout << setw(5) << 0;
		cout << setw(5) << -1;
		for (int j = i + 1; j < PP_N; j++)
			cout << setw(5) << 0;
		cout << setw(5) << 0 << endl;
	};
#endif // PP_MATRIX_OUTPUT
	cout << "Start Point: "; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(7) << parameter.approximation[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;
	cout << "Direction: "; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(7) << PD_direction[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;
	cout << "-------------------------------------------" << endl;
};

void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int newJobCase) {
	cout << "-------------------- " << PP_BSF_iterCounter << " -------------------" << endl;
	//cout << "Reduce Result:"; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(10) << reduceResult->point [j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;/**/
	cout << "Shift:\t\t\t"; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(12) << parameter.shift[j]; cout << (PP_OUTPUT_LIMIT < PP_N - 1 ? "..." : "") << endl;/**/
	cout << "Approximation:\t\t"; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(12) << parameter.approximation[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;/**/
	cout << "Shift - Approximation :\t"; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(12) << parameter.shift[j] - parameter.approximation[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;
	cout << "PointIn = " << setw(1) << reduceResult->pointIn << endl;
};

void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int newJobCase)
{
	cout << "------------------ " << PP_BSF_iterCounter << " ------------------" << endl;
	/* not used */
};

void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int newJobCase)
{
	cout << "------------------ " << PP_BSF_iterCounter << " ------------------" << endl;
	/* not used */
};

void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {// Output Function
	cout << "=============================================" << endl;
	cout << "Time: " << t << endl;
	cout << "Iterations: " << PP_BSF_iterCounter << endl;
	cout << "Solution: "; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(12) << parameter.approximation[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;/**/
	cout << "Solution without Shift: "; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(12) << parameter.approximation[j] - PD_previousShift[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;
	/*float coordinate_sum = 0;
	for (int j = 0; j < PP_N; j++) 
		coordinate_sum += parameter.approximation[j] - PD_previousShift[j];
	cout << "Sum of coordindtes: " << coordinate_sum << endl;/**/
};

void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {// Output Function
	// optional filling
};

void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double t) {// Output Function
	// optional filling
};

//----------------------------- User functions -----------------------------
static bool PointIn(PT_point_T point, PT_inequality_T inequality, float wayed_b) { // If the point satisfies to the inequality 
	float sum = 0;

/*#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
#pragma omp parallel for num_threads(PP_BSF_NUM_THREADS) reduction(+:sum)
#else
#pragma omp parallel for reduction(+:sum)
#endif // PP_BSF_NUM_THREADS
#endif // PP_BSF_OMP/**/
	for (int j = 0; j < PP_N; j++)
		sum += inequality[j] * point[j];
	if (sum > wayed_b + PP_EPS_IN)
		return false;
	else 
		return true;
};

static float DotProduct(PT_point_T x, PT_point_T y) {
	float sum = 0;
/*#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
#pragma omp parallel for num_threads(PP_BSF_NUM_THREADS) reduction(+:sum)
#else
#pragma omp parallel for reduction(+:sum)
#endif // PP_BSF_NUM_THREADS
#endif // PP_BSF_OMP/**/
	for (int j = 0; j < PP_N; j++)
		sum += x[j] * y[j];
	return sum;
};

static float NormSquare(PT_point_T x) {
	float sum = 0;
/*#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
#pragma omp parallel for num_threads(PP_BSF_NUM_THREADS) reduction(+:sum)
#else
#pragma omp parallel for reduction(+:sum)
#endif // PP_BSF_NUM_THREADS
#endif // PP_BSF_OMP/**/
	for (int j = 0; j < PP_N; j++)
		sum += x[j] * x[j];
	return sum;
};

static float Way_b(PT_inequality_T inequality, PT_shift_T shift) {
	float sum = inequality[PP_N];

	for (int j = 0; j < PP_N; j++)
		sum += inequality[j] * shift[j];
	return sum;
};

static bool ExitCondition(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter) {
	static float shift_approx_prev = FLT_MAX, shift_approx_next;

	if (reduceResult->pointIn) {
		//*debug*/cout << "pointIn=" << reduceResult->pointIn << endl;
		return true;
	};

	shift_approx_next = parameter->shift[0] - parameter->approximation[0];
	if (PP_BSF_iterCounter > 100)
		if (shift_approx_prev + PP_EPS_LAG < shift_approx_next) {
			cout << ">>>>>>>>>>>>>>>>> Process began to lag!!! <<<<<<<<<<<<<<<<<<<<<<" << endl;
			cout << "shift_approx_next - shift_approx_prev = " << shift_approx_next - shift_approx_prev << endl;
			return true;
		};
	shift_approx_prev = shift_approx_next;

#ifdef PP_MAX_ITER_COUNT
	if (PP_BSF_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = " << PP_MAX_ITER_COUNT << endl;
		return true;
	};
#endif // PP_MAX_ITER_COUNT

	float sum = 0;
/*#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
#pragma omp parallel for num_threads(PP_BSF_NUM_THREADS) reduction(+:sum)
#else
#pragma omp parallel for reduction(+:sum)
#endif // PP_BSF_NUM_THREADS
#endif // PP_BSF_OMP/**/
	for (int j = 0; j < PP_N; j++) {
		float s;
		s = (float)(PP_LAMBDA * reduceResult->point[j] / reduceCounter);
		sum += s * s;
	};

	if (sum < PP_EPS_RELAX) 
		return true;
	
	return false;
};

static void GetShift(PT_shift_T shift) {
	static double t0 = MPI_Wtime();
	double delta_t;

	delta_t = MPI_Wtime() - t0;

	for (int j = 0; j < PP_N; j++) 
		shift[j] = (float)(PP_VELOCITY * PD_direction[j] * delta_t);/**/
		/*debug*///shift[j] = 0;
};