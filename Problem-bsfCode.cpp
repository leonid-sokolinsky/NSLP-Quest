/*==============================================================================
Project: CoFePro
Theme: Projection Algorithm for Solving Convex Feasibility Problems
Module: Problem-bsfCode.cpp (Implementation of the Problem)
Prefix: PI
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	for (int j = 0; j < PP_N; j++) // Generating initial approximation
		parameter->x[j] = PD_exteriorPoint[j];
};

void PC_bsf_Init(bool* success) {

	// ------------- Load inequality system -------------------
	PD_inequalitiesFile = PP_PATH;
	PD_inequalitiesFile += PP_LPP_FILE;
	const char* inequalitiesFile = PD_inequalitiesFile.c_str();

	FILE* stream;
	float buf;
	int m, n;
	stream = fopen(inequalitiesFile, "r");

	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << inequalitiesFile << "'.\n";
		*success = false; return;
	}

	if (fscanf(stream, "%d%d", &m, &n) == 0) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Unexpected end of file" << endl;
		*success = false;
		return;
	}

	if (n != PP_N || m != PP_M) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Error in input data '" << inequalitiesFile << "': PP_N != n and/or PP_M != m (PP_N = "
			<< PP_N << ", n = " << n << "; PP_M = " << PP_M << ", m = " << m << ").\n";
		*success = false; return;
	}

	for (int i = 0; i < PP_M; i++) {
		for (int j = 0; j < PP_N; j++) {
			if (fscanf(stream, "%f", &buf) == 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster) cout << "Unexpected end of file" << endl; *success = false; return;
			};
			PD_A[i][j] = buf;
		}
		if (fscanf(stream, "%f", &buf) == 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster) cout << "Unexpected end of file" << endl; *success = false; return;
		};
		PD_b[i] = buf;
	}
	fclose(stream);

	// --------------- Load exterior point ---------------
	PD_exteriorPointFile = PP_PATH;
	PD_exteriorPointFile += PP_EXTERIOR_POINT_FILE;
	const char* exteriorPointFile = PD_exteriorPointFile.c_str();
	stream = fopen(exteriorPointFile, "r");
	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << exteriorPointFile << "'.\n";
		*success = false; return;
	}

	if (fscanf(stream, "%d", &n) == 0) { cout << "Unexpected end of file" << endl; *success = false; return; }
	if (n != PP_N) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Error in input data '" << exteriorPointFile << "': PP_N != n (PP_N = " << PP_N << ", n = " << n << ").\n";
		*success = false; return;
	}

	for (int j = 0; j < PP_N; j++) {
		if (fscanf(stream, "%f", &buf) == 0) { if (BSF_sv_mpiRank == BSF_sv_mpiMaster) cout << "Unexpected end of file" << endl; *success = false; return; }
		PD_exteriorPoint[j] = buf;
	}

	bool pointIn = true;
	for (int i = 0; i < PP_M; i++) 
		if (!PointIn(PD_exteriorPoint, PD_A[i], PD_b[i])) {
			pointIn = false;
			break;
		}

	if (pointIn) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "The point in the file '" 
			<< PD_exteriorPointFile << "' is the interior point of the inequality system from the file '"
			<< PD_inequalitiesFile << "'.\n";
		*success = false; return;
	}

	fclose(stream);

	*success = true;

	cout << "The calculations have started, please wait..." << endl;

}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = PP_M;
}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	elem->a = PD_A[i];
	elem->b = &(PD_b[i]);

	// Calculating norm square
	elem->normSquare = Vector_NormSquare(PD_A[i]);
}

void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success // 1 - reduceElem was produced successfully; 0 - otherwise
){
	*success = Vector_ProjectOnHalfspace(BSF_sv_parameter.x, mapElem->a, *mapElem->b, reduceElem->projection);
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

void PC_bsf_JobDispatcher(
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* job,
	bool* exit
) {
	// Optional filling. Do not delete!
}

void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int jobCase) {
	cout << "------------------ " << BSF_sv_iterCounter << " ------------------" << endl;
	// optional filling

};

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	cout << "=================================================== Quest ====================================================" << endl;
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

#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix A & Column b -------" << endl;
	for (int i = 0; i < PP_M; i++) {
		cout << i << ")";
		for (int j = 0; j < PP_N; j++)
			cout << setw(PP_SETW) << PD_A[i][j];
		cout << "\t<=" << setw(PP_SETW) << PD_b[i] << endl;
	}
#endif // PP_MATRIX_OUTPUT
	cout << "Exterior Point: "; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(12) << PD_exteriorPoint[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;
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
	cout << "Elapsed time: " << t << endl;
	cout << "Iterations: " << BSF_sv_iterCounter << endl;
	cout << "Interior point: "; for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PP_N); j++) cout << setw(12) << parameter.x[j]; cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;

	PD_interiorPointFile = PP_PATH;
	PD_interiorPointFile += PP_INTERIOR_POINT_FILE;
	if (SaveSolution(parameter.x, PD_interiorPointFile))
		cout << "Solution is saved into the file '" << PD_interiorPointFile << "'." << endl;
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
void PC_bsfAssignMpiMaster(int value) { BSF_sv_mpiMaster = value; };
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; };
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; };
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; };
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; };

//----------------------- Problem functions ---------------------------
static double Vector_DotProduct(PT_vector_T x, PT_vector_T y) {
	double sum = 0;
	for (int j = 0; j < PP_N; j++)
		sum += x[j] * y[j];
	return sum;
};

static double Vector_NormSquare(PT_vector_T x) {
	double sum = 0;
	for (int j = 0; j < PP_N; j++)
		sum += x[j] * x[j];
	return sum;
};

// Point projection onto Half-space <a,x> <= b
inline bool // true if the point does not belong to the half-space and false otherwise 
Vector_ProjectOnHalfspace(PT_vector_T point, PT_vector_T a, PT_float_T b, PT_vector_T projection) {
	double factor;
	double aNormSquare = Vector_NormSquare(a);

	if (aNormSquare < PP_EPS_ZERO)
		return false;

	factor = (b - Vector_DotProduct(point, a)) / aNormSquare;

	if (factor > PP_EPS_ZERO)
		return false;

	for (int j = 0; j < PP_N; j++) {
		projection[j] = factor * a[j];
	}

	return true;
}

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

static bool SaveSolution(PT_vector_T x, string solutionFile) {
	FILE* stream;

	const char* file = solutionFile.c_str();
	stream = fopen(file, "w");
	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << solutionFile << "'.\n";
		return false;
	}

	fclose(stream);
	return true;
}

inline bool PointIn(PT_vector_T x, PT_vector_T a, PT_float_T b) { // If the point belonges to the Halfspace <a,x> <= b
	PT_float_T dotProduct_a_x = Vector_DotProduct(a, x);

	if (dotProduct_a_x < b + PP_EPS_ZERO)
		return true;
	else
		return false;
}