#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

using namespace std;

int numDP = 1000; // Vietoviu skaicius (demand points, max 10000)
int numPF = 10;	  // Esanciu objektu skaicius (preexisting facilities)
int numCL = 25;	  // Kandidatu naujiems objektams skaicius (candidate locations)
int numX = 3;	  // Nauju objektu skaicius

double **demandPoints; // Geografiniai duomenys

//=============================================================================

double getTime();
void loadDemandPoints();
double HaversineDistance(double *a, double *b);
double evaluateSolution(int *X, int tag);
int increaseX(int *X, int index, int maxindex);

//=============================================================================

int main()
{

	double ts = getTime(); // Algoritmo vykdymo pradzios laikas

	loadDemandPoints(); // Nuskaitomi duomenys

	// Sudarom pradini sprendini: [0, 1, 2, 3, ...]
	int *X = new int[numX];
	int *bestX = new int[numX];
	for (int i = 0; i < numX; i++)
	{
		X[i] = i;
		bestX[i] = i;
	}
	MPI_Init(NULL, NULL);
	int tag = 0;
	double u = evaluateSolution(X, tag++);
	double bestU = u;
	int r;

	//----- Pagrindinis ciklas ------------------------------------------------
	while (true)
	{
		if (increaseX(X, numX - 1, numCL))
		{
			u = evaluateSolution(X, tag++);
			if (u > bestU)
			{
				bestU = u;
				for (int i = 0; i < numX; i++)
					bestX[i] = X[i];
			}
		}
		else
			break;
	}
	//----- Rezultatu spausdinimas --------------------------------------------

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
	{

		double tf = getTime(); // Skaiciavimu pabaigos laikas

		cout << "Geriausias sprendinys: ";
		for (int i = 0; i < numX; i++)
			cout << bestX[i] << " ";
		cout << "(" << bestU << ")" << endl;
		cout << "Skaiciavimo trukme: " << tf - ts << endl;
	}

	MPI_Finalize();
}

//=============================================================================

void loadDemandPoints()
{
	FILE *f;
	f = fopen("demandPoints.dat", "r");
	demandPoints = new double *[numDP];
	for (int i = 0; i < numDP; i++)
	{
		demandPoints[i] = new double[3];
		fscanf(f, "%lf%lf%lf", &demandPoints[i][0], &demandPoints[i][1], &demandPoints[i][2]);
	}
	fclose(f);
}

//=============================================================================

double HaversineDistance(double *a, double *b)
{
	double dlon = fabs(a[0] - b[0]);
	double dlat = fabs(a[1] - b[1]);
	double aa = pow((sin((double)dlon / (double)2 * 0.01745)), 2) + cos(a[0] * 0.01745) * cos(b[0] * 0.01745) * pow((sin((double)dlat / (double)2 * 0.01745)), 2);
	double c = 2 * atan2(sqrt(aa), sqrt(1 - aa));
	double d = 6371 * c;
	return d;
}

//=============================================================================

double getTime()
{
	struct timeval laikas;
	gettimeofday(&laikas, NULL);
	double rez = (double)laikas.tv_sec + (double)laikas.tv_usec / 1000000;
	return rez;
}

//=============================================================================

double evaluateSolution(int *X, int tag)
{
	double U = 0;
	int bestPF;
	int bestX;
	double d;

	int rank;
	int size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int iters = numDP / size;
	int i = rank * iters;
	int iMax = i + iters;

	for (; i < iMax; i++)
	{
		bestPF = 1e5;
		for (int j = 0; j < numPF; j++)
		{
			d = HaversineDistance(demandPoints[i], demandPoints[j]);
			if (d < bestPF)
				bestPF = d;
		}
		bestX = 1e5;
		for (int j = 0; j < numX; j++)
		{
			d = HaversineDistance(demandPoints[i], demandPoints[X[j]]);
			if (d < bestX)
				bestX = d;
		}
		if (bestX < bestPF)
			U += demandPoints[i][2];
		else if (bestX == bestPF)
			U += 0.3 * demandPoints[i][2];
	}

	double UBuf[1];
	if (rank == 0)
	{
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(UBuf, 1, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			U += UBuf[0];
		}
		UBuf[0] = U;
	}
	else
	{
		UBuf[0] = U;
		MPI_Send(UBuf, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
	}
	MPI_Bcast(UBuf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return UBuf[0];
}

//=============================================================================

int increaseX(int *X, int index, int maxindex)
{
	if (X[index] + 1 < maxindex - (numX - index - 1))
	{
		X[index]++;
	}
	else
	{
		if ((index == 0) && (X[index] + 1 == maxindex - (numX - index - 1)))
		{
			return 0;
		}
		else
		{
			if (increaseX(X, index - 1, maxindex))
				X[index] = X[index - 1] + 1;
			else
				return 0;
		}
	}
	return 1;
}