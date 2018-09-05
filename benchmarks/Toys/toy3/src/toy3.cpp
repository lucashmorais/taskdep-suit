#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>

#include "../../common/Utils.h"

int N = 100;
int M = 2000;

void product(int* a, int* b, int* c) {
	for (int i=0; i<M; i++) {
		*c += (a[i] * b[i]);
	}
}

int main() {
	int matrix[N][M];
	int result[N][N];


	srand(100);
	for(int i=0;i<N;i++) {
		for (int j=0; j<N; j++)
			result[i][j] = 0;

		for (int j=0; j<M; j++)
			matrix[i][j] = rand() % (N*M);
	}

	double t_start, t_end;
            printf("Starting parallel code...\n");
            t_start = rtclock();

	#pragma omp parallel
	#pragma omp single
	{
		int a, b, c;

		for(int i=0; i<N; i++) {
			for (int j=0; j<N; j++) {
				int *a = matrix[i];
				int *b = matrix[j];
				int *c = &result[i][j];
				#pragma omp task depend(in:a, b) depend(out:c)
				product(a, b, c);
			}
		}
	}

	t_end = rtclock();
//
//	for(int i=0; i<N; i++) {
//		for (int j=0; j<N; j++) {
//
//			for (int k=0; k<M; k++) { printf("%02d ", matrix[i][k]); } printf("\n");
//			for (int k=0; k<M; k++) { printf("%02d ", matrix[j][k]); } printf("\n");
//			printf("------------------> %04d\n", result[i][j]);
//
//		}
//	}
//
	printf("Finishing.\n");
	int i, j;
	FILE *file;

#ifdef _OPENMP
	file = fopen("../../output/par", "w");
#else
	file = fopen("../../output/seq", "w");
#endif
	    for(i=0;i<N;i++) {
	    	for(j=0;i<N;i++) {
	          fprintf(file, "%d ", result[i][j]);
	      }
	      fprintf(file, "\n");
	    }
	    fclose(file);

	fprintf(stdout, "Time: %0.6lf\n", t_end - t_start);
	return 0;
}
