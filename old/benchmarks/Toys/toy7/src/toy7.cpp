/**
 *
 * The ideia of this program is that it will create a task graph where there is
 * no dependence between the tasks because each task read from *three* different
 * address.
 *
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>

#include "../../common/Utils.h"
void three_read(int *a, int *b, int *c) {
	*a = *b = *c = 10;
///	printf("[three_read] a = %d, b = %d, c = %d\n", *a, *b, *c);
}

int main() {
	const int N = 100000;
	int *u = (int *) malloc(sizeof(int) * (3 * N + N));

	for(int i=1; i<=3*N; i++) { u[i] = 111;}
	double t_start, t_end;
            printf("Starting parallel code...\n");
            t_start = rtclock();
	#pragma omp parallel
	#pragma omp single
	{
		for(int i=1; i<=3*N; i+=3) {
			int *a = &u[i];
			int *b = &u[i+1];
			int *c = &u[i+2];

			#pragma omp task depend(in:a, b, c)
			three_read(a, b, c);
		}
	}
	t_end = rtclock();

	srand(100);
	printf("Finished. %d %d %d\n", u[rand() % N], u[rand() % N], u[rand() % N]);
	int s;
             FILE *file;

	for(s=0;s<(3 * N + N);s++) {
			fprintf(stdout, "%d ", u[s]);
	}
	fprintf(stdout, "Time: %0.6lf\n", t_end - t_start);
	fprintf(stderr, "Time: %0.6lf\n", t_end - t_start);
	return 0;
}
