/**
 *
 * The ideia of this program is that it will create a task graph where there is
 * no dependence between the tasks because each task read from a different
 * address.
 *
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>

#include "../../common/Utils.h"
void one_read(int *a) {
	*a = 100;
	//printf("[one_read] a = %d\n", *a);
}

int main() {
	const int N = 100000;
	int u[N];


	for(int i=0;i<N;i++) { u[i] = 111;}

	double t_start, t_end;
            printf("Starting parallel code...\n");
            t_start = rtclock();
	#pragma omp parallel
	#pragma omp single
	{
		for(int i=0; i<N; i++) {
			int *a = &u[i];

			#pragma omp task depend(in:a)
			one_read(a);
		}
	}
	t_end = rtclock();
	printf("Finished. %d\n", u[rand() % N]);
	int s;
	for(s=0;s<N;s++) {
			fprintf(stdout, "%d ", u[s]);
	}
	fprintf(stdout, "Time: %0.6lf\n", t_end - t_start);
	fprintf(stderr, "Time: %0.6lf\n", t_end - t_start);
	return 0;
}
