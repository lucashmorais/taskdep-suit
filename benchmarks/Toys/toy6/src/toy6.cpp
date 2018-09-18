/**
 *
 * The ideia of this program is that it will create a task graph where there is
 * no dependence between the tasks because each task read from *two* different
 * address.
 *
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>

#include "../../common/Utils.h"
void two_read(int *a, int *b) {
	*a = *b;
///	printf("[two_read] a = %d, b = %d\n", *a, *b);
}

int main() {
	const int N = 100000;
	int* u = (int *) malloc(sizeof(int) * (2*N + N));

	for(int i=1; i<=2*N; i++) { u[i] = 111;}

	double t_start, t_end;
            printf("Starting parallel code...\n");
            t_start = rtclock();

	#pragma omp parallel
	#pragma omp single
	{
		for(int i=1; i<=2*N; i+=2) {
                                        int *a = &u[i];
                                        int *b = &u[i+1];
			#pragma omp task depend(in:a, b)
			two_read(a, b);
		}
	}
	t_end = rtclock();

	srand(100);
	printf("Finished. %d %d\n", u[rand() %N], u[rand() % N]);
	int s;

	for(s=0;s<(2*N + N);s++) {
			fprintf(stdout, "%d ", u[s]);
	}
	fprintf(stdout, "Time: %0.6lf\n", t_end - t_start);
	fprintf(stderr, "Time: %0.6lf\n", t_end - t_start);
	return 0;
}
