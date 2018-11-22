/**
 *
 * The ideia of this program the task graph will be like a linked list, that is,
 * each task will read from its predecessor.
 *
 */
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>

#include "../../common/Utils.h"
void one_read_write(int *a, int *b) {
	*b = *a;
	///printf("[one_read_write] a = %d, b = %d\n", *a, *b);
}

int main() {
	const int N = 100000;
	int* u = (int*) malloc(sizeof(int) * N);
	for(int i=0;i<N;i++) { u[i] = 111;}
            double t_start, t_end;
            printf("Starting parallel code...\n");
            t_start = rtclock();
	#pragma omp parallel
	#pragma omp single
	{
		for(int i=1; i<N; i++) {
                                        int *a = &u[i-1];
                                        int *b = &u[i];
			#pragma omp task depend(in:a) depend(out:b)
			one_read_write(a, b);
		}
	}
            t_end = rtclock();
	srand(100);
	printf("Finished. %d\n", u[rand() % N]);
            int s;

    FILE *file;

#ifdef _OPENMP
    file = fopen("../../output/par", "w");
#else
    file = fopen("../../output/seq", "w");
#endif
        for(s=0;s<N;s++) {
              fprintf(file, "%d ", u[s]);
        }
        fclose(file);

        fprintf(stdout, "Time: %0.6lf\n", t_end - t_start);
	return 0;
}
