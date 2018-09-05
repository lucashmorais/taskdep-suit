/**
 *
 * The ideia of this program is that it will create a task graph where there is
 * no dependence between the tasks because each task write to *three* different
 * address.
 *
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>

#include "../../common/Utils.h"
void three_write(int *a, int *b, int *c) {
	*a = 10;
	*b = 20;
	*c = 30;
	///printf("[three_write] a = %d, b = %d, c = %d\n", *a, *b, *c);
}

int main() {
	const int N = 100000;
	int* u = (int *) malloc(sizeof(int) * (3*N + N));

	srand(100);
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
			#pragma omp task depend(out:a, b, c)
			three_write(a, b, c);
		}
	}
	t_end = rtclock();
	printf("Finished. %d %d %d\n", u[rand() % N], u[rand() % N], u[rand() % N]);
	int s;

    FILE *file;

#ifdef _OPENMP
    file = fopen("../../output/par", "w");
#else
    file = fopen("../../output/seq", "w");
#endif
        for(s=0;s<(3 * N + N);s++) {
              fprintf(file, "%d ", u[s]);
        }
        fclose(file);
        fprintf(stdout, "Time: %0.6lf\n", t_end - t_start);
	return 0;
}
