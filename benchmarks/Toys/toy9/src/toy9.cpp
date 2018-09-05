/**
 *
 * The ideia of this program is that it will create a task graph where there is
 * no dependence between the tasks because each task write to *two* different
 * address.
 *
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
 #include <cstdlib>

#include "../../common/Utils.h"

void two_write(int *a, int *b) {
	*a = 10;
	*b = 20;
	//printf("[two_write] a = %d, b = %d\n", *a, *b);
}

int main() {
	const int N = 100000;
	int* u = (int *) malloc(sizeof(int) * (2*N + N));
	srand(100);
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
			#pragma omp task depend(out:a, b)
			two_write(a, b);
		}
	}
	t_end = rtclock();
	printf("Finished. %d %d\n", u[rand() % N], u[rand() % N]);
	int s;
             FILE *file;

#ifdef _OPENMP
    file = fopen("../../output/par", "w");
#else
    file = fopen("../../output/seq", "w");
#endif
        for(s=0;s<(2*N + N);s++) {
              fprintf(file, "%d ", u[s]);
        }
        fclose(file);
        fprintf(stdout, "Time: %0.6lf\n", t_end - t_start);
	return 0;
}
