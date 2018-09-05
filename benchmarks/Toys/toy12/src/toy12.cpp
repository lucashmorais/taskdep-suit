/**
 *
 * The ideia of this program is that .each task in the task graph will depend
 * on two previous nodes.  Of course with exception of the initials.
 *
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>

#include "../../common/Utils.h"

void two_r__two_w(int *a, int *b, int *c, int *d) {
	*c = *a;
	*d = *b;
	//printf("[two_r__two_w]   R = [a = %d, b = %d]   W = [c = %d, d = %d]\n", *a, *b, *c, *d);
}

int main() {
	const int N = 100000;
	int* u = (int*) malloc(sizeof(int) * (2*N));
	srand(100);
	for(int i=0;i<2*N;i++) { u[i] = 111;}
	double t_start, t_end;
            printf("Starting parallel code...\n");
            t_start = rtclock();

	#pragma omp parallel
	#pragma omp single
	{
		for(int i=2; i<2*N; i+=2) {
                                        int *a = &u[i-2];
                                        int *b = &u[i-1];
                                        int *c = &u[i];
                                        int *d = &u[i+1];
			#pragma omp task depend(in:a, b) depend(out:c, d)
			two_r__two_w(a, b, c, d);
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
        for(s=0;s<(2*N);s++) {
              fprintf(file, "%d ", u[s]);
        }
        fclose(file);
        fprintf(stdout, "Time: %0.6lf\n", t_end - t_start);
	return 0;
}
