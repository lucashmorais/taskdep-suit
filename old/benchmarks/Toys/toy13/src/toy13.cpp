/**
 *
 * The ideia of this program is that .each task in the task graph will depend
 * on three previous nodes.  Of course with exception of the initials.
 *
 */

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>

#include "../../common/Utils.h"
void three_r__three_w(int *a, int *b, int *c, int *d, int *e, int *f) {
	*d = *a;
	*e = *b;
	*f = *c;
	///printf("[three_r__three_w]   R = [a = %d, b = %d, c = %d]   W = [d = %d, e = %d, f = %d]\n", *a, *b, *c, *d, *e, *f);
}

int main() {
	const int N = 100000;
	int* u = (int *) malloc(sizeof(int) * (3*N + N));

	for(int i=0;i<3*N;i++) { u[i] = 111;}

	double t_start, t_end;
            printf("Starting parallel code...\n");
            t_start = rtclock();
	#pragma omp parallel
	#pragma omp single
	{
		for(int i=3; i<3*N; i+=3) {
                                        int *a = &u[i-3];
                                        int *b = &u[i-2];
                                        int *c = &u[i-1];
                                        int *d = &u[i];
                                        int *e = &u[i+1];
                                        int *f = &u[i+2];
			#pragma omp task depend(in:a, b, c) depend(out:d, e, f)
			three_r__three_w(a, b, c, d, e, f);
		}
	}
	t_end = rtclock();
	srand( time(NULL) );
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
