/**
 *
 * The ideia of this program is that it will create a task graph where there is
 * no dependence between the tasks because each task write to a different
 * address.
 *
 */
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>

#include "../../common/Utils.h"
void one_write(int *a) {
	*a = 10;
	///printf("[one_write] a = %d\n", *a);
}

int main() {
	const int N = 100000;
	int* u = (int *) malloc(sizeof(int) * N);

	srand(time(NULL));

	for(int i=0;i<N;i++) { u[i] = 111; }
	double t_start, t_end;
            printf("Starting parallel code...\n");
            t_start = rtclock();
	#pragma omp parallel
	#pragma omp single
	{
		for(int i=0; i<N; i++) {
			int *a = &u[i];
			#pragma omp task depend(out:a)
			one_write(a);
		}
	}
	t_end = rtclock();
	printf("Finished. %d \n", u[rand() % N]);
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
