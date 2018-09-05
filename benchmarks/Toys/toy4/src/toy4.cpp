#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>

#include "../../common/Utils.h"

void fn1(int * a, int * b) {
	*b = *a;
//	printf("fun1.\n");
}

void fn2(int * a, int * b) {
	*b = *a;
//	printf("fun2.\n");
}

int main() {
    int i;
    int N = 100;
    int u[N*100], v[N*100];


    srand(100);
    for(i=0;i<N*100;i++) {
    	u[i] = 222;
    	v[i] = 222;
    }
    u[0] = -1;
    v[0] = 111;
    double t_start, t_end;
    printf("Starting parallel code...\n");
    t_start=rtclock();
    #pragma omp parallel
    #pragma omp single
	{
		int i, j, k;

		for(i=1, j=1;i<N;i++) {
			int *a = &v[i-1];
			int *b = &v[i];

			#pragma omp task depend(in:a) depend(out:b)
			fn1(a, b);

			for (k = 0;k<i;k++, j++) {
				int *c = &u[j];

				 #pragma omp task depend(in:b) depend(out:c)
				fn2(b, c);
			}
		}
	}
    t_end=rtclock();

//	printf("u: [");
//	for(i=0;i<N;i++)
//		printf("%d ", u[i]);
//	printf("]\n");
//
//	printf("v: [");
//	for(i=0;i<N;i++)
//		printf("%d ", v[i]);
//	printf("]\n");

        printf("Finishing.\n");

    int s;
    FILE *file;

#ifdef _OPENMP
    file = fopen("../../output/par", "w");
#else
    file = fopen("../../output/seq", "w");
#endif
        for(s=0;s<N*100;s++) {
              fprintf(file, "%d ", u[s]);
        }
        fprintf(file, "\n");
        for(s=0;s<N*100;s++) {
              fprintf(file, "%d ", v[s]);
        }
        fclose(file);
        fprintf(stdout, "Time: %0.6lf\n", t_end - t_start);

	return 0;
}
