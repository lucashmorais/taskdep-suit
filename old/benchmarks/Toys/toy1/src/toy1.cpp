#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>

#include "../../common/Utils.h"

void fn1(int * a, int * b) {
	*b = *a;
	// printf("Fun1-Task) a=%p b=%p\n", a, b);
}

void fn2(int * a, int * b) {
	*b = *a;
	// printf("Fun2-Task) a=%p b=%p\n", a, b);
}

int main() {
	int i;
	int N = 50;
	int u[N*100], v[N*100];

	//srand(time(NULL));
	srand(100);
	for(i=0;i<N*100;i++) {
		u[i] = rand();
		v[i] = rand();
	}
	u[0] = -1;
	v[0] = 666;

            double t_start, t_end;
            printf("Starting parallel code...\n");
            t_start = rtclock();
            int j, k;
	#pragma omp parallel
	#pragma omp single
	{


		for(i=1, j=1;i<N;i++) {
			int *a = &v[i-1];
			int *b = &v[i];
			#pragma omp task depend(in:a) depend(out:b)
			fn1(a, b);

			for(k = 0;k<i;k++, j++) {
				int *c = &u[j];
				#pragma omp task depend(in:b) depend(out:c)
				fn2(b, c);
			}

		}
	}
			#pragma omp barrier

            t_end = rtclock();
	printf("u: [");
	for(i=0;i<N;i++)
		printf("%d ", u[i]);
	printf("]\n");

	printf("v: [");
	for(i=0;i<N;i++)
		printf("%d ", v[i]);
	printf("]\n");

	printf("Finishing.\n");

	for(i=0;i<N*100;i++) {
			printf("%d ", u[i]);
	}

	printf("\n");

	for(i=0;i<N*100;i++) {
			printf("%d ", v[i]);
	}

	fprintf(stdout, "Time: %0.6lf\n", t_end - t_start);
	fprintf(stderr, "Time: %0.6lf\n", t_end - t_start);
	return 0;
}
