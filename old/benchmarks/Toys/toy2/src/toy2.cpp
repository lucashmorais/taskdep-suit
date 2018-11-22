#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <cstdlib>

#include "../../common/Utils.h"
void minimum(int* a, int* b, int* c) {
	int tmp = 0;

	if (*a > *b) {
		tmp = *b;
		*b = *a;
		*a = tmp;
	}

	if (*a > *c) {
		tmp = *c;
		*c = *a;
		*a = tmp;
	}

	if (*b > *c) {
		tmp = *c;
		*c = *b;
		*b = tmp;
	}
}

int main() {
	int N = 10;
	int u[N], v[N], w[N];


	for(int i=0;i<N;i++) {
		u[i] = rand() % N;
		v[i] = rand() % N;
		w[i] = rand() % N;
	}

	 double t_start, t_end;
            printf("Starting parallel code...\n");
            t_start = rtclock();

	#pragma omp parallel
	#pragma omp single
	{
		int a, b, c;
		srand(100);
		for(int i=0; i<N*N; i++) {
			a = rand() % N;
			b = rand() % N;
			c = rand() % N;
			int *d = &u[a];
			int *e = &v[b];
			int *f = &w[c];

			#pragma omp task depend(inout:d, e, f)
			minimum(d, e, f);
		}
	}

	//for(int i=0;i<N;i++) {
	//	printf("%02d) %02d %02d %02d", i, u[i], v[i], w[i]);
	//	if (u[i] > v[i] || v[i] > w[i]) printf(" ***");
	//	printf("\n");
	//}
	t_end = rtclock();

	printf("Finishing.\n");
		int i;

	/*for(i=0;i<N;i++) {
			fprintf(file, "%d ", u[i]);
	}
	fprintf(file, "\n");
	for(i=0;i<N;i++) {
			fprintf(file, "%d ", v[i]);
	}
	fprintf(file, "\n");
	for(i=0;i<N;i++) {
			fprintf(file, "%d ", w[i]);
	}*/
	fwrite(u, sizeof(int), N, stdout);
	fwrite(v, sizeof(int), N, stdout);
	fwrite(w, sizeof(int), N, stdout);

	fprintf(stderr, "Time: %0.6lf\n", t_end - t_start);
	fprintf(stdout, "Time: %0.6lf\n", t_end - t_start);
	return 0;
}
