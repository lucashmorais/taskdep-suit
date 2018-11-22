#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX 10000

int * randSample(int size);
int avgSample(int * sample, int size);
int maxSample(int * sample, int size);
int minSample(int * sample, int size);


int main() {
	srand(time(NULL));

	int * sample, size = 100000000, avg, max, min;
	#pragma omp parallel
	#pragma omp single
	{
		#pragma omp task depend(out: sample)
		sample = randSample(size);

		#pragma omp task depend(in: sample)
		avg = avgSample(sample, size);
		#pragma omp task depend(in: sample)
		max = maxSample(sample, size);
		#pragma omp task depend(in: sample)
		min = minSample(sample, size);

		#pragma omp taskwait
		printf("avg is  %d, max is %d, min is %d\n", avg, max, min);
	}
	return 0;
}

int * randSample(int size) {
	int * t = malloc(sizeof(int) * size);
	for(int i = 0; i < size; i++) t[i] = (rand() % (MAX + 1));
	return t;
}

int avgSample(int * sample, int size) {
	long long int t = 0;
	long long int holder = 0;
	long long int s = 0;
	long long int q = 0;
	for(int i = 0; i < size; i++, holder++) {
		t += sample[i];
		if(holder > 1024) {
			holder = 0;
			s++;
			q += t/1024;
			t = 0;
		}
	}
	return (int) q/holder;
}


int maxSample(int * sample, int size) {
	int t = 0;
	for(int i = 0; i < size; i++) t = (t < sample[i]) ? sample[i] : t;
	return t;
}

int minSample(int * sample, int size) {
	int t = MAX;
	for(int i = 0; i < size; i++) t = (t > sample[i]) ? sample[i] : t;
	return t;
}

