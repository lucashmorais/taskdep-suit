CC=gcc
FLAGS_OMP=-fopenmp -O3 -lm
FLAGS_SEQ=-O3 -lm
G=

bench:
	$(CC) $(G) -c src/bench.c -o obj/bench.o
	$(CC) $(G) -c src/bench.c -o obj/bench_omp.o -fopenmp
	$(CC) $(G) -c src/main.c -o obj/main.o
	$(CC) $(G) -c src/poisson.c -o obj/poisson.o
	$(CC) $(G) -c $(FLAGS_OMP) src/main.c -o obj/main-omp.o
	$(CC) $(G) -c $(FLAGS_OMP) src/poisson.c -o obj/poisson-omp.o
	


build: bench
	$(CC) $(G) src/jacobi-seq.c obj/main.o obj/bench.o obj/poisson.o $(FLAGS_SEQ) -o bin/jacobi-seq
	$(CC) $(G) src/jacobi-task.c src/jacobi-seq.c obj/main-omp.o obj/bench_omp.o obj/poisson-omp.o $(FLAGS_OMP) -o bin/jacobi-task
	$(CC) $(G) src/jacobi-task-dep.c src/jacobi-seq.c obj/main-omp.o obj/bench_omp.o obj/poisson-omp.o $(FLAGS_OMP) -o bin/jacobi-task-dep

run:
	echo "{ \"out\" : ["
	./bin/jacobi-seq -i 1 -n 8192 -b 128 -r 8
	echo ','
	./bin/jacobi-task -i 16 -n 8192 -b 128 -r 8
	echo ','
	./bin/jacobi-task-dep -i 16 -n 8192 -b 128 -r 8
	echo "]}"
