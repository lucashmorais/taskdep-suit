CC=gcc
FLAGS_OMP=-fopenmp -O3 -lm
FLAGS_SEQ=-O3 -lm
G=

bench:
	$(CC) $(G) -c src/bench.c -o obj/bench.o
	$(CC) $(G) -c src/bench.c -o obj/bench-omp.o -fopenmp
	$(CC) $(G) -c src/main.c -o obj/main.o
	$(CC) $(G) -c src/sparselu.c -o obj/sparselu.o
	$(CC) $(G) -c $(FLAGS_OMP) src/main.c -o obj/main-omp.o
	$(CC) $(G) -c $(FLAGS_OMP) src/sparselu.c -o obj/sparselu-omp.o
	


build: bench
	$(CC) $(G) src/sparselu-seq.c obj/main.o obj/bench.o obj/sparselu.o $(FLAGS_SEQ) -o bin/sparselu-seq
	$(CC) $(G) src/sparselu-task.c src/sparselu-seq.c obj/main-omp.o obj/bench-omp.o obj/sparselu-omp.o $(FLAGS_OMP) -o bin/sparselu-task
	$(CC) $(G) src/sparselu-task-dep.c src/sparselu-seq.c obj/main-omp.o obj/bench-omp.o obj/sparselu-omp.o $(FLAGS_OMP) -o bin/sparselu-task-dep

run:
	echo "{ \"out\" : ["
	./bin/sparselu-seq -i 3 -n 128 -m 64
	echo ','
	./bin/sparselu-task -i 16 -n 128 -m 64
	echo ','
	./bin/sparselu-task-dep -i 16 -n 128 -m 64
	echo "]}"
