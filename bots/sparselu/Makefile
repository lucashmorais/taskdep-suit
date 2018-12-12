CC=gcc
FLAGS_OMP=-fopenmp -O3 -lm
FLAGS_SEQ=-O3 -lm
G=

bench:
	$(CC) -c src/common/bots_common.c -o obj/bots_common.o
	$(CC) -DBENCH_MODE=SEQ -I src/c/ -I src/serial/ -c src/common/bots_main.c -o obj/bots_main_seq.o
	$(CC) -DBENCH_MODE=OPENMP -I src/c/ -I src/omp/sparselu_for/ -c src/common/bots_main.c -o obj/bots_main_omp.o -fopenmp
	$(CC) -c src/c/bench.c -o obj/bench.o
	$(CC) -c src/c/bench.c -o obj/bench_omp.o -fopenmp

build: bench
	$(CC) -I src/common/ -I src/c/ -c src/serial/sparselu.c -o obj/sparselu_seq.o
	$(CC) obj/sparselu_seq.o obj/bots_common.o obj/bots_main_seq.o obj/bench.o $(FLAGS_SEQ) -o bin/sparselu_seq

	$(CC) -I src/common/ -I src/c/ -c src/omp/sparselu_single/sparselu.c -o obj/sparselu_omp_single.o $(FLAGS_OMP)
	$(CC) obj/sparselu_omp_single.o obj/bots_common.o obj/bots_main_omp.o obj/bench_omp.o $(FLAGS_OMP) -o bin/sparselu_omp_single

	$(CC) -I src/common/ -I src/c/ -c src/omp/sparselu_for/sparselu.c -o obj/sparselu_omp_for.o $(FLAGS_OMP)
	$(CC) obj/sparselu_omp_for.o obj/bots_common.o obj/bots_main_omp.o obj/bench_omp.o $(FLAGS_OMP) -o bin/sparselu_omp_for


run:
	echo "{ \"out\" : ["
	./bin/sparselu_seq -v 0
	echo ','
	./bin/sparselu_omp_single -v 0
	echo ','
	./bin/sparselu_omp_for -v 0
	echo "]}"

