CC=gcc
FLAGS_OMP=-fopenmp -O3 -lm
FLAGS_SEQ=-O3 -lm
G=

bench:
	$(CC) -c src/common/bots_common.c -o obj/bots_common.o
	$(CC) -DBENCH_MODE=SEQ -I src/c/ -I src/serial/ -c src/common/bots_main.c -o obj/bots_main_seq.o
	$(CC) -DBENCH_MODE=OPENMP -I src/c/ -I src/omp/ -c src/common/bots_main.c -o obj/bots_main_omp.o -fopenmp
	$(CC) -c src/c/bench.c -o obj/bench.o
	$(CC) -c src/c/bench.c -o obj/bench_omp.o -fopenmp

build: bench
	$(CC) -I src/common/ -I src/c/ -c src/serial/sort.c -o obj/sort_seq.o
	$(CC) obj/sort_seq.o obj/bots_common.o obj/bots_main_seq.o obj/bench.o $(FLAGS_SEQ) -o bin/sort_seq

	$(CC) -I src/common/ -I src/c/ -c src/omp/sort.c -o obj/sort_omp.o $(FLAGS_OMP)
	$(CC) obj/sort_omp.o obj/bots_common.o obj/bots_main_omp.o obj/bench_omp.o $(FLAGS_OMP) -o bin/sort_omp



run:
	echo "{ \"out\" : ["
	./bin/sort_seq -v 0 -n 3000000
	echo ','
	./bin/sort_omp -v 0 -n 3000000
	echo "]}"

