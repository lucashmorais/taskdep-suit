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
	$(CC) -I src/common/ -I src/c/ -c src/serial/strassen.c -o obj/strassen_seq.o
	$(CC) obj/strassen_seq.o obj/bots_common.o obj/bots_main_seq.o obj/bench.o $(FLAGS_SEQ) -o bin/strassen_seq

	$(CC) -I src/common/ -I src/c/ -c src/omp/strassen.c -o obj/strassen_omp.o $(FLAGS_OMP)
	$(CC) obj/strassen_omp.o obj/bots_common.o obj/bots_main_omp.o obj/bench_omp.o $(FLAGS_OMP) -o bin/strassen_omp

	$(CC) -DIF_CUTOFF -DBENCH_MODE=OPENMP -I src/c/ -I src/omp/ -c src/common/bots_main.c -o obj/bots_main_omp_coff.o -fopenmp
	$(CC) -DIF_CUTOFF -I src/common/ -I src/c/ -c src/omp/strassen.c -o obj/strassen_omp_coff.o $(FLAGS_OMP)
	$(CC) obj/strassen_omp_coff.o obj/bots_common.o obj/bots_main_omp_coff.o obj/bench_omp.o $(FLAGS_OMP) -o bin/strassen_omp_coff


run:
	echo "{ \"out\" : ["
	./bin/strassen_seq -v 0 -n 2048
	echo ','
	./bin/strassen_omp -v 0 -n 2048
	echo ','
	./bin/strassen_omp_coff -v 0 -n 2048
	echo "]}"

