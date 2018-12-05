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
	$(CC) -I src/common/ -I src/c/ -c src/serial/fft.c -o obj/fft_seq.o
	$(CC) obj/fft_seq.o obj/bots_common.o obj/bots_main_seq.o obj/bench.o $(FLAGS_SEQ) -o bin/fft_seq

	$(CC) -I src/common/ -I src/c/ -c src/omp/fft.c -o obj/fft_omp.o $(FLAGS_OMP)
	$(CC) obj/fft_omp.o obj/bots_common.o obj/bots_main_omp.o obj/bench_omp.o $(FLAGS_OMP) -o bin/fft_omp

run:
	echo "{ \"out\" : ["
	./bin/fft_seq -v 0
	echo ','
	./bin/fft_omp -v 0
	echo "]}"

