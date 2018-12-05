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
	$(CC) -I src/common/ -I src/c/ -c src/serial/fib.c -o obj/fib_seq.o
	$(CC) obj/fib_seq.o obj/bots_common.o obj/bots_main_seq.o obj/bench.o $(FLAGS_SEQ) -o bin/fib_seq

	$(CC) -I src/common/ -I src/c/ -c src/omp/fib.c -o obj/fib_omp.o $(FLAGS_OMP)
	$(CC) obj/fib_omp.o obj/bots_common.o obj/bots_main_omp.o obj/bench_omp.o $(FLAGS_OMP) -o bin/fib_omp

	$(CC) -DIF_CUTOFF -DBENCH_MODE=OPENMP -I src/c/ -I src/omp/ -c src/common/bots_main.c -o obj/bots_main_omp_coff.o -fopenmp
	$(CC) -DIF_CUTOFF -I src/common/ -I src/c/ -c src/omp/fib.c -o obj/fib_omp_coff.o $(FLAGS_OMP)
	$(CC) obj/fib_omp_coff.o obj/bots_common.o obj/bots_main_omp_coff.o obj/bench_omp.o $(FLAGS_OMP) -o bin/fib_omp_coff


	$(CC) -DFINAL_CUTOFF -DBENCH_MODE=OPENMP -I src/c/ -I src/omp/ -c src/common/bots_main.c -o obj/bots_main_omp_fin.o -fopenmp
	$(CC) -DFINAL_CUTOFF -I src/common/ -I src/c/ -c src/omp/fib.c -o obj/fib_omp_fin.o $(FLAGS_OMP)
	$(CC) obj/fib_omp_fin.o obj/bots_common.o obj/bots_main_omp_fin.o obj/bench_omp.o $(FLAGS_OMP) -o bin/fib_omp_fin


	$(CC) -DMANUAL_CUTOFF -DBENCH_MODE=OPENMP -I src/c/ -I src/omp/ -c src/common/bots_main.c -o obj/bots_main_omp_man.o -fopenmp
	$(CC) -DMANUAL_CUTOFF -I src/common/ -I src/c/ -c src/omp/fib.c -o obj/fib_omp_man.o $(FLAGS_OMP)
	$(CC) obj/fib_omp_man.o obj/bots_common.o obj/bots_main_omp_man.o obj/bench_omp.o $(FLAGS_OMP) -o bin/fib_omp_man



run:
	echo "{ \"out\" : ["
	./bin/fib_seq -v 0 -n 30
	echo ','
	./bin/fib_omp -v 0 -n 30
	echo ','
	./bin/fib_omp_coff -v 0 -n 30
	echo ','
	./bin/fib_omp_fin -v 0 -n 30
	echo ','
	./bin/fib_omp_man -v 0 -n 30
	echo "]}"

