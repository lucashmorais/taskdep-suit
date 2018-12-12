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
	$(CC) -I src/common/ -I src/c/ -c src/serial/health.c -o obj/health_seq.o
	$(CC) obj/health_seq.o obj/bots_common.o obj/bots_main_seq.o obj/bench.o $(FLAGS_SEQ) -o bin/health_seq

	$(CC) -I src/common/ -I src/c/ -c src/omp/health.c -o obj/health_omp.o $(FLAGS_OMP)
	$(CC) obj/health_omp.o obj/bots_common.o obj/bots_main_omp.o obj/bench_omp.o $(FLAGS_OMP) -o bin/health_omp

	$(CC) -DIF_CUTOFF -DBENCH_MODE=OPENMP -I src/c/ -I src/omp/ -c src/common/bots_main.c -o obj/bots_main_omp_coff.o -fopenmp
	$(CC) -DIF_CUTOFF -I src/common/ -I src/c/ -c src/omp/health.c -o obj/health_omp_coff.o $(FLAGS_OMP)
	$(CC) obj/health_omp_coff.o obj/bots_common.o obj/bots_main_omp_coff.o obj/bench_omp.o $(FLAGS_OMP) -o bin/health_omp_coff


	$(CC) -DFINAL_CUTOFF -DBENCH_MODE=OPENMP -I src/c/ -I src/omp/ -c src/common/bots_main.c -o obj/bots_main_omp_fin.o -fopenmp
	$(CC) -DFINAL_CUTOFF -I src/common/ -I src/c/ -c src/omp/health.c -o obj/health_omp_fin.o $(FLAGS_OMP)
	$(CC) obj/health_omp_fin.o obj/bots_common.o obj/bots_main_omp_fin.o obj/bench_omp.o $(FLAGS_OMP) -o bin/health_omp_fin


	$(CC) -DMANUAL_CUTOFF -DBENCH_MODE=OPENMP -I src/c/ -I src/omp/ -c src/common/bots_main.c -o obj/bots_main_omp_man.o -fopenmp
	$(CC) -DMANUAL_CUTOFF -I src/common/ -I src/c/ -c src/omp/health.c -o obj/health_omp_man.o $(FLAGS_OMP)
	$(CC) obj/health_omp_man.o obj/bots_common.o obj/bots_main_omp_man.o obj/bench_omp.o $(FLAGS_OMP) -o bin/health_omp_man



run:
	echo "{ \"out\" : ["
	./bin/health_seq -v 0 -f input/test.input
	echo ','
	./bin/health_omp -v 0 -f input/test.input
	echo ','
	./bin/health_omp_coff -v 0 -f input/test.input
	echo ','
	./bin/health_omp_fin -v 0 -f input/test.input
	echo ','
	./bin/health_omp_man -v 0 -f input/test.input
	echo ','
	./bin/health_seq -v 0 -f input/medium.input
	echo ','
	./bin/health_omp -v 0 -f input/medium.input
	echo ','
	./bin/health_omp_coff -v 0 -f input/medium.input
	echo ','
	./bin/health_omp_fin -v 0 -f input/medium.input
	echo ','
	./bin/health_omp_man -v 0 -f input/medium.input
	echo ','
	./bin/health_seq -v 0 -f input/small.input
	echo ','
	./bin/health_omp -v 0 -f input/small.input
	echo ','
	./bin/health_omp_coff -v 0 -f input/small.input
	echo ','
	./bin/health_omp_fin -v 0 -f input/small.input
	echo ','
	./bin/health_omp_man -v 0 -f input/small.input
	echo ','
	./bin/health_seq -v 0 -f input/large.input
	echo ','
	./bin/health_omp -v 0 -f input/large.input
	echo ','
	./bin/health_omp_coff -v 0 -f input/large.input
	echo ','
	./bin/health_omp_fin -v 0 -f input/large.input
	echo ','
	./bin/health_omp_man -v 0 -f input/large.input
	echo "]}"

