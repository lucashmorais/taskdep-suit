CC=gcc
FLAGS_OMP=-fopenmp -O3 -lm
FLAGS_SEQ=-O3 -lm
G=

bench:
	$(CC) -c src/common/bots_common.c -o obj/bots_common.o
	$(CC) -DBENCH_MODE=SEQ -I src/c/ -I src/serial/ -c src/common/bots_main.c -o obj/bots_main_seq.o
	$(CC) -DBENCH_MODE=OPENMP -I src/c/ -I src/omp/alignment_single/ -c src/common/bots_main.c -o obj/bots_main_omp_sing.o -fopenmp
	$(CC) -DBENCH_MODE=OPENMP -I src/c/ -I src/omp/alignment_for/ -c src/common/bots_main.c -o obj/bots_main_omp_for.o -fopenmp
	$(CC) -c src/c/bench.c -o obj/bench.o
	$(CC) -c src/c/bench.c -o obj/bench_omp.o -fopenmp

build: bench
	$(CC) -I src/common/ -I src/c/ -c src/serial/sequence.c -o obj/sequence_seq.o
	$(CC) -I src/common/ -I src/c/ -c src/serial/alignment.c -o obj/alignment_seq.o
	$(CC) obj/alignment_seq.o obj/sequence_seq.o obj/bots_common.o obj/bots_main_seq.o obj/bench.o $(FLAGS_SEQ) -o bin/alignment_seq

	$(CC) -I src/common/ -I src/c/ -c src/omp/alignment_single/sequence.c -o obj/sequence_omp_sing.o $(FLAGS_OMP)
	$(CC) -I src/common/ -I src/c/ -c src/omp/alignment_single/alignment.c -o obj/alignment_omp_sing.o $(FLAGS_OMP)
	$(CC) obj/alignment_omp_sing.o obj/sequence_omp_sing.o obj/bots_common.o obj/bots_main_omp_sing.o obj/bench_omp.o $(FLAGS_OMP) -o bin/alignment_omp_sing

	$(CC) -I src/common/ -I src/c/ -c src/omp/alignment_for/sequence.c -o obj/sequence_omp_for.o $(FLAGS_OMP)
	$(CC) -I src/common/ -I src/c/ -c src/omp/alignment_for/alignment.c -o obj/alignment_omp_for.o $(FLAGS_OMP)
	$(CC) obj/alignment_omp_for.o obj/sequence_omp_for.o obj/bots_common.o obj/bots_main_omp_for.o obj/bench_omp.o $(FLAGS_OMP) -o bin/alignment_omp_for

run:
	echo "{ \"out\" : ["
	./bin/alignment_seq -f input/prot.20.aa -v 0
	echo ','
	./bin/alignment_omp_sing -f input/prot.20.aa -v 0
	echo ','
	./bin/alignment_omp_for -f input/prot.20.aa -v 0
	echo ','
	./bin/alignment_seq -f input/prot.100.aa -v 0
	echo ','
	./bin/alignment_omp_sing -f input/prot.100.aa -v 0
	echo ','
	./bin/alignment_omp_for -f input/prot.100.aa -v 0
	echo "]}"

