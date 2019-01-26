CC=gcc
FLAGS_OMP=-fopenmp -lm
FLAGS_SEQ= -lm
FLAGS_PTH= -lm -lpthread
THREADS=4
G=

build:
	$(CC) src/bs-seq.c src/c/bench.c $(FLAGS_SEQ) -o bin/bs-seq
	$(CC) src/bs-pt.c src/c/bench.c $(FLAGS_PTH) -o bin/bs-pt
	$(CC) src/bs-omp.c src/c/bench.c $(FLAGS_OMP) -o bin/bs-omp
	$(CC) src/bs-td.c src/c/bench.c $(FLAGS_OMP) -o bin/bs-td
	$(CC) src/inputgen/inputgen.c -o bin/ig
	./bin/ig  1000 input/small
	./bin/ig  5000 input/medium
	./bin/ig 25000 input/big

run:
	echo "{ \"out\" : ["
	./bin/bs-seq 1 input/small output/seq-small
	echo ","
	./bin/bs-seq 1 input/medium output/seq-medium
	echo ","
	./bin/bs-seq 1 input/big output/seq-big
	echo ","
	./bin/bs-pt $(THREADS) input/small output/seq-small
	echo ","
	./bin/bs-pt $(THREADS) input/medium output/seq-medium
	echo ","
	./bin/bs-pt $(THREADS) input/big output/seq-big
	echo ","
	./bin/bs-omp $(THREADS) input/small output/seq-small
	echo ","
	./bin/bs-omp $(THREADS) input/medium output/seq-medium
	echo ","
	./bin/bs-omp $(THREADS) input/big output/seq-big
	echo ","
	OMP_NUM_THREADS=$(THREADS) ./bin/bs-td 2 input/small output/seq-small 1000
	echo ","
	OMP_NUM_THREADS=$(THREADS) ./bin/bs-td 2 input/medium output/seq-medium 1000
	echo ","
	OMP_NUM_THREADS=$(THREADS) ./bin/bs-td 2 input/big output/seq-big 1000
	echo "]}"

