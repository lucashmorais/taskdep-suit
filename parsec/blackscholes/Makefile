include ../../conf
JACOBI=-DCHK_ERR

build:
	$(CC) src/bs-seq.c $(BENCH_PATH)bench.c $(FLAGS_BENCH) $(JACOBI) -lm -o bin/bs-seq
	$(MCC) src/bs-ompss.c $(BENCH_PATH)bench.c $(FLAGS_BENCH) $(JACOBI) -lm --ompss -D_OMPSS -o bin/bs-ompss
	$(MCC) src/bs-td.c $(BENCH_PATH)bench.c $(FLAGS_BENCH) $(JACOBI) -lm --openmp -o bin/bs-omp
run:
	echo "{ \"out\" : ["
	./bin/bs-seq 1 input/small output/seq-small
	echo ','
	./bin/bs-seq 1 input/small output/seq-small
	echo ','
	./bin/bs-seq 1 input/medium output/seq-medium
	echo ','
	./bin/bs-seq 1 input/medium output/seq-medium
	echo ','
	./bin/bs-seq 1 input/big output/seq-big
	echo ','
	./bin/bs-seq 1 input/big output/seq-big
	echo ','
	./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	./bin/bs-ompss 1 input/big output/ompss-big 500
	echo "]}"
