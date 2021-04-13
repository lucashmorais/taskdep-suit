include ../../conf
JACOBI=-DCHK_ERR

build:
	$(PREFIX)$(CC) $(ARCH_CCARGS) src/bs-seq.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm -o bin/bs-seq
	$(PREFIX)$(MCC) $(ARCH_CCARGS) src/bs-ompss.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm --ompss -D_OMPSS -o bin/bs-ompss
	#$(MCC) src/bs-td.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm --openmp -o bin/bs-omp
run:
	echo "{ \"out\" : ["
	BENCH_SILENT= ./bin/bs-seq 1 input/small output/seq-small
	echo ','
	BENCH_SILENT= ./bin/bs-seq 1 input/small output/seq-small
	echo ','
	BENCH_SILENT= ./bin/bs-seq 1 input/medium output/seq-medium
	echo ','
	BENCH_SILENT= ./bin/bs-seq 1 input/medium output/seq-medium
	echo ','
	BENCH_SILENT= ./bin/bs-seq 1 input/big output/seq-big
	echo ','
	BENCH_SILENT= ./bin/bs-seq 1 input/big output/seq-big
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/small output/ompss-small 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/medium output/ompss-medium 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/big output/ompss-big 500
	echo ','
	BENCH_SILENT= ./bin/bs-ompss 1 input/big output/ompss-big 500
	echo "]}"
