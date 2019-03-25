include ../../conf
JACOBI=-DMSIZE -DSMSIZE

build:
	$(PREFIX)$(CC) $(ARCH_CCARGS) src/sparselu.c src/main.c src/sparselu-seq.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm -o bin/sparselu-seq
	$(PREFIX)$(MCC) $(ARCH_CCARGS) src/sparselu-ompss.c src/sparselu.c src/main.c src/sparselu-seq.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm --ompss -D_OMPSS -o bin/sparselu-ompss
	#$(MCC) src/sparselu-task-dep.c src/sparselu.c src/main.c src/sparselu-seq.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm --openmp -o bin/sparselu-omp
run:
	echo "{ \"out\" : ["
	BENCH_SILENT= ./bin/sparselu-seq -n 64 -m 8
	echo ','
	BENCH_SILENT= ./bin/sparselu-seq -n 64 -m 8
	echo ','
	BENCH_SILENT= ./bin/sparselu-omp -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-omp -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-omp -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-omp -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-omp -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-omp -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-omp -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-omp -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-ompss -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-ompss -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-ompss -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-ompss -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-ompss -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-ompss -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-ompss -n 64 -m 8 
	echo ','
	BENCH_SILENT= ./bin/sparselu-ompss -n 64 -m 8
	echo "]}"
