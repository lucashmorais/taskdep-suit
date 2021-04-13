include ../../conf
JACOBI=-DMSIZE -DCUTOFF_SIZE -DCUTOFF_DEPTH

build:
	$(CC) src/strassen.c src/main.c src/strassen-seq.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm -o bin/strassen-seq
	$(MCC) src/strassen-ompss.c src/strassen.c src/main.c src/strassen-seq.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm --ompss -D_OMPSS -o bin/strassen-ompss
	$(MCC) src/strassen-task-dep.c src/strassen.c src/main.c src/strassen-seq.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm --openmp -o bin/strassen-omp
run:
	echo "{ \"out\" : ["
	BENCH_SILENT= ./bin/strassen-seq -n 512 -s 64 -d 5
	echo ','
	BENCH_SILENT= ./bin/strassen-seq -n 512 -s 64 -d 5
	echo ','
	BENCH_SILENT= ./bin/strassen-omp -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-omp -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-omp -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-omp -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-omp -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-omp -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-omp -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-omp -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-ompss -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-ompss -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-ompss -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-ompss -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-ompss -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-ompss -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-ompss -n 512 -s 64 -d 5 
	echo ','
	BENCH_SILENT= ./bin/strassen-ompss -n 512 -s 64 -d 5
	echo "]}"
