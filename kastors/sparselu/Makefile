include ../../conf
JACOBI=-DMSIZE -DSMSIZE

build:
	$(CC) src/sparselu.c src/main.c src/sparselu-seq.c $(BENCH_PATH)bench.c $(FLAGS_BENCH) $(JACOBI) -lm -o bin/sparselu-seq
	$(MCC) src/sparselu-ompss.c src/sparselu.c src/main.c src/sparselu-seq.c $(BENCH_PATH)bench.c $(FLAGS_BENCH) $(JACOBI) -lm --ompss -D_OMPSS -o bin/sparselu-ompss
	$(MCC) src/sparselu-task-dep.c src/sparselu.c src/main.c src/sparselu-seq.c $(BENCH_PATH)bench.c $(FLAGS_BENCH) $(JACOBI) -lm --openmp -o bin/sparselu-omp
run:
	echo "{ \"out\" : ["
	./bin/sparselu-seq -n 128 -m 16
	echo ','
	./bin/sparselu-seq -n 128 -m 16
	echo ','
	./bin/sparselu-omp -n 128 -m 16 
	echo ','
	./bin/sparselu-omp -n 128 -m 16 
	echo ','
	./bin/sparselu-omp -n 128 -m 16 
	echo ','
	./bin/sparselu-omp -n 128 -m 16 
	echo ','
	./bin/sparselu-omp -n 128 -m 16 
	echo ','
	./bin/sparselu-omp -n 128 -m 16 
	echo ','
	./bin/sparselu-omp -n 128 -m 16 
	echo ','
	./bin/sparselu-omp -n 128 -m 16 
	echo ','
	./bin/sparselu-ompss -n 128 -m 16 
	echo ','
	./bin/sparselu-ompss -n 128 -m 16 
	echo ','
	./bin/sparselu-ompss -n 128 -m 16 
	echo ','
	./bin/sparselu-ompss -n 128 -m 16 
	echo ','
	./bin/sparselu-ompss -n 128 -m 16 
	echo ','
	./bin/sparselu-ompss -n 128 -m 16 
	echo ','
	./bin/sparselu-ompss -n 128 -m 16 
	echo ','
	./bin/sparselu-omps -n 128 -m 16
	echo "]}"
