CC=gcc
FLAGS_OMP=-fopenmp -O3 -lm
FLAGS_SEQ=-O3 -lm
G=

bench:
	$(CC) $(G) -c src/bench.c -o obj/bench.o
	$(CC) $(G) -c src/bench.c -o obj/bench_omp.o -fopenmp
	$(CC) $(G) -c src/main.c -o obj/main.o
	$(CC) $(G) -c src/strassen.c -o obj/strassen.o
	$(CC) $(G) -c $(FLAGS_OMP) src/main.c -o obj/main-omp.o
	$(CC) $(G) -c $(FLAGS_OMP) src/strassen.c -o obj/strassen-omp.o
	


build: bench
	$(CC) $(G) src/strassen-seq.c obj/main.o obj/bench.o obj/strassen.o $(FLAGS_SEQ) -o bin/strassen-seq
	$(CC) $(G) src/strassen-task.c obj/main-omp.o obj/bench_omp.o obj/strassen-omp.o $(FLAGS_OMP) -o bin/strassen-task
	$(CC) $(G) src/strassen-task-dep.c obj/main-omp.o obj/bench_omp.o obj/strassen-omp.o $(FLAGS_OMP) -o bin/strassen-task-dep

run:
	echo "{ \"out\" : ["
	./bin/strassen-seq -i 3 -n 512 -s 64 -d 5
	echo ','
	./bin/strassen-task -i 16 -n 512 -s 64 -d 5
	echo ','
	./bin/strassen-task-dep -i 16 -n 512 -s 64 -d 5
	echo "]}"
