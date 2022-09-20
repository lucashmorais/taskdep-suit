include ../../conf

MCC=/home/lucas/Downloads/Installations/riscv/mercurium/bin/riscv64-buildroot-linux-gnu-mcc
JACOBI=-DMSIZE -DSMSIZE
#FLAGS_BENCH=-DZERO_PACKETS_V2
FLAGS_BENCH=

n6vers:
	$(MCC) $(ARCH_CCARGS) src/sparselu-ompss2.c src/sparselu.c src/main.c src/sparselu-seq.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm --ompss-2 -D_OMPSS2 -o bin/sparselu-ompss2

build:
	$(PREFIX)$(CC) $(ARCH_CCARGS) src/sparselu.c src/main.c src/sparselu-seq.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm -o bin/sparselu-seq -Wno-write-strings 
	$(PREFIX)$(CXX) $(ARCH_CCARGS) src/sparselu-phentos-fast.cpp src/sparselu.c src/main.c src/sparselu-seq.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm -o bin/sparselu-phentos-fast -Wno-write-strings -lpthread -D_PHENTOS 
	$(PREFIX)$(CXX) $(ARCH_CCARGS) src/sparselu-phentos-ord.cpp src/sparselu.c src/main.c src/sparselu-seq.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm -o bin/sparselu-phentos-ord -Wno-write-strings -lpthread -D_PHENTOS 
	#$(PREFIX)$(MCC) $(ARCH_CCARGS) src/sparselu-ompss.c src/sparselu.c src/main.c src/sparselu-seq.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm --ompss -D_OMPSS -o bin/sparselu-ompss -Wno-write-strings
	#$(MCC) src/sparselu-task-dep.c src/sparselu.c src/main.c src/sparselu-seq.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm --openmp -o bin/sparselu-omp -Wno-write-strings
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
copy-to-overlay:
	cp ./bin/* ~/Downloads/Compressed/buildroot-2021.08.1/overlay_dir/root/kastors/sparselu/
