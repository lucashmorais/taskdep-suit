include ../../conf
JACOBI=-DMSIZE -DBSIZE -DTITER

CC=riscv64-buildroot-linux-gnu-g++
CXX=riscv64-buildroot-linux-gnu-g++
MCC=/home/lucas/Downloads/Installations/riscv/mercurium/bin/riscv64-buildroot-linux-gnu-mcc
#FLAGS_BENCH=-DZERO_PACKETS_V2
FLAGS_BENCH=

n6vers:
	$(MCC) $(ARCH_CCARGS) src/jacobi-ompss2.c src/jacobi-main.c src/main.c src/jacobi-seq.c src/poisson.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm --ompss-2 -D_OMPSS2 -o bin/jacobi-ompss2

build:
	$(PREFIX)$(CC) $(ARCH_CCARGS) src/jacobi-main.c src/main.c src/jacobi-seq.c src/poisson.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm -o bin/jacobi-seq
	$(PREFIX)$(CXX) $(ARCH_CCARGS) src/jacobi-phentos-ord.cpp src/jacobi-seq.c src/jacobi-main.c src/main.c src/poisson.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm -lpthread -DPHENTOS -o bin/jacobi-phentos-ord -w
	$(PREFIX)$(CXX) $(ARCH_CCARGS) src/jacobi-phentos-fast.cpp src/jacobi-seq.c src/jacobi-main.c src/main.c src/poisson.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm -lpthread -DPHENTOS -o bin/jacobi-phentos-fast -w
	#$(PREFIX)$(MCC) $(ARCH_CCARGS) src/jacobi-ompss.c src/jacobi-main.c src/main.c src/jacobi-seq.c src/poisson.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm --ompss -D_OMPSS -o bin/jacobi-ompss
#	$(MCC) src/jacobi-task-dep.c src/jacobi-main.c src/main.c src/jacobi-seq.c src/poisson.c $(BENCH_PATH)bench.c $(BENCH_PATH)sha256.c $(FLAGS_BENCH) $(JACOBI) -lm --openmp -o bin/jacobi-omp
run:
	echo "{ \"out\" : ["
	BENCH_SILENT= ./bin/jacobi-seq -n 2048 -b 512
	echo ','
	BENCH_SILENT= ./bin/jacobi-seq -n 2048 -b 512
	echo ','
	# BENCH_SILENT= ./bin/jacobi-omp -n 2048 -b 512 
	# echo ','
	# BENCH_SILENT= ./bin/jacobi-omp -n 2048 -b 512 
	# echo ','
	# BENCH_SILENT= ./bin/jacobi-omp -n 2048 -b 512 
	# echo ','
	# BENCH_SILENT= ./bin/jacobi-omp -n 2048 -b 512 
	# echo ','
	# BENCH_SILENT= ./bin/jacobi-omp -n 2048 -b 512 
	# echo ','
	# BENCH_SILENT= ./bin/jacobi-omp -n 2048 -b 512 
	# echo ','
	# BENCH_SILENT= ./bin/jacobi-omp -n 2048 -b 512 
	# echo ','
	# BENCH_SILENT= ./bin/jacobi-omp -n 2048 -b 512 
	# echo ','
	BENCH_SILENT= ./bin/jacobi-ompss -n 2048 -b 512 
	echo ','
	BENCH_SILENT= ./bin/jacobi-ompss -n 2048 -b 512 
	echo ','
	BENCH_SILENT= ./bin/jacobi-ompss -n 2048 -b 512 
	echo ','
	BENCH_SILENT= ./bin/jacobi-ompss -n 2048 -b 512 
	echo ','
	BENCH_SILENT= ./bin/jacobi-ompss -n 2048 -b 512 
	echo ','
	BENCH_SILENT= ./bin/jacobi-ompss -n 2048 -b 512 
	echo ','
	BENCH_SILENT= ./bin/jacobi-ompss -n 2048 -b 512 
	echo ','
	BENCH_SILENT= ./bin/jacobi-ompss -n 2048 -b 512
	echo "]}"
copy-to-overlay:
	cp ./bin/* ~/Downloads/Compressed/buildroot-2021.08.1/overlay_dir/root/kastors/jacobi/
