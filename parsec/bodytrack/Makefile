build:
	cd src; make; mv dedup-serial ../bin; make clean;
	cd src; version=ompss make; mv dedup-ompss ../bin; version=ompss make clean;

run:
	echo "{ \"out\" : ["
	BENCH_SILENT= ./bin/dedup-serial -c -w none -i input/00.txt output/seq-00.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-serial -c -w none -i input/01.txt output/seq-01.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-serial -c -w none -i input/02.txt output/seq-02.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-ompss -c -w none -i input/00.txt output/omps-00.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-ompss -c -w none -i input/01.txt output/ompss-01.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-ompss -c -w none -i input/02.txt output/ompss-02.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-ompss -c -w none -i input/00.txt output/omps-00.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-ompss -c -w none -i input/01.txt output/ompss-01.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-ompss -c -w none -i input/02.txt output/ompss-02.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-ompss -c -w none -i input/00.txt output/omps-00.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-ompss -c -w none -i input/01.txt output/ompss-01.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-ompss -c -w none -i input/02.txt output/ompss-02.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-ompss -c -w none -i input/00.txt output/omps-00.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-ompss -c -w none -i input/01.txt output/ompss-01.ddp
	echo ","
	BENCH_SILENT= ./bin/dedup-ompss -c -w none -i input/02.txt output/ompss-02.ddp
	echo "]}"
