build:
	cd src; make; mv swaptions-serial ../bin; make clean;
	cd src; version=ompss make; mv swaptions-ompss ../bin; version=ompss make clean;

run:
	echo "{ \"out\" : ["
	BENCH_SILENT= BENCH_OUT=output/seq-32 ./bin/swaptions-serial -ns 32 -sm 512 -nt 1 -sd 1
	echo ","
	BENCH_SILENT= BENCH_OUT=output/seq-128 ./bin/swaptions-serial -ns 128 -sm 2056 -nt 1 -sd 1
	echo ","
	BENCH_SILENT= BENCH_OUT=output/ompss-32 ./bin/swaptions-ompss -ns 32 -sm 512 -sd 1
	echo ","
	BENCH_SILENT= BENCH_OUT=output/ompss-128 ./bin/swaptions-ompss -ns 128 -sm 2056 -sd 1
	echo ","
	BENCH_SILENT= BENCH_OUT=output/ompss-32 ./bin/swaptions-ompss -ns 32 -sm 512 -sd 1
	echo ","
	BENCH_SILENT= BENCH_OUT=output/ompss-128 ./bin/swaptions-ompss -ns 128 -sm 2056 -sd 1
	echo ","
	BENCH_SILENT= BENCH_OUT=output/ompss-32 ./bin/swaptions-ompss -ns 32 -sm 512 -sd 1
	echo ","
	BENCH_SILENT= BENCH_OUT=output/ompss-128 ./bin/swaptions-ompss -ns 128 -sm 2056 -sd 1
	echo ","
	BENCH_SILENT= BENCH_OUT=output/ompss-32 ./bin/swaptions-ompss -ns 32 -sm 512 -sd 1
	echo ","
	BENCH_SILENT= BENCH_OUT=output/ompss-128 ./bin/swaptions-ompss -ns 128 -sm 2056 -sd 1
	echo ","
	BENCH_SILENT= BENCH_OUT=output/ompss-32 ./bin/swaptions-ompss -ns 32 -sm 512 -sd 1
	echo ","
	BENCH_SILENT= BENCH_OUT=output/ompss-128 ./bin/swaptions-ompss -ns 128 -sm 2056 -sd 1
	echo "]}"
