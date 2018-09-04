FLAGS = -g -DMKL_ILP64 -m64 -I${MKLROOT}/include  ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -lgomp -lpthread -lm -ldl -fopenmp


cholesky:
	gcc cholesky.c -o cholesky $(FLAGS)

