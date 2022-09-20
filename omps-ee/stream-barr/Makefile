PROGRAM=stream-barr

#TARGETS=$(PROGRAM)-p $(PROGRAM)-i $(PROGRAM)-d
TARGETS=$(PROGRAM)-phentos-ord-${TARGET_N} $(PROGRAM)-phentos-fast-${TARGET_N}

JOB_SCHED_VERSION=-smp

BASE_DIR=..
include $(BASE_DIR)/common-files/Makefile

# CC = mcc
CC = riscv64-linux-gnu-gcc
CXX = riscv64-linux-gnu-g++

# CFLAGS = --ompss
CFLAGS =
CFLAGS_P =
CFLAGS_I = --instrument
CFLAGS_D = --debug

LIBS  = -fopenmp
INCS  = 

EXTRA = -std=c99 -O3 -Wall -Wno-unused -mcmodel=large 
EXTRA = -std=c99 -O3 -Wall -Wno-unused -DN=${TARGET_N}
#EXTRA_PHENTOS = -O3 -Wall -Wno-unused -DPHENTOS -DZERO_PACKETS_V2 -DN=${TARGET_N}
EXTRA_PHENTOS = -O3 -Wall -Wno-unused -DPHENTOS -DN=${TARGET_N}

$(PROGRAM)-phentos-ord-${TARGET_N}: $(PROGRAM)-phentos-ord.cpp
	echo $(PROGRAM)-phentos-ord-${TARGET_N}
	$(CXX) $(GFLAGS) $(CFLAGS) $(CFLAGS_P) $(EXTRA_PHENTOS) $(INCS) -o $@ $< $(LIBS)

$(PROGRAM)-phentos-fast-${TARGET_N}: $(PROGRAM)-phentos-fast.cpp
	echo $(PROGRAM)-phentos-fast-${TARGET_N}
	$(CXX) $(GFLAGS) $(CFLAGS) $(CFLAGS_P) $(EXTRA_PHENTOS) $(INCS) -o $@ $< $(LIBS)

$(PROGRAM)-p-${TARGET_N}: $(PROGRAM).c
	$(CC) $(GFLAGS) $(CFLAGS) $(CFLAGS_P) $(EXTRA) $(INCS) -o $@ $< $(LIBS)

$(PROGRAM)-i: $(PROGRAM).c
	$(CC) $(GFLAGS) $(CFLAGS) $(CFLAGS_I) $(EXTRA) $(INCS) -o $@ $< $(LIBS)

$(PROGRAM)-d: $(PROGRAM).c
	$(CC) $(GFLAGS) $(CFLAGS) $(CFLAGS_D) $(EXTRA) $(INCS) -o $@ $< $(LIBS)

clean:
	rm -f $(CC)_* *.o *~ $(TARGETS) 

