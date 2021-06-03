PROGRAM=stream-barr

#TARGETS=$(PROGRAM)-p $(PROGRAM)-i $(PROGRAM)-d
TARGETS=$(PROGRAM)-phentos-ord $(PROGRAM)-phentos-fast

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
EXTRA = -std=c99 -O3 -Wall -Wno-unused
EXTRA_PHENTOS = -O3 -Wall -Wno-unused -DPHENTOS


$(PROGRAM)-phentos-ord: $(PROGRAM)-phentos-ord.cpp
	$(CXX) $(GFLAGS) $(CFLAGS) $(CFLAGS_P) $(EXTRA_PHENTOS) $(INCS) -o $@ $< $(LIBS)

$(PROGRAM)-phentos-fast: $(PROGRAM)-phentos-fast.cpp
	$(CXX) $(GFLAGS) $(CFLAGS) $(CFLAGS_P) $(EXTRA_PHENTOS) $(INCS) -o $@ $< $(LIBS)

$(PROGRAM)-p: $(PROGRAM).c
	$(CC) $(GFLAGS) $(CFLAGS) $(CFLAGS_P) $(EXTRA) $(INCS) -o $@ $< $(LIBS)

$(PROGRAM)-i: $(PROGRAM).c
	$(CC) $(GFLAGS) $(CFLAGS) $(CFLAGS_I) $(EXTRA) $(INCS) -o $@ $< $(LIBS)

$(PROGRAM)-d: $(PROGRAM).c
	$(CC) $(GFLAGS) $(CFLAGS) $(CFLAGS_D) $(EXTRA) $(INCS) -o $@ $< $(LIBS)

clean:
	rm -f $(CC)_* *.o *~ $(TARGETS) 

