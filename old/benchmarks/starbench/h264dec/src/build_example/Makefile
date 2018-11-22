BENCH_DIR=../..
BENCH_NAME=h264dec

PLATFORM_MK=$(BENCH_DIR)/../../../common/c.mk
include $(PLATFORM_MK)

#BENCH_MK=../../Makefile
BENCH_MK=$(BENCH_DIR)/../Makefile
include $(BENCH_MK)

all:
	cd ..;\
	autoreconf -i -f; \
	cd build_example; \

$(BENCH_DIR)/build/$(BENCH_NAME):
	make compile

$(BENCH_DIR)/build:
	mkdir $(BENCH_DIR)/build

$(BENCH_DIR)/log:
	mkdir $(BENCH_DIR)/log

$(BENCH_DIR)/build/mtsp:
	mkdir $(BENCH_DIR)/build/mtsp

$(BENCH_DIR)/build/iomp:
	mkdir $(BENCH_DIR)/build/iomp

$(BENCH_DIR)/build/gomp:
	mkdir $(BENCH_DIR)/build/gomp

$(BENCH_DIR)/build/seq:
	mkdir $(BENCH_DIR)/build/seq

$(BENCH_DIR)/build/tioga:
	mkdir $(BENCH_DIR)/build/tioga

compile: $(BENCH_DIR)/build $(BENCH_DIR)/log all
ifeq ($(OMP_LIB),gomp)
	make compileGOMP
endif
ifeq ($(OMP_LIB),iomp)
	make compileIOMP
endif
ifeq ($(OMP_LIB),mtsp)
	make compileIOMP
endif
ifeq ($(OMP_LIB),tioga)
	make compileIOMP
endif
ifeq ($(OMP_LIB),seq)
	make compileSEQ
endif

compileSEQ: $(BENCH_DIR)/build $(BENCH_DIR)/log $(BENCH_DIR)/build/seq
ifeq ($(CPP_FLAG),YES)
	$(eval CC = $(CXX))
endif
	cp ../libavcodec/MakefileSource ../libavcodec/Makefile.am; \
	cd $(BENCH_DIR)/build/seq; \
	../../src/configure CC=$(CC) --enable-openmp; \
	make
	

compileGOMP: $(BENCH_DIR)/build $(BENCH_DIR)/log $(BENCH_DIR)/build/gomp
ifeq ($(CPP_FLAG),YES)
	$(eval CC_GOMP = $(CXX_GOMP))
endif
	cp ../libavcodec/MakefileSource ../libavcodec/Makefile.am; \
	cd $(BENCH_DIR)/build/gomp;\
	../../src/configure CC=$(CC_GOMP) --enable-openmp; \
	make

compileIOMP: $(BENCH_DIR)/build $(BENCH_DIR)/log $(BENCH_DIR)/build/iomp
ifeq ($(CPP_FLAG),YES)
	$(eval CC_IOMP = $(CXX_IOMP))
endif
	cd ../libavcodec;\
	head -n 48 ./MakefileSource > tmp; \
	echo "AM_CFLAGS = -I$(INCLUDE_PATH)" >> tmp; \
	mv tmp Makefile.am
	cd $(BENCH_DIR)/build/iomp;\
	echo "LIBRARY_PATH=$(LIBRARY_PATH) ../src/configure CC=$(CC_IOMP) --enable-openmp --enable-iomp";\
	LIBRARY_PATH=$(LIBRARY_PATH) ../../src/configure CC=$(CC_IOMP) --enable-openmp --enable-iomp --host=$(HOST_ARCHITECTURE); \
	make
	

