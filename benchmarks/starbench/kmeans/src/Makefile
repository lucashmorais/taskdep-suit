CC=gcc
SRC_DIR=$(BENCH_DIR)/src
SRC_OBJS=$(SRC_DIR)/wtime.c $(SRC_DIR)/file_io.c $(SRC_DIR)/kmeans_ompss.c $(SRC_DIR)/main.c -lm
BENCH_FLAGS=-g -O2 -std=gnu99
#INPUT_FLAGS=-i ../../input/edge - 16 -b -o ../../output/par
INPUT_FLAGS=-i $(INPUT_FOLDER)/starbench/kmeans/edge - 16 -b -o ../../output/par
INPUT_FLAGS_SEQUENTIAL=-i $(INPUT_FOLDER)/starbench/kmeans/edge - 16 -b -o ../../output/seq -s
