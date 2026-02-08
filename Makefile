# spmv-omp-phi — Sparse Matrix-Vector multiplication with OpenMP
# Build with: make [CFLAGS=-O3] [CC=gcc]

CC     ?= gcc
CFLAGS ?= -O2 -Wall
OMP    = -fopenmp
INC    = -I include

SRC_DIR = src
SRC = $(SRC_DIR)/mmio.c $(SRC_DIR)/util.c $(SRC_DIR)/dense.c $(SRC_DIR)/sparse.c $(SRC_DIR)/blocked.c $(SRC_DIR)/reorder.c
OBJ = $(SRC_DIR)/mmio.o $(SRC_DIR)/util.o $(SRC_DIR)/dense.o $(SRC_DIR)/sparse.o $(SRC_DIR)/blocked.o $(SRC_DIR)/reorder.o
LIB = libspmv.a

.PHONY: all clean

all: $(LIB)

$(LIB): $(OBJ)
	$(AR) rcs $@ $(OBJ)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(OMP) $(INC) -c -o $@ $<

clean:
	rm -f $(OBJ) $(LIB)
