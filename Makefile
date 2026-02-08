# spmv-omp-phi — Sparse Matrix-Vector multiplication with OpenMP
# Build with: make [CFLAGS=-O3] [CC=gcc]

CC     ?= gcc
CFLAGS ?= -O2 -Wall
OMP    = -fopenmp
INC    = -I include

SRC = mmio.c util.c dense.c sparse.c blocked.c reorder.c
OBJ = $(SRC:.c=.o)
LIB = libspmv.a

.PHONY: all clean

all: $(LIB)

$(LIB): $(OBJ)
	$(AR) rcs $@ $(OBJ)

%.o: %.c
	$(CC) $(CFLAGS) $(OMP) $(INC) -c -o $@ $<

clean:
	rm -f $(OBJ) $(LIB)
