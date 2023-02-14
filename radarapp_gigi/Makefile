#EPI_TOOLCHAIN=/apps/llvm-toolchain/EPI-0.7/development
EPI_TOOLCHAIN=/apps/llvm-toolchain/EPI/development
EPI_CC=${EPI_TOOLCHAIN}/bin/clang
CC=${EPI_CC}
#EXTRAE_HOME=/apps/BSCTOOLS/extrae/3.8.3
#CFLAGS = -mepi -O2 -mllvm -combiner-store-merging=false -mllvm -vectorizer-use-vp-strided-load-store -fno-slp-vectorize -fno-vectorize -ffast-math -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize -I/apps/vehave/EPI-0.7/development/include/vehave/ -mcpu=avispado #-Rpass-missed=loop-vectorize 

CFLAGS = -Wl,-rpath -Wl,$(EXTRAE_HOME)/lib -L$(EXTRAE_HOME)/lib -I/$(EXTRAE_HOME)/include -lseqtrace -mepi -O2 -mllvm -combiner-store-merging=false -mllvm -vectorizer-use-vp-strided-load-store -fno-slp-vectorize -fno-vectorize -ffast-math -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize -I/apps/riscv/vehave/EPI-0.7/development/include/vehave/ /apps/riscv/fpga-sdv/ila2prv/include/trigger.h -mcpu=avispado  #-Rpass-missed=loop-vectorize 
LDFLAGS = -lm

EXEC = main
SRC = main.c
SRC += covariance_v2.c
SRC += hessenberg.c
#SRC += eigen.c
SRC += eigen_v2.c
SRC += unitary.c
#SRC += qr_iterration.c
SRC += qr_iterrationv2.c
#SRC += MUSIC_algorithm.c
SRC += MUSIC_algorithm_v3.c
.PHONY: all 

OBJS = $(SRC:%.c=%.o)
	
$(EXEC): $(OBJS)

all: dyn st vec

dyn: $(SRC)
	$(CC) $(CFLAGS) $(SRC)  $(LDFLAGS) -o main

st: $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) -o main-st -static

vec: $(SRC)
	$(CC) $(CFLAGS) $(SRC) $(LDFLAGS) -o main-vec 


.PHONY: clean
clean:
	rm -f *.o
	rm -f main
	rm -f main-st
	rm -f main-vec
