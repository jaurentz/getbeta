#this is the home directory of where the project lives. 
HOMEDIR := /media/FC74-9970/getBeta

# C compiler
CC := nvcc
CFLAGS := -O3 -fopenmp

# CUDA compiler
CUC := nvcc
CUFLAGS := -O3 --compiler-options '$(CFLAGS)' -gencode arch=compute_30,code=sm_30

#source file directory 
SRCDIR := $(HOMEDIR)/src

#example file directory 
EXDIR := $(HOMEDIR)/examples

#wildcard picks all the files with .c extension 
SRCS := $(wildcard $(SRCDIR)/*.c) $(wildcard $(SRCDIR)/*.cu) 

#will replace the .c files with the .o files once they have been compiled. 
OBJS := $(patsubst $(SRCDIR)/%.cu,$(SRCDIR)/%.o,$(wildcard $(SRCDIR)/*.cu)) $(patsubst $(SRCDIR)/%.c,$(SRCDIR)/%.o,$(wildcard $(SRCDIR)/*.c))

#wildcard picks all of the files in example directory with .c file extension. 
CEXS := $(wildcard $(EXDIR)/*.c)
CUEXS := $(wildcard $(EXDIR)/*.cu)

#patsubst will replace the .c files with .o 
EXS := $(patsubst $(EXDIR)/%.c,$(EXDIR)/% ,$(CEXS)) $(patsubst $(EXDIR)/%.cu,$(EXDIR)/% ,$(CUEXS))

# cucheb
CUCHEBDIR := /usr/local/cucheb

#calling LAPACKE LAPACK & BLAS 
LIBS := -L $(CUCHEBDIR)/lib -lcucheb -lcudart -lcuda -lcusparse -lcublas -lcurand -lcufft -lm -llapacke -llapack -lblas

#calling the header file. 
INCS := -I $(HOMEDIR)/include/ -I $(CUCHEBDIR)/include -I /usr/include/lapacke/ 
