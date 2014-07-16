#this is the home directory of where the project lives. 
HOMEDIR := /home/johnny/project/getBeta

# C compiler
CC := gcc
CFLAGS := #-O3

# CUDA compiler
CUC := nvcc
CUFLAGS := --compiler-options '$(CFLAGS)' -gencode arch=compute_30,code=sm_30



#source file directory 
SRCDIR := $(HOMEDIR)/src

#example file directory 
EXDIR := $(HOMEDIR)/examples

#wildcard picks all the files with .c extension 
SRC := $(wildcard $(SRCDIR)/*.c)

#will replace the .c files with the .o files once they have been compiled. 
OBJS := $(patsubst $(SRCDIR)/%.c,$(SRCDIR)/%.o,$(SRC))

#wildcard picks all of the files in example directory with .c file extension. 
CEXS := $(wildcard $(EXDIR)/*.c)
CUEXS := $(wildcard $(EXDIR)/*.cu)

#patsubst will replace the .c files with .o 
EXS := $(patsubst $(EXDIR)/%.c,$(EXDIR)/% ,$(CEXS)) $(patsubst $(EXDIR)/%.cu,$(EXDIR)/% ,$(CUEXS))

# cucheb
CUCHEBDIR := /usr/local/cucheb

#calling LAPACKE LAPACK & BLAS 
LIBS := -L $(CUCHEBDIR)/lib -lcucheb -lcuda -lcusparse -lcublas -lcurand -lcufft -lm -llapacke -llapack -lblas

#calling the header file. 
INCS := -I $(HOMEDIR)/include/ -I $(CUCHEBDIR)/include -I /usr/include/lapacke/ 
