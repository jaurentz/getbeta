# CUCHEB home directory
CUCHEBDIR := /usr/local/cucheb

# directory for CUDA library
CUDADIR := #/usr/local/cuda

# CUDA compiler and flags
CUC := nvcc
CUFLAGS := -O3 --compiler-options '-fPIC -fopenmp' -arch=sm_30

# directory for CULA library
CULADIR := /usr/local/cula

# directory for CULASPARSE library
CULASPARSEDIR := /usr/local/culasparse

# directory for BLAS library
BLASDIR :=

# directory for LAPACK library
LAPACKDIR :=

# linking directories
INCDIRS := -I $(CUCHEBDIR)/include -I $(CUDADIR)/include -I $(CULADIR)/include -I $(CULASPARSEDIR)/include
LIBDIRS := -L $(CUCHEBDIR)/lib -L $(CULADIR)/lib64 -L $(CULASPARSEDIR)/lib64 -L $(BLASDIR) -L $(LAPACKDIR)
LIBS := -lcucheb -lcuda -lcudart -lcublas -lcusparse -lcufft -lcurand -lcula_lapack -lcula_sparse -llapack -lblas -lm 
