# CUDA compiler
CUC := /usr/local/cuda-6.0/bin/nvcc
CUFLAGS := --compiler-options '-O3 -fopenmp' -gencode arch=compute_35,code=sm_35 -gencode arch=compute_30,code=sm_30

# cucheb
CUCHEBDIR := /home/aurentz/Mylibs/cucheb

#calling LAPACKE LAPACK & BLAS 
LIBS := -L$(CUCHEBDIR)/lib -lcucheb -lcudart -lcuda -lcusparse -lcublas -lcurand -lcufft \
-L/home/aurentz/Mylibs/lapack-3.5.0 -llapacke -llapack -lblas -lm

#calling the header file. 
INCS := -I$(CUCHEBDIR)/include -I/home/aurentz/Mylibs/lapack-3.5.0/lapacke/include 
