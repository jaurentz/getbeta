#this is the home directory of where the project lives. 
HOMEDIR := /home/johnny/projects/getBeta

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

#patsubst will replace the .c files with .o 
EXS := $(patsubst $(EXDIR)/%.c,$(EXDIR)/ % ,$(CEXS) )

#calling LAPACKE LAPACK & BLAS 
LIBS := -lm -llapacke -llapack -lblas

#calling the header file. 
INCS := -I ../include/  -I /usr/include/lapacke/
