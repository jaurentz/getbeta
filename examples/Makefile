include ../make.inc

#EXS := cedar

#OBJS := ../src/trapz.o ../src/innerproduct.o ../src/laplacian.o ../src/polarization.o  ../src/potential.o ../src/x_operator.o  

#LIBS := -lm -llapacke -llapack -lblas

#INCS := -I ../include/ 


all: $(EXS) 

%:: %.c $(OBJS)
	gcc -o $@ $<  $(LIBS) $(INCS) $(OBJS)

$(OBJS) : 
	make -C ../src/

clean:
	-rm $(EXS)


