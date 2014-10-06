include ../make.inc

#will replace the .c files with the .o files once they have been compiled. 
OBJS := $(patsubst ../src/%.cu,../src/%.o,$(wildcard ../src/*.cu)) $(patsubst ../src/%.c,../src/%.o,$(wildcard ../src/*.c))

#patsubst will replace the .c and .cu files with .o 
EXS := $(patsubst ./%.c,./% ,$(wildcard ./*.c)) $(patsubst ./%.cu,./% ,$(wildcard ./*.cu))

all: $(EXS) 

%:: %.c
	$(CUC) $(CUFLAGS) -o $@ $< $(OBJS) -I../include/ $(INCS) $(LIBS) 

%:: %.cu
	$(CUC) $(CUFLAGS) -o $@ $< $(OBJS) -I../include/ $(INCS) $(LIBS) 

clean:
	-rm $(EXS)


