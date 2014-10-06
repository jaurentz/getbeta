include ../make.inc

#will replace the .c files with the .o files once they have been compiled. 
OBJS := $(patsubst ./%.cu,./%.o,$(wildcard ./*.cu)) $(patsubst ./%.c,./%.o,$(wildcard ./*.c))

all: $(OBJS)

%.o: %.c
	$(CUC) $(CUFLAGS) -c $< -I../include/ $(INCS) 

%.o: %.cu
	$(CUC) $(CUFLAGS) -c $< -I../include/ $(INCS)

clean:
	-rm *.o


