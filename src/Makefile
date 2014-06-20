include ../make.inc

all: $(OBJS) 

%.o: %.c
	$(CC) $(CFlAGS) -c $< $(LIBS) $(INCS)

%.o: %.cu
	$(CUC) $(CUFlAGS) -c $< $(LIBS) $(INCS)

clean:
	-rm $(OBJS)


