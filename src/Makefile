include ../make.inc

all: $(OBJS)

%.o: %.c
	$(CUC) $(CUFLAGS) -c $< $(LIBS) $(INCS)

%.o: %.cu
	$(CUC) $(CUFLAGS) -c $< $(LIBS) $(INCS)

clean:
	-rm $(OBJS)


