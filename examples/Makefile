include ../make.inc

all: $(EXS) 

%:: %.c $(OBJS)
	$(CUC) $(CUFLAGS) -o $@ $< $(OBJS) $(LIBS) $(INCS) 

%:: %.cu $(OBJS)
	$(CUC) $(CUFLAGS) -o $@ $< $(OBJS) $(LIBS) $(INCS) 

$(OBJS): $(SRCS)
	make -C $(SRCDIR)

$(SRCS):

clean:
	-rm $(EXS)


