include ../make.inc

all: $(EXS) 

%:: %.c $(OBJS)
	$(CC) $(CFlAGS) -o $@ $< $(OBJS) $(LIBS) $(INCS) 

%:: %.cu $(OBJS)
	$(CUC) $(CUFlAGS) -o $@ $< $(OBJS) $(LIBS) $(INCS) 

$(OBJS): $(SRCS)
	make -C $(SRCDIR)

$(SRCS):

clean:
	-rm $(EXS)


