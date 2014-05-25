include ../environment

all: testcucheb

testcucheb: testcucheb.cu
	$(CUC) $(CUFLAGS) -o testcucheb testcucheb.cu $(INCDIRS) $(LIBDIRS) $(LIBS)
	
clean:
	-rm testcucheb
