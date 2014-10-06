include ./make.inc

#will replace the .c files with the .o files once they have been compiled. 
OBJS := $(patsubst ./src/%.cu,./src/%.o,$(wildcard ./src/*.cu)) $(patsubst ./src/%.c,./src/%.o,$(wildcard ./src/*.c))

SRCS := $(wildcard ./src/*.cu) $(wildcard ./src/*.c)

#patsubst will replace the .c and .cu files with .o 
EXS := $(patsubst ./examples/%.c,./examples/% ,$(wildcard ./examples/*.c)) $(patsubst ./examples/%.cu,./examples/% ,$(wildcard ./examples/*.cu))

EXSRCS := $(wildcard ./examples/*.c) $(wildcard ./examples/*.cu)

all: $(EXS)

$(EXS): $(OBJS) $(EXSRCS)
	make -C ./examples
	
$(EXSRCS):
	
$(OBJS): $(SRCS)
	make -C ./src
	
$(SRCS):

clean:
	make clean -C ./examples &&\
	make clean -C ./src


