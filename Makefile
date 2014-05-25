include ./environment

EXAMPLES := $(wildcard ./examples/*.cu)

all: examples
	make -C ./examples
	
examples: $(EXAMPLES)
	make -C ./examples

$(EXAMPLES):
	
clean:
	make clean -C ./examples
