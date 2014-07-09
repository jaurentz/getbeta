include make.inc

all:
	make -C $(EXDIR)

clean:
	make clean -C $(EXDIR) &&\
	make clean -C $(SRCDIR)


