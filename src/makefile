COMP=gfort
include make_$(COMP).inc

all: graps.so graps
graps.so:
	$(FF) $(FFLAGS) -shared -fpic path_vars.f90 definitions.f90 My_variables.f90 ffsqp.f qld.f init.f90 all_simul.f90 -o../lib/graps.so

graps:
	$(FF) $(FFLAGS) path_vars.f90 definitions.f90 My_variables.f90 ffsqp.f qld.f init.f90 all_simul.f90 multireservoir.f90 -o ../bin/graps

clean:
	rm -f *.o standard_* *.out *.mod ../lib/graps.so ../bin/graps
