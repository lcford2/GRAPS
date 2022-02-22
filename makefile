
all:
	+$(MAKE) -C src

graps.so:
	+$(MAKE) -C src graps.so

graps:
	+$(MAKE) -C src graps

clean:
	+$(MAKE) -C src clean
	# rm -f lib/graps.so bin/graps src/*.o src/*.out src/*.mod
