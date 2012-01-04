
foaw: foaw.o

run: foaw
	./foaw | ./plot.py

velocity.so: velocity.c
	gcc -O3 -shared -fPIC -I/usr/include/python2.6 -o $@ $< -lpython2.6

velocity.c: velocity.pyx
	cython $<

cvelocity.so: cvelocity.c
	gcc -Wall -Werror -O3 -shared -fPIC cvelocity.so cvelocity.c
