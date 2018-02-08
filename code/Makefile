CC = gcc
LIBODE = libode.a
XLIBS = -L. -lode
LOADLIBES =  -L/usr/local/lib -lgsl -lgslcblas -lm
CFLAGS = -g -I/usr/local/include
OBJS = gsl-interface.o   models.o
DEPS = frequencies.h #$(LIBODE) Makefile

PIC = -fPIC

#gaussian: gaussian.c
#	gcc $(CFLAGS) gaussian.c -o gaussian.o $(LOADLIBES)

#%.o:%.c
#	$(GCC) $(CFLAGS) -c  $(PIC) -Wall $(INCLUDE)  $< -o $@
phonons: phonons.c $(DEPS)
	$(CC) $(CFLAGS)  $(PIC) -Wall $(INCLUDE)  $<   -o $@ $(LOADLIBES)

%:%.c #	 $(DEPS)
	$(CC) $(CFLAGS)  $(PIC) -Wall $(INCLUDE)  $<   -o $@ $(LOADLIBES)
#%.so: %.o
	#$(GCC) $(CBITSFLAGS) -shared -Wl,-soname,$@ -o $@ $<

#$(LIBODE): $(OBJS)
#	ar -cr $@  $?

clean:
	rm -rf *.so *.dylib *.o *~ *.a
