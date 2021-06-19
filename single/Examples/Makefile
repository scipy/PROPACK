#
#     (C) Rasmus Munk Larsen, Stanford University, 2004
#
include ../../make.inc
IPATH = -I..
LPATH = $(LIBPATH) -L.. -L.
LIBS = ../libspropack_$(PLAT).a  ../libslapack_util_$(PLAT).a $(BLAS)  ../libslapack_util_$(PLAT).a 
#LIBS = -lspropack_$(PLAT)  -lslapack_util_$(PLAT) $(BLAS)  -lslapack_util_$(PLAT)
#LIBS = -lspropack_$(PLAT) -lmkl_lapack $(BLAS)  


all: example.x example_irl.x compare.x

.F.o:
	$(F77) $(FFLAGS)  $(IPATH) -c -o $*.o $*.F
.f.o:
	$(F77) $(FFLAGS)  $(IPATH) -c -o $*.o $*.f
.c.o:
	$(CC)  $(CFLAGS)  $(IPATH) -c -o $*.o $*.c

### Dependencies ###

matvec.o example.o example_irl.o: matvec.h

### Linking information for executables. ###

lib:
	@( cd ..; \
	$(MAKE) lib; \
	cd Examples )

example.x: example.$(PLAT).x

example_irl.x: example_irl.$(PLAT).x

compare.x: compare.$(PLAT).x

example.$(PLAT).x: lib example.o matvec.o
	$(LINKER) $(LINKFLAGS) -o example.$(PLAT).x example.o matvec.o $(LPATH) $(LIBS)

example_irl.$(PLAT).x: lib example_irl.o matvec.o
	$(LINKER) $(LINKFLAGS) -o example_irl.$(PLAT).x example_irl.o matvec.o $(LPATH) $(LIBS)

compare.$(PLAT).x: lib compare.o 
	$(LINKER) $(LINKFLAGS) -o compare.$(PLAT).x compare.o $(LPATH) $(LIBS)

test:
	@( echo "==Testing REAL version without implicit restart.=="; \
	example.$(PLAT).x <  example.in; \
	echo "==Testing REAL version with implicit restart.=="; \
	example_irl.$(PLAT).x <  example_irl.in)

verify: 
	@( echo "==Comparing REAL results to reference.=="; \
	compare.$(PLAT).x < compare.in; \
	echo "==Comparing REAL IRL results to reference.=="; \
	compare.$(PLAT).x < compare_irl.in )

clean:	
	rm -f  *.o *.il
	rm -rf rii_files

cleanall: 
	rm -f  *.o *.a *.x *~ work* *.il  *.ascii
	rm -rf rii_files


