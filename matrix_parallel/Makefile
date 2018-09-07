ALLSOURCE=$(wildcard *.f90)
LIBDEVICE_OBJ=$(patsubst %.f90,%.o,$(ALLSOURCE))
DEPENDALL=../staticparameters.f90
MODFILES=mtmod.mod
include ../Makefile.inc

libmparallel.a: ${MODFILES} ${LIBDEVICE_OBJ}
	${AR} ${ARFLAGS} $@ $^
	${RANLIB} $@

mtmod.mod: mt19937.f90
	${FC} -c ${FCFLAGS} $<

clean:
	echo cleaning up in hostsource.
	$(RM) -f *.o *.a *.mod
