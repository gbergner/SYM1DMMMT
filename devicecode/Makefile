ALLSOURCE=$(wildcard *.f90)
LIBDEVICE_OBJ=$(patsubst %.f90,%.o,$(ALLSOURCE))
# be careful with the order of the module files
BASEMODFILES=compiletimeconstants.mod outputstreamnumbers.mod gammamatrix.mod cublasinterface.mod
MODFILES=timer.mod adjust_margins.mod dirac_operator.mod cgm_solver.mod fourier_transform.mod rand_fermion_fields.mod  utils_measurements.mod lattice_action.mod hmc_hamiltonean.mod hmc_force.mod hmc_molecular_dynamics.mod rhmc_updater.mod
DEPENDALL=${BASEMODFILES}
include ../Makefile.inc

libdevice.a: ${BASEMODFILES} ${MODFILES}
	${AR} ${ARFLAGS} $@ ${LIBDEVICE_OBJ}
	${RANLIB} $@

# The mod-files are listed explicitely.
compiletimeconstants.mod: compiletimeconstants.f90 ../staticparameters.f90
	${FC} -c ${FCFLAGS} $<

gammamatrix.mod: gammamatrix.f90
	${FC} -c ${FCFLAGS} $<

outputstreamnumbers.mod: outputstreamnumbers.f90 ../unit_number.inc
	${FC} -c ${FCFLAGS} $<

cublasinterface.mod: cublasinterface.f90
	${FC} -c ${FCFLAGS} $<

timer.mod: timer.f90
	${FC} -c ${FCFLAGS} $<
	
adjust_margins.mod: adjust_margins.f90 ${BASEMODFILES}
	${FC} -c ${FCFLAGS} $<

hmc_molecular_dynamics.mod: hmc_molecular_dynamics.f90 ${BASEMODFILES}
	${FC} -c ${FCFLAGS} $<

lattice_action.mod: lattice_action.f90 ${BASEMODFILES}
	${FC} -c ${FCFLAGS} $<

hmc_hamiltonean.mod:  hmc_hamiltonean.f90 ${BASEMODFILES}
	${FC} -c ${FCFLAGS} $<

hmc_force.mod: hmc_force.f90 ${BASEMODFILES}
	${FC} -c ${FCFLAGS} $<

dirac_operator.mod: dirac_operator.f90 ${BASEMODFILES}
	${FC} -c ${FCFLAGS} $<

cgm_solver.mod: cgm_solver.f90 ${BASEMODFILES}
	${FC} -c ${FCFLAGS} $<

fourier_transform.mod: fourier_transform.f90 ${BASEMODFILES}
	${FC} -c ${FCFLAGS} $<

rand_fermion_fields.mod: rand_fermion_fields.f90 ${BASEMODFILES}
	${FC} -c ${FCFLAGS} $<

utils_measurements.mod: utils_measurements.f90 ${BASEMODFILES}
	${FC} -c ${FCFLAGS} $<

rhmc_updater.mod: rhmc_updater.f90 ${BASEMODFILES}
	${FC} -c ${FCFLAGS} $<

clean:
	echo cleaning up in devicesource.
	$(RM) -f *.o *.a *.mod
