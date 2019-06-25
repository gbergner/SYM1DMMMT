
DEPENDALL=staticparameters.f90
include Makefile.inc
BASELIB=-llapack -lblas
ifeq ($(NODEVICE),1)
ifeq ($(MPARALLEL),1)
 HDIR=matrix_parallel
 DIRS=matrix_parallel
 DLIB=
 HLIB=matrix_parallel/libmparallel.a
 EXE=bfss_mparallel bfss_mparallel_input
else
 HDIR=hostcode
 DIRS=hostcode
 DLIB=
 HLIB=hostcode/libhost.a $(BASELIB)
 EXE=bfss_host_serial
endif
else
 HDIR=hostcode
 DIRS=devicecode hostcode
 DLIB=devicecode/libdevice.a
 HLIB=hostcode/libhost.a $(BASELIB)
 EXE=testprogram bfss_host_serial_read test_timing bfss_host_serial bfss_device_serial test_timing_inverter bfss_device_serial_hyst
endif

all: $(EXE)

${DLIB}: force_look
	echo looking into devicecode : $(MAKE) $(MFLAGS)
	cd devicecode; $(MAKE) $(MFLAGS)

${HLIB}: force_look
	echo looking into $(HDIR) : $(MAKE) $(MFLAGS)
	cd $(HDIR); $(MAKE) $(MFLAGS)

compare_host_device.mod: compare_host_device.f90
	${FC} -c ${FCFLAGS} $<

bfss_host_serial: ${HLIB} main_host_serial.o 
	${FC} ${FCFLAGS} main_host_serial.o ${HLIB} -o $@

bfss_host_serial_read: ${HLIB} main_host_serial_readintermediate.o 
	${FC} ${FCFLAGS} main_host_serial_readintermediate.o ${HLIB} -o $@

bfss_device_serial: ${HLIB} ${DLIB} measure_host_device.o main_device_serial.o 
	${FC} ${FCFLAGS} main_device_serial.o measure_host_device.o ${DLIB} ${HLIB}  -o $@

bfss_device_serial_hyst: ${HLIB} ${DLIB} measure_host_device.o main_device_serial_hyst.o 
	${FC} ${FCFLAGS} main_device_serial_hyst.o measure_host_device.o ${DLIB} ${HLIB}  -o $@

testprogram:  ${HLIB} ${DLIB} compare_host_device.mod main_testprogram.o compare_host_device.o
	${FC} ${FCFLAGS}  main_testprogram.o compare_host_device.o ${HLIB} ${DLIB}  -o $@

test_timing: main_time.o ${HLIB} ${DLIB}  
	${FC} ${FCFLAGS} main_time.o  ${HLIB} ${DLIB} -o $@

test_timing_inverter: main_time_inverter.o ${HLIB} ${DLIB}  
	${FC} ${FCFLAGS} main_time_inverter.o  ${HLIB} ${DLIB} -o $@

bfss_mparallel: ${HLIB} main_parallel.o
	${FC} ${FCFLAGS} main_parallel.o ${DLIB} ${HLIB}  -o $@

bfss_mparallel_input: ${HLIB} main_parallel_input.o
	${FC} ${FCFLAGS} main_parallel_input.o ${DLIB} ${HLIB}  -o $@

.PHONY: clean

clean:
	echo cleaning up in .
	$(RM) -f $(EXE) *.o *.mod
	for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true


