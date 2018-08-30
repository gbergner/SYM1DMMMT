 DIRS=devicecode hostcode
EXE=testprogram test_timing
DLIB=devicecode/libdevice.a
HLIB=hostcode/libhost.a
EXE=testprogram test_timing bfss_host_serial bfss_device_serial
DEPENDALL=staticparameters.f90
include Makefile.inc

all: $(EXE)

${DLIB}: force_look
	echo looking into devicecode : $(MAKE) $(MFLAGS)
	cd devicecode; $(MAKE) $(MFLAGS)

${HLIB}: force_look
	echo looking into hostcode : $(MAKE) $(MFLAGS)
	cd hostcode; $(MAKE) $(MFLAGS)
	
compare_host_device.mod: compare_host_device.f90
	${FC} -c ${FCFLAGS} $<
	
bfss_host_serial: ${HLIB} ${DLIB} main_host_serial.o 
	${FC} ${FCFLAGS} main_host_serial.o ${HLIB} ${DLIB}  -o $@
	
bfss_device_serial: ${HLIB} ${DLIB} main_device_serial.o 
	${FC} ${FCFLAGS} main_device_serial.o ${DLIB} ${HLIB}  -o $@

testprogram:  ${HLIB} ${DLIB} compare_host_device.mod main_testprogram.o compare_host_device.o
	${FC} ${FCFLAGS}  main_testprogram.o compare_host_device.o ${HLIB} ${DLIB}  -o $@

test_timing: main_time.o ${HLIB} ${DLIB}  
	${FC} ${FCFLAGS} main_time.o  ${HLIB} ${DLIB} -o $@

.PHONY: clean

clean:
	echo cleaning up in .
	$(RM) -f $(EXE) *.o *.mod
	for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true


