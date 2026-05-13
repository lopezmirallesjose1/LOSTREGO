################################################################################
# LOSTREGO Makefile
################################################################################

#GCC
H5FC = h5pfc
FC = ifort
H5PFC = h5pfc
MPIFC = mpiifort

#XLF
#IFC = xlf
#H5PFC = h5pfc
#FC = h5pfc

#inputs para hacer pruebas
HDF=1
OPENMP=1
MPI=1
DBG=0
MN=0

# Options

ifeq ($(HDF),1)
	UHDF = HDF
        ifeq ($(MPI),1)
                FC = $(H5PFC)
		UPARALELO=PARALELO
	else
                FC = $(H5FC)
		UPARALELO=NOPARALELO	
        endif
else
	UHDF = NOHDF
        ifeq ($(MPI),1)
                FC = $(MPIFC)
		UPARALELO=PARALELO
        else
		ifeq ($(IBM),1)
			FC = $(IFC)
		endif
                UPARALELO=NOPARALELO
        endif
endif


ifeq ($(OPENMP),1)
        UOPENMP = OPENMP
else
        UOPENMP = NOOPENMP
endif

	UCASERO=CASERO
	ifeq ($(OPENMP),1)
           OMPF = -fopenmp
	endif

	ifeq ($(DBG),1)
            FCOPT = -O0 -check all -g -debug all -fpe0 -fpe-all=0 -m64 -implicitnone -ftrapuv -static-libgcc -traceback #-fdump-core -fbacktrace -fcheck=bounds -Wfatal-errors  -Wextra -fdefault-real-8  #-Wall  
        else
            #FCOPT = -m64 -O3 -real-size 64 -fimf-precision=high             
            #FCOPT = -O3 -march=core-avx2 -align array64byte -fma -ftz -fomit-frame-pointer -m64 -real-size 64 -fimf-precision=high 
            FCOPT = -O3 -axSSSE3,CASCADELAKE,COMMON-AVX512,CORE-AVX2,SSE4.2,AVX -align array64byte -fma -ftz -fomit-frame-pointer -m64 -real-size 64 -fimf-precision=high

        endif
        DEFINES =  -D$(UOPENMP) -D$(UHDF) -D$(UPARALELO) -D$(UCASERO) # -DNONSTCPP
        FLAGS = -free #-x f95-cpp-input
        FFLAGS = $(FCOPT) $(OMPF) $(FLAGS) $(DEFINES)


################################################################################

SRC = modulos.F  main3bu.F bndrybu.F eos.F filnam.F geom.F        \
      getprfq3dNR.F gridbu.F mhydrow3.F initbu.F input.F pltout.F   \
      tstep_2D.F tstep_3D.F paralelo.F paraleloCT.F paraleloCTT.F sweeps_2D.F sweeps_3D.F misc.F HLLDD.F plm.F sources.F    \
      restrt_2D.F restrt_3D.F rec_hdf5.F rst_hdf5.F paraleloi.F paraleloj.F paralelok.F sphere.F  

OBJ = modulos.o  main3bu.o bndrybu.o eos.o filnam.o geom.o        \
      getprfq3dNR.o gridbu.o mhydrow3.o initbu.o input.o pltout.o   \
      tstep_2D.o tstep_3D.o paralelo.o paraleloCT.o paraleloCTT.o sweeps_2D.o sweeps_3D.o misc.o HLLDD.o plm.o sources.o    \
      restrt_2D.o restrt_3D.o rec_hdf5.o rst_hdf5.o paraleloi.o paraleloj.o paralelok.o sphere.o

HEADERS = types.h


################################################################################
RATPENAT: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^
#	mv $@ ../bin

%.o %.mod : %.F $(HEADERS)
	$(FC) $(FFLAGS) -c $<

################################################################################
clean:
	rm -rf *.bck
	rm -rf *.lst
	rm -rf *.o
	rm -rf *.mod
	rm -rf *.f
	rm -rf *.err *.out
	rm -rf RATPENAT
