FC= mpiifort
# --- RELEASE MODE ----
# ON ITASCA, LOAD INTEL/2013.5 IMPI/INTEL MKL/11.0.5.192
# ON MESABI, LOAD INTEL/2013.5 IMPI/INTEL MKL/11.0.5.192
# MAYBE ON MESABI, LOAD INTEL/2015/UPDATE2 IMPI/INTEL 
# --- RELEAE MODE OPTION 1 GOOD
# FCFLAGS = -g -O3 -openmp -xCORE-AVX2 -mt_mpi -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include   
# --- RELEAE MODE OPTION 2 GOOD    
# FCFLAGS = -g -O3 -openmp -axAVX -mt_mpi -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include       

# --- DEBUG MODE ---
FCFLAGS = -g -check all -traceback -check noarg_temp_created -openmp -axAVX -mt_mpi -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include  
     
LDFLAGS = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
      $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a \
      $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a \
      -Wl,--start-group \
      $(MKLROOT)/lib/intel64/libmkl_cdft_core.a \
      $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
      $(MKLROOT)/lib/intel64/libmkl_core.a \
      $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
      -Wl,--end-group \
      $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
      -lpthread -lm
PROGRAM = main #executable filename      

# "make" builds all       
all: ${PROGRAM}
#	OMP_NUM_THREADS=24 mpirun -perhost 1 -np 1 ./main | tee output.txt
#	./main | tee output.txt
main: equilibrium.o model.o variable.o toolbox.o
#main: variable.o toolbox.o
main.o: equilibrium.o
#main.o: variable.o
equilibrium.o: model.o
model.o: variable.o
variable.o: toolbox.o

# general rule
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS) 
%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(LDFLAGS)
#-I$(MODDIR)

# utility targets	
.PHONY: clean 
clean:
	find . -maxdepth 1 -type f -name "*.txt" ! -name "_*.txt" -exec rm '{}' \;
# rm -f ${PROGRAM}
	rm -f *.o
	rm -f *.mod
	rm -f *~
	rm -f *.gz
	rm -f fort.*
	rm main


