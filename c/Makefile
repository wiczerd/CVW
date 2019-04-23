FC 	= gfortran
IC	= gcc
FF77	= gfortran

SOLDIR	= /home/david/Computation/DFBOLS/
CDIR	= /home/david/workspace/CVW/c/
FDIR	= /home/david/workspace/CVW/fortran/
MKLROOT = /opt/intel/mkl/

CLIBS	= -lgsl -lgslcblas -lm -lgfortran -lgomp
IADV	= -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread 
ILOC	= -L/usr/local/lib/ -L/home/david/Computation/NLopt/lib -L/opt/intel/lib/intel64/
IINC	= -I/home/david/Computation/DFBOLS


ICFLAGS = -O3 -fopenmp -D_DFBOLS_USE=1
ICGFLAGS= -g -D_DFBOLS_USE=1
FCFLAGS = -c
F77FLAGS= -c

SOLOBJS = altmov.o rescue_h.o update.o prelim_h.o bobyqa_h.o bobyqb_h.o trsbox_h.o dfovec.o


main: $ $(CDIR)main.c $(SOLOBJS) 
	$(IC)  $(ICFLAGS)  $(CDIR)main.c $(SOLOBJS)  $(IINC) $(ILOC) $(CLIBS) -o CVW.out

debug: $ $(CDIR)main.c $(SOLOBJS) 
	$(IC)  $(ICGFLAGS)  $(CDIR)main.c $(SOLOBJS)  $(IINC) $(ILOC) $(CLIBS) -o CVW_dbg.out


# All the solver objects:
altmov.o :  $(SOLDIR)altmov.f
	$(FF77) $(F77FLAGS) $(SOLDIR)altmov.f

rescue_h.o :  $(SOLDIR)rescue_h.f
	$(FF77) $(F77FLAGS) $(SOLDIR)rescue_h.f

update.o :  $(SOLDIR)update.f
	$(FF77) $(F77FLAGS) $(SOLDIR)update.f
	
prelim_h.o :  $(SOLDIR)prelim_h.f
	$(FF77) $(F77FLAGS) $(SOLDIR)prelim_h.f	
	
bobyqa_h.o :  $(SOLDIR)bobyqa_h.f
	$(FF77) $(F77FLAGS) $(SOLDIR)bobyqa_h.f

bobyqb_h.o :  $(SOLDIR)bobyqb_h.f
	$(FF77) $(F77FLAGS) $(SOLDIR)bobyqb_h.f

trsbox_h.o :  $(SOLDIR)trsbox_h.f
	$(FF77) $(F77FLAGS) $(SOLDIR)trsbox_h.f
	
dfovec.o : $(FDIR)dfovec.f95
	$(FC) $(FCFLAGS) $(FDIR)dfovec.f95
