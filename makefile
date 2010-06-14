# Makefile for s2dw library
# Jason McEwen


# ======== COMPILER ========

#FC      = f95
FC      = gfortran
#FC      = g95

ifeq ($(FC),f95)
  OPTF95 = -w=x95
endif

OPT = $(OPTF95)


# ======== LINKS ========

PROGDIR = /Users/jasonmcewen/Prog

S2DWDIR  = $(PROGDIR)/s2dw-1.0
S2DWLIB  = $(S2DWDIR)/lib
S2DWLIBNM= s2dw
S2DWINC  = $(S2DWDIR)/include
S2DWSRC  = $(S2DWDIR)/src/mod
S2DWPROG = $(S2DWDIR)/src/prog
S2DWBIN  = $(S2DWDIR)/bin
S2DWDOC  = $(S2DWDIR)/doc

FFTWLIB      = $(PROGDIR)/fftw-3.1.2/lib
FFTWLIBNM    = fftw3

HPIXDIR = $(PROGDIR)/Healpix_2.01
HPIXLIB = $(HPIXDIR)/lib
HPIXLIBNM= healpix
HPIXINC = $(HPIXDIR)/include

S2DIR  = $(PROGDIR)/s2-1.0
S2LIB  = $(S2DIR)/lib
S2LIBNM= s2
S2INC  = $(S2DIR)/include
S2DOC  = $(S2DIR)/doc

CFITSIOLIB   = $(PROGDIR)/cfitsio/lib
CFITSIOLIBNM = cfitsio


# ======== FFFLAGS ========

FFLAGS  = -I$(S2DWINC)
FFLAGSPROG = -I$(HPIXINC) -I$(S2INC)


# ======== LDFLAGS ========

LDFLAGS = -L$(S2DWLIB) -l$(S2DWLIBNM) \
          -L$(FFTWLIB) -l$(FFTWLIBNM) \
          -L$(CFITSIOLIB) -l$(CFITSIOLIBNM) 

LDFLAGSPROG = -L$(S2LIB) -l$(S2LIBNM) \
           -L$(HPIXLIB) -l$(HPIXLIBNM) 


# ======== PPFLAGS ========

ifeq ($(FC),f95)
  PPFLAGS = -fpp $(OPT)
else ifeq ($(FC),g95)
  PPFLAGS = -cpp $(OPT)
else ifeq ($(FC),gfortran)
  PPFLAGS = -x f95-cpp-input $(OPT)
endif


# ======== OBJECT FILES TO MAKE ========

S2DWOBJ = $(S2DWINC)/s2dw_types_mod.o    \
          $(S2DWINC)/s2dw_error_mod.o   \
          $(S2DWINC)/s2dw_dl_mod.o      \
          $(S2DWINC)/s2dw_fileio_mod.o  \
          $(S2DWINC)/s2dw_core_mod.o   


# ======== MAKE RULES ========

default: lib

all:     lib prog test

lib:	 $(S2DWLIB)/lib$(S2DWLIBNM).a

test:    $(S2DWBIN)/s2dw_test

runtest: test
	$(S2DWBIN)/s2dw_test 64

prog:    $(S2DWBIN)/s2dw_wav2sky $(S2DWBIN)/s2dw_analysis $(S2DWBIN)/s2dw_synthesis $(S2DWBIN)/s2dw_wavplot $(S2DWBIN)/s2dw_mat2fits $(S2DWBIN)/s2dw_fits2mat

$(S2DWINC)/%.o: $(S2DWSRC)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 
	mv *.mod $(S2DWINC)

$(S2DWINC)/s2dw_test.o:     $(S2DWPROG)/s2dw_test.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 

$(S2DWINC)/%.o: $(S2DWPROG)/%.f90
	$(FC) $(FFLAGS) $(FFLAGSPROG) $(PPFLAGS) -I$(S2INC) -c $< -o $@ 


# Library

$(S2DWLIB)/lib$(S2DWLIBNM).a: $(S2DWOBJ)
	ar -r $(S2DWLIB)/lib$(S2DWLIBNM).a $(S2DWOBJ)


# Documentation

docs:
	./f90doc_fpp $(S2DWSRC)/*.f90
	./f90doc_fpp $(S2DWPROG)/*.f90
	./ln_multi $(S2DOC)/s2_*
	./ln_multi $(S2DOC)/index_s2.html
	mv *.html $(S2DWDOC)/.
	./addstyle $(S2DWDOC)/s2dw_*

cleandocs:
	rm -f $(S2DWDOC)/s2dw_*.html
	rm -f $(S2DWDOC)/gasdev2_dp.html $(S2DWDOC)/ran2_dp.html
	rm -f $(S2DWDOC)/s2_*.html $(S2DWDOC)/index_s2.html

# Cleaning up

clean:	tidy
	rm -f $(S2DWINC)/*.mod
	rm -f $(S2DWINC)/*.o
	rm -f $(S2DWLIB)/lib$(S2DWLIBNM).a
	rm -f $(S2DWBIN)/*

tidy:
	rm -f *.mod
	rm -f $(S2DWSRC)/*~ 
	rm -f $(S2DWPROG)/*~ 


# Module dependencies

$(S2DWINC)/s2dw_types_mod.o: $(S2DWSRC)/s2dw_types_mod.f90
$(S2DWINC)/s2dw_error_mod.o: $(S2DWSRC)/s2dw_error_mod.f90  \
                           $(S2DWINC)/s2dw_types_mod.o
$(S2DWINC)/s2dw_dl_mod.o:    $(S2DWSRC)/s2dw_dl_mod.f90     \
                           $(S2DWINC)/s2dw_types_mod.o
$(S2DWINC)/s2dw_fileio_mod.o:  $(S2DWSRC)/s2dw_fileio_mod.f90   \
                           $(S2DWINC)/s2dw_types_mod.o    \
                           $(S2DWINC)/s2dw_error_mod.o    \
                           $(S2DWINC)/s2dw_core_mod.o    
$(S2DWINC)/s2dw_core_mod.o:  $(S2DWSRC)/s2dw_core_mod.f90   \
                           $(S2DWINC)/s2dw_types_mod.o    \
                           $(S2DWINC)/s2dw_error_mod.o    \
                           $(S2DWINC)/s2dw_dl_mod.o       


# Program dependencies and compilation

$(S2DWINC)/s2dw_test.o:     $(S2DWPROG)/s2dw_test.f90 lib
$(S2DWBIN)/s2dw_test:       $(S2DWINC)/s2dw_test.o
	$(FC)                                          \
	-o $(S2DWBIN)/s2dw_test                          \
	$(S2DWINC)/s2dw_test.o $(LDFLAGS)                            

$(S2DWINC)/s2dw_wav2sky.o:     $(S2DWPROG)/s2dw_wav2sky.f90 lib
$(S2DWBIN)/s2dw_wav2sky:       $(S2DWINC)/s2dw_wav2sky.o
	$(FC)                                          \
	-o $(S2DWBIN)/s2dw_wav2sky                       \
	$(S2DWINC)/s2dw_wav2sky.o $(LDFLAGSPROG) $(LDFLAGS)                           

$(S2DWINC)/s2dw_analysis.o:     $(S2DWPROG)/s2dw_analysis.f90 lib
$(S2DWBIN)/s2dw_analysis:       $(S2DWINC)/s2dw_analysis.o
	$(FC)                                          \
	 -o $(S2DWBIN)/s2dw_analysis                     \
	$(S2DWINC)/s2dw_analysis.o $(LDFLAGSPROG) $(LDFLAGS)                         

$(S2DWINC)/s2dw_synthesis.o:     $(S2DWPROG)/s2dw_synthesis.f90 lib
$(S2DWBIN)/s2dw_synthesis:       $(S2DWINC)/s2dw_synthesis.o
	$(FC)                                          \
	-o $(S2DWBIN)/s2dw_synthesis                     \
	$(S2DWINC)/s2dw_synthesis.o $(LDFLAGSPROG) $(LDFLAGS)                          

$(S2DWINC)/s2dw_wavplot.o:     $(S2DWPROG)/s2dw_wavplot.f90 lib
$(S2DWBIN)/s2dw_wavplot:       $(S2DWINC)/s2dw_wavplot.o
	$(FC)                                          \
	-o $(S2DWBIN)/s2dw_wavplot                       \
	$(S2DWINC)/s2dw_wavplot.o $(LDFLAGSPROG) $(LDFLAGS)                          

$(S2DWINC)/s2dw_mat2fits.o:     $(S2DWPROG)/s2dw_mat2fits.f90 lib
$(S2DWBIN)/s2dw_mat2fits:       $(S2DWINC)/s2dw_mat2fits.o
	$(FC)                                          \
	-o $(S2DWBIN)/s2dw_mat2fits                       \
	$(S2DWINC)/s2dw_mat2fits.o $(LDFLAGSPROG) $(LDFLAGS)  

$(S2DWINC)/s2dw_fits2mat.o:     $(S2DWPROG)/s2dw_fits2mat.f90 lib
$(S2DWBIN)/s2dw_fits2mat:       $(S2DWINC)/s2dw_fits2mat.o
	$(FC)                                          \
	-o $(S2DWBIN)/s2dw_fits2mat                       \
	$(S2DWINC)/s2dw_fits2mat.o $(LDFLAGSPROG) $(LDFLAGS) 


