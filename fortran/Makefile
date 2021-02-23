FC=ifort
LDIR=./lib/lapack-3.9.0/
IDIR=./FFTs/
LIBS= -llapack -lrefblas
OBJDIR=obj/
MODDIR=mod/
DATDIR=dat/
DEPS=$(IDIR)fftpack5.o

SRC=./src/
MYEXEC=harms

# Normal, opt flags
FFLAGS = -O3 -funroll-loops -fomit-frame-pointer -DNDEBUG

# Can specify debug flags with make dbg=1
ifeq ($(dbg),1)
	FFLAGS = -O0 -g3 -traceback -check all -check bounds -debug all -gen-interfaces -warn interfaces
endif

# Profiler flags
ifeq ($(prof),1)
	FFLAGS = -profile-functions -profile-loops=all -profile-loops-report=2
endif

# Functions
MYFUN = UTIL.o\
		HARMONICS.o\
		SHAPE_noali.o\
		main.o

# Modules
MYMOD = harmmod.mod\
		utilmod.mod\
		shapemod.mod

OBJ = $(patsubst %,$(OBJDIR)%,$(MYFUN))
MOD = $(patsubst %,$(MODDIR)%,$(MYMOD))

$(OBJDIR)%.o : $(SRC)%.f90 | $(OBJDIR) $(MODDIR) $(DATDIR)
	$(FC) $(FFLAGS) -module $(MODDIR) -o $@ -c $<

$(MYEXEC): $(OBJ)
	$(FC) $(FFLAGS) -module $(MODDIR) -o $@ $^ $(DEPS) -L$(LDIR) $(LIBS)
	@echo "-----------------------"
	@echo "Compiled successfully!!"
	@echo "-----------------------"

$(OBJDIR):
	mkdir -p $@

$(MODDIR):
	mkdir -p $@

$(DATDIR):
	mkdir -p $@

# Dependencies
$(OBJDIR)UTIL.o : $(SRC)UTIL.f90

$(OBJDIR)HARMONICS.o : $(SRC)HARMONICS.f90 $(SRC)UTIL.f90

$(OBJDIR)SHAPE_noali.o : $(SRC)SHAPE_noali.f90 $(SRC)HARMONICS.f90 $(SRC)UTIL.f90

$(OBJDIR)main.o : $(SRC)main.f90 $(SRC)SHAPE_noali.f90 $(SRC)HARMONICS.f90 $(SRC)UTIL.f90

.PHONY : clean
clean:
	rm -f $(OBJDIR)*.o
	rm -f $(MODDIR)*.mod
	rm -f $(MYEXEC)

ech:
	@echo "OBJ = $(OBJ)"