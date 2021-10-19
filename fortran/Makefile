FC=mpiifort
LDIR=./lib/lapack-3.9.0/
IDIR=./FFTs/
LIBS= -llapack -lrefblas
OBJDIR=obj/
MODDIR=mod/
DATDIR=dat/
DEPS=$(IDIR)fftpack5.o
MPI_DIR=./lib/dummy_mpi

SRC=./src/
MYEXEC=harms

# Normal, opt flags
# Removing fast can make the last line of the compilation faster, whereas other way shape is slow
FFLAGS = -O3 -funroll-loops -fomit-frame-pointer -DNDEBUG #-fast

# Can specify debug flags with make dbg=1
ifeq ($(dbg),1)
	FFLAGS = -O0 -g3 -traceback -check all -check bounds -debug all -gen-interfaces -warn interfaces
endif

# Profiler flags
# -profile-functions -profile-loops=all -profile-loops-report=2
ifeq ($(prof),1)
#	FFLAGS = -prof-gen
	FFLAGS = -prof-use -ipo
endif

# Compile in sequenctial
ifeq ($(seq),1)
	FC = ifort
	MPI_LIBS = -L$(MPI_DIR) -lmpi
	MPI_INCS = -I$(MPI_DIR)
else
	MPI_LIBS =
	MPI_INCS =
endif

# Functions
MYFUN = UTIL.o\
		COMU.o\
		HARMONICS.o\
		SHAPE_interp.o\
		main.o

# Modules
MYMOD = harmmod.mod\
		cmmod.mod\
		utilmod.mod\
		shapemod.mod

OBJ = $(patsubst %,$(OBJDIR)%,$(MYFUN))
MOD = $(patsubst %,$(MODDIR)%,$(MYMOD))

$(OBJDIR)%.o : $(SRC)%.f90 | $(OBJDIR) $(MODDIR) $(DATDIR)
	$(FC) $(FFLAGS) $(MPI_INCS) -module $(MODDIR) -o $@ -c $<

$(MYEXEC): $(OBJ)
	$(FC) $(FFLAGS) -module $(MODDIR) -o $@ $^ $(DEPS) $(MPI_LIBS) -L$(LDIR) $(LIBS)
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

$(OBJDIR)COMU.o : $(SRC)COMU.f90 $(SRC)UTIL.f90

$(OBJDIR)HARMONICS.o : $(SRC)HARMONICS.f90 $(SRC)COMU.f90 $(SRC)UTIL.f90

$(OBJDIR)SHAPE_interp.o : $(SRC)SHAPE_interp.f90 $(SRC)HARMONICS.f90 $(SRC)COMU.f90 $(SRC)UTIL.f90

$(OBJDIR)main.o : $(SRC)main.f90 $(SRC)SHAPE_interp.f90 $(SRC)HARMONICS.f90 $(SRC)COMU.f90 $(SRC)UTIL.f90

.PHONY : clean
clean:
	rm -f $(OBJDIR)*.o
	rm -f $(MODDIR)*.mod
	rm -f $(MYEXEC)

ech:
	@echo "OBJ = $(OBJ)"