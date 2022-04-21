FC=mpiifort
LIBDIR= ./lib/
LA_LDIR= $(LIBDIR)lapack/
IDIR= $(LIBDIR)FFT/
LA_LIBS= -llapack
OBJDIR=obj/
MODDIR=mod/
DATDIR=dat/
DEPS=$(IDIR)fftpack5.o
MPI_DIR=$(LIBDIR)dummy_mpi
SRC=./src/
MYEXEC=harms

# Normal, opt flags
# Removing fast can make the last line of the compilation faster, whereas other way shape is slow
FFLAGS = -O3 -funroll-loops -fomit-frame-pointer -DNDEBUG #-fast -traceback 

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

# Compile in sequential
ifeq ($(seq),1)
	FC = ifort
	MPI_LIBS = -L$(MPI_DIR) -lmpi
else
	MPI_LIBS =
endif

# Functions
MYFUN = UTIL.o\
		COMU.o\
		HARMONICS.o\
		SHARED.o\
		SHAPE.o\
		PERIODIC.o\
		PROBLEM.o\
		main.o

# Modules
MYMOD = harmmod.mod\
		cmmod.mod\
		utilmod.mod\
		sharedmod.mod\
		periodicmod.mod\
		probmod.mod\
		shapemod.mod

OBJ = $(patsubst %,$(OBJDIR)%,$(MYFUN))
MOD = $(patsubst %,$(MODDIR)%,$(MYMOD))

$(OBJDIR)%.o : $(SRC)%.f90 | $(OBJDIR) $(MODDIR) $(DATDIR)
	$(FC) $(FFLAGS) -module $(MODDIR) -o $@ -c $<

$(MYEXEC): $(OBJ)
	$(FC) $(FFLAGS) -module $(MODDIR) -o $@ $^ $(DEPS) $(MPI_LIBS) -L$(LA_LDIR) $(LA_LIBS)
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

$(OBJDIR)COMU.o : $(SRC)COMU.f90

$(OBJDIR)HARMONICS.o : $(SRC)HARMONICS.f90 $(OBJDIR)UTIL.o

$(OBJDIR)SHARED.o : $(SRC)SHARED.f90 $(OBJDIR)HARMONICS.o

$(OBJDIR)SHAPE.o : $(SRC)SHAPE.f90 $(OBJDIR)SHARED.o

$(OBJDIR)PERIODIC.o : $(SRC)PERIODIC.f90 $(OBJDIR)SHARED.o

$(OBJDIR)PROBLEM.o : $(SRC)PROBLEM.f90 $(OBJDIR)SHAPE.o $(OBJDIR)COMU.o  $(OBJDIR)PERIODIC.o

$(OBJDIR)main.o : $(SRC)main.f90 $(OBJDIR)PROBLEM.o 

.PHONY : clean
clean:
	rm -f $(OBJDIR)*.o
	rm -f $(MODDIR)*.mod
	rm -f $(MYEXEC)

ech:
	@echo "OBJ = $(OBJ)"