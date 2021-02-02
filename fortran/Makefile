FC=ifort
LDIR=/lib/lapack-3.9.0/
IDIR=./inc/
LIBS= -llapack -lrefblas
OBJDIR=obj/
MYEXEC=harms

ifeq ($(dbg),1)
	FFLAGS = -O0 -g -traceback -check all -check bounds -debug all -gen-interfaces -warn interfaces
else
	FFLAGS = -O3 -funroll-loops -fomit-frame-pointer
endif

MYFUN = UTIL.o \
		HARMONICS.o \
		SHAPE_noali.o \
		main.o \

OBJ = $(patsubst %,$(OBJDIR)%,$(MYFUN))

#hw7: $(OBJ)
	#$(CXX) $(CXXFLAGS) -I$(HEADDIR) -I$(EIGEN) -o $(BINDIR)hw7.exe $(OBJ)
	# mv hw7.exe $(BINDIR)hw7.exe

$(MYEXEC): $(OBJ)
	$(FC) $(FFLAGS) -I$(IDIR) -o $@ $(OBJ) -L$(LDIR)$(LIBS)

UTIL.o : UTIL.f90
	$(FC) $(FFLAGS) -I$(IDIR) -o $@ $(OBJ) -L$(LDIR)$(LIBS)


# $(OBJDIR)%.o : ./%.f90
# 	$(FC) $(FFLAGS) -c -I$(IDIR) $(OBJ) $< -o $@ -L$(LDIR)$(LIBS)

# $(OBJDIR)%.o : $(CPPDIR)%.cpp $(HEADDIR)
# 	$(CXX) -c -I$(HEADDIR) -I$(EIGEN) $(CXXFLAGS) $< -o $@

# $(OBJDIR)hw7.o: $(HWDIR)hw7.cpp
# 	$(CXX) -c -I$(HEADDIR) -I$(EIGEN) $(CXXFLAGS) $< -o $@

.PHONY : clean
clean:
	rm -f $(OBJDIR)*
	rm -f $(MYEXEC)
