#     Copyright, 2013
#     Mahdi Esmaily Moghadam

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#--------------------------------------------------------------------
#
# include ../../Makefile.in
FORTRAN=ifort
FFLAGS=-O3
AR 		= ar rv
RANLIB 	= ar -ts

libmpi.a: mpi.o
	$(AR) $@ mpi.o
	$(RANLIB) $@

mpi.o : mpi.f
	$(FORTRAN) $(FFLAGS) -c $< -o $@

clean:	
	rm -r -f libmpi.a mpi.o
