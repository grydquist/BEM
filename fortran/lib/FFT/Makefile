FC	   = ifort
FFLAGS = -O3

fftpack5.o : fftpack5.f90
	$(FC) $(FFLAGS) -o $@ -c $<

.PHONY : clean

clean:
	rm -f fftpack5.o
