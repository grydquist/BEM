.PHONY: clean

default:
	(cd lapack ; make )
	(cd dummy_mpi ; make )
	(cd FFT ; make )
	(cd scalapack ; make )
	@echo "-----------------------"
	@echo "Compiled successfully!!"
	@echo "-----------------------"

clean:
	(cd lapack ; make clean)
	(cd dummy_mpi ; make clean)
	(cd FFT ; make clean)
	(cd scalapack ; make clean)