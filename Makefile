.PHONY : clean cleanall default

ifeq ($(seq),1)
default:
	(cd fortran/lib ; make ) 
	(cd fortran ; make seq=1 )
	@echo "-----------------------"
	@echo "Compiled successfully!!"
	@echo "-----------------------"
else
default:
	(cd fortran/lib ; make ) 
	(cd fortran ; make ) 
	@echo "-----------------------"
	@echo "Compiled successfully!!"
	@echo "-----------------------"
endif

clean: 
	(cd fortran ; make clean)

cleanall:
	(cd fortran/lib ; make clean)
	(cd fortran ; make clean)
