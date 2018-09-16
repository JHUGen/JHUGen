ifeq ($(origin FC),environment)
    $(info Inheriting FC from environment: $(FC))
else
    FC = gfortran
endif
FFLAGS 	= -g

recurdir:
	(cd smallG; make) 
	(cd smallY; make) 
	(cd smallP; make) 
	(cd smallF; make) 

# -----------------------------------------------------------------------------

clean:
	- rm -f *.a */*.o

# -----------------------------------------------------------------------------

# DO NOT DELETE


