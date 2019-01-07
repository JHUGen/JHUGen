FC = gfortran
#FFLAGS 	= -g
FFLAGS 	=  -fno-f2c -fPIC -g -O2

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


