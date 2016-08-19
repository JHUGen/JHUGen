      subroutine mcfmsub(r,er)
c--- This is an entry point into MCFM (usually called by mcfm program)    
      implicit none
      double precision r,er
      character*72 inputfile,workdir
      call determinefilenames(inputfile,workdir)
      call mcfmmain(inputfile,workdir,r,er)
      return
      end
