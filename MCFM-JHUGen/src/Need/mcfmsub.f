      subroutine mcfmsub(r,er)
      implicit none
      include 'types.f'
c--- This is an entry point into MCFM (usually called by mcfm program)

      real(dp):: r,er
      character*72 inputfile,workdir
      call determinefilenames(inputfile,workdir)
      call mcfmmain(inputfile,workdir,r,er)
      return
      end
