      subroutine spinoruorz(N,p,za,zb)
      implicit none
      include 'types.f'
c--- This routine just provides an easy way of switching between
c--- spinoru (normal MCFM running) and spinorz (checks of virtual)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4)
      integer:: N,j,k

      call spinoru(N,p,za,zb)
c      call spinorz(N,p,za,zb)
      
      return
      end
