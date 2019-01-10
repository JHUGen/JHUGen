      subroutine spinoruorz(N,p,za,zb)
c--- This routine just provides an easy way of switching between
c--- spinoru (normal MCFM running) and spinorz (checks of virtual)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      double precision p(mxpart,4)
      integer N,j,k

      call spinoru(N,p,za,zb)
c      call spinorz(N,p,za,zb)
      
      return
      end
