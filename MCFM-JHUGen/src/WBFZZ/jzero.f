      subroutine jzero(n2,n1,zab,zba,j0)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      complex(dp):: zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & j0(4,2)
C---The one Z-current multiplied by i
C---order of indices Lorentz,quark-line helicity

      integer:: n1,n2,nu
      do nu=1,4
      j0(nu,1)=zab(n2,nu,n1)
      j0(nu,2)=zba(n2,nu,n1)
      enddo
      return
      end
