      subroutine jzero(n2,n1,zab,zba,j0)
      implicit none
      include 'constants.f'
      include 'spinzerohiggs_anomcoupl.f'
      double complex zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),j0(4,2)
C---The one Z-current multiplied by i
C---order of indices Lorentz,quark-line helicity

      integer n1,n2,nu
      do nu=1,4
      j0(nu,1)=alpha_SMEW*zab(n2,nu,n1)
      j0(nu,2)=alpha_SMEW*zba(n2,nu,n1)
      enddo
      return
      end
