      subroutine ubarspinor0(q,i,f)
      implicit none
      include 'types.f'
c-----subroutine for massless ubar spinor
!     Majorana representation
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'swapxz.f'
      integer:: i
      complex(dp):: p(4),q(4),f(4),fc
      real(dp)::  Ep,px,py,pz,phase,rtEon2
      logical,save::first
      data first/.true./
     
      if (first) then
      write(6,*) 'ubarspinor0:swapxz=',swapxz
      first=.false.
      endif

      if (swapxz) then
C performing the swap (x<->z),(y->-y)
      p(1)=q(4)
      p(2)=q(3)
      p(3)=-q(2)
      p(4)=q(1)
      else
      p(1)=q(4)
      p(2)=q(1)
      p(3)=q(2)
      p(4)=q(3)
      endif

      Ep=+real(p(1))
      px=+real(p(2))
      py=+real(p(3))
      pz=+real(p(4))

      fc=sqrt(p(1)+p(4))
      if (real(p(1)+p(4)) > 0) then
        phase=1._dp
      else
        phase=-1._dp
      endif
      
      if (abs(fc)>1.e-8_dp) then 
      if (i == 1) then 
      f(1)=cplx2(Ep+pz-py,-px)/(2._dp*fc)
      f(2)=+im*f(1)
      f(3)=cplx2(Ep+pz+py,px)/(2._dp*fc)
      f(4)=-im*f(3)
      elseif (i == -1) then 
      f(1)=cplx2(px,-Ep-pz+py)/(2._dp*fc)
      f(2)=-im*f(1)
      f(3)=cplx2(-px,-Ep-pz-py)/(2._dp*fc)
      f(4)=+im*f(3)
      endif

      else
      rtEon2=sqrt(Ep/2._dp)
      if (i == 1) then 
        f(1)=cplx2(0._dp,-rtEon2)
        f(2)=+im*f(1)
        f(3)=cplx2(0._dp,+rtEon2)
        f(4)=-im*f(3)
       elseif (i == -1) then 
        f(1)=cplx2(rtEon2,0._dp)
        f(2)=-im*f(1)
        f(3)=cplx2(-rtEon2,0._dp)
        f(4)=+im*f(3)
      endif 
      endif 
      return
      end

