      subroutine uspinor0(q,i,f)
      implicit none
      include 'types.f'
      include 'cplx.h'
!-----subroutine for massless spinor
!     Weyl representation
      
      include 'swapxz.f'
      integer:: i
      complex(dp):: p(4),q(4),f(4),fc,czip
      real(dp)::  px,py,phase
      parameter(czip=(0._dp,0._dp))
      logical,save::first
      data first/.true./
     
      if (first) then
      write(6,*) 'uspinor0:swapxz=',swapxz
      first=.false.
      endif
 
C     Translate from MCFM notation (E=q(4)), 
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
      px=+real(p(2))
      py=+real(p(3))

      fc=sqrt(p(1)+p(4))
      if (real(p(1)+p(4)) > 0) then
        phase=1._dp
      else
        phase=-1._dp
      endif

      if (abs(fc)>1.e-8_dp) then 

      if (i==+1) then 
        f(1)=+fc
        f(2)=cplx2(px,py)/fc
        f(3)=czip
        f(4)=czip
      elseif (i==-1) then 
        f(1)=czip
        f(2)=czip
        f(3)=+cplx2(px,-py)/fc
        f(4)=-fc
      endif 
      else
      if (i==+1) then 
        f(1)=czip
        f(2)=-phase*sqrt(2._dp*p(1))
        f(3)=czip
        f(4)=czip
      elseif (i==-1) then 
        f(1)=czip
        f(2)=czip
        f(3)=phase*sqrt(2._dp*p(1))
        f(4)=czip
      endif 
      endif 
         
      return
      end

