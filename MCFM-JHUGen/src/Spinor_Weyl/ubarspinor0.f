      subroutine ubarspinor0(q,i,f)
c-----subroutine for massless ubar spinor
!     Weyl representation
      implicit none
      include 'swapxz.f'
      integer i
      double complex p(4),q(4),f(4),fc,czip
      double precision  px,py,phase
      parameter(czip=(0d0,0d0))
      logical,save::first
      data first/.true./
     
      if (first) then
      write(6,*) 'ubarspinor0:swapxz=',swapxz
      first=.false.
      endif

C     Translate from MCFM notation 
      p(1)=q(4)
      p(2)=q(3)
      p(3)=-q(2)
      p(4)=q(1)
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

      px=+dreal(p(2))
      py=+dreal(p(3))
c      pz=+dreal(p(4))

      fc=sqrt(p(1)+p(4))
      if (dble(p(1)+p(4)) .gt. 0) then
        phase=1d0
      else
        phase=-1d0
      endif
      


      if (abs(fc).gt.1D-8) then 
      if (i .eq. 1) then 
      f(1)=czip
      f(2)=czip
      f(3)=+fc
      f(4)=dcmplx(px,-py)/fc
      elseif (i .eq. -1) then 
      f(1)=dcmplx(px,py)/fc
      f(2)=-fc
      f(3)=czip
      f(4)=czip
      endif

      else
      if (i .eq. 1) then 
        f(1)=czip
        f(2)=czip
        f(3)=czip
        f(4)=-phase*sqrt(2d0*p(1))
      elseif (i .eq. -1) then 
        f(1)=phase*sqrt(2d0*p(1))
        f(2)=czip
        f(3)=czip
        f(4)=czip
      endif 
      endif 
      return
      end

