      subroutine uspinor0(q,i,f)
!-----subroutine for massless spinor
!     Majorana representation
      implicit none
      include 'constants.f'
      include 'swapxz.f'
      integer i
      double complex p(4),q(4),f(4),fc
      double precision  Ep,px,py,pz,phase,rtEon2
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
      Ep=dreal(p(1))
      px=+dreal(p(2))
      py=+dreal(p(3))
      pz=dreal(p(4))

      fc=sqrt(p(1)+p(4))
      if (dble(p(1)+p(4)) .gt. 0) then
        phase=1d0
      else
        phase=-1d0
      endif

      if (abs(fc).gt.1D-8) then 

      if (i.eq.+1) then 
        f(1)=dcmplx(pz+py+Ep,-px)/(2d0*fc)
        f(2)=+im*f(1)
        f(3)=-dcmplx(pz-py+Ep,+px)/(2d0*fc)
        f(4)=-im*f(3)
      elseif (i.eq.-1) then 
        f(1)=-dcmplx(-px,pz+py+Ep)/(2d0*fc)
        f(2)=-im*f(1)
        f(3)=+dcmplx(px,pz-py+Ep)/(2d0*fc)
        f(4)=+im*f(3)
      endif 
      else
      rtEon2=sqrt(Ep/2d0)
      if (i.eq.+1) then 
        f(1)=dcmplx(0d0,-rtEon2)
        f(2)=+im*f(1)
        f(3)=dcmplx(0d0,-rtEon2)
        f(4)=-im*f(3)
      elseif (i.eq.-1) then 
        f(1)=dcmplx(rtEon2,0d0)
        f(2)=-im*f(1)
        f(3)=dcmplx(rtEon2,0d0)
        f(4)=+im*f(3)
      endif 
      endif 
         
      return
      end

