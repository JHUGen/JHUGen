      subroutine gen2a(r,p,wt2,*)
c----1+2 --> 3+4
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'phasemin.f'
      include 'x1x2.f'
      integer nu
      double precision r(mxdim),p(mxpart,4),wt2
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      include 'energy.f'

      p(:,:)=0d0
      xx(1)=xmin+(1d0-xmin)*r(1)
      xx(2)=xmin+(1d0-xmin)*r(2)

      p1(4)=-0.5d0*xx(1)*sqrts
      p1(1)=0d0
      p1(2)=0d0
      p1(3)=-0.5d0*xx(1)*sqrts
      
      p2(4)=-0.5d0*xx(2)*sqrts
      p2(1)=0d0
      p2(2)=0d0
      p2(3)=+0.5d0*xx(2)*sqrts

      call phase2(r,p1,p2,p3,p4,wt2,*999)

c---5,6,7=dummies -- just so nothing crashes
      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      enddo
      wt2=wt2*(1d0-xmin)**2
      return

 999  continue
      wt2=0d0
      p(:,:)=0d0
      return
      end
