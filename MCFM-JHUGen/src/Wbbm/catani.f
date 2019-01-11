      subroutine catani(mhsq,k1,k2,k3,k4,k5,k6,coeff2,coeff1)
      implicit none
      include 'constants.f'
      include 'momwbbm.f'
      include 'scale.f'
      integer k1,k2,k3,k4,k5,k6,nu
      double precision mhsq,p2(4),p3(4),p1Dp2,p3Dp4
      double complex xl12,xl34,coeff2,coeff1
      
      do nu=1,4
      p2(nu)=bp*mom(k2,nu)+bm*mom(k3,nu)
      p3(nu)=bp*mom(k3,nu)+bm*mom(k2,nu)
      enddo
      p1Dp2=mom(k1,4)*p2(4)-mom(k1,1)*p2(1)
     &     -mom(k1,2)*p2(2)-mom(k1,3)*p2(3)
      p3Dp4=mom(k4,4)*p3(4)-mom(k4,1)*p3(1)
     &     -mom(k4,2)*p3(2)-mom(k4,3)*p3(3)
      xl12=dcmplx(dlog(abs(2d0*p1Dp2)/musq))
      xl34=dcmplx(dlog(abs(2d0*p3Dp4)/musq))
      
      coeff2=-cone
      coeff1=xl12+xl34+dcmplx(8d0/3d0-dlog(mhsq/musq))
      
      return
      end


c--- Note: this is not Catani's formula, this is my reinterpretation 
c---       of the result that we have found      
      subroutine catanisl(mhsq,k1,k2,k3,k4,k5,k6,coeff2,coeff1)
      implicit none
      include 'constants.f'
      include 'momwbbm.f'
      include 'scale.f'
      integer k1,k2,k3,k4,k5,k6,nu
      double precision mhsq,p2(4),p3(4),be,p2Dp3,p1Dp4
      double complex xl14,coeff2,coeff1,Lnrat
      
      do nu=1,4
      p2(nu)=bp*mom(k2,nu)+bm*mom(k3,nu)
      p3(nu)=bp*mom(k3,nu)+bm*mom(k2,nu)
      enddo
      p2Dp3=p2(4)*p3(4)-p2(1)*p3(1)
     &     -p2(2)*p3(2)-p2(3)*p3(3)
      p1Dp4=mom(k1,4)*mom(k4,4)-mom(k1,1)*mom(k4,1)
     &     -mom(k1,2)*mom(k4,2)-mom(k1,3)*mom(k4,3)
      xl14=Lnrat(musq,-2d0*p1Dp4)

      coeff2=-cone
      be=bp-bm
      coeff1=-cone-xl14-dcmplx(0.5d0*(1d0+be**2)/be)
     & *Lnrat(mhsq,-2d0*(p2Dp3+mhsq)*bp**2)
      
      return
      end
     
 
