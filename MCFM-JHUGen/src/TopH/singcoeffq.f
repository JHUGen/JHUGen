      subroutine singcoeffq(p,j1,j2,j3,j4,cnab,cqed)
      implicit none
      include 'types.f'
      
C     Singular parts for process t(j2) t~(j1) q(p3) q~(p4)+H
C     cnab is the coefficient of d_(i2,i4)*d_(i3,i1)
C     cqed is the coefficient of d_(i2,i1)*d_(i3,i4)/xn
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'scale.f'
      real(dp):: p(mxpart,4),v12,htheta,
     & x,s12,s13,s14,s23,s24,s34,
     & xlog12,xl13,xl14,xl23,xl24,xl34,xlm13,xlm14,xlm23,xlm24,xlm34
      complex(dp):: cnab(-2:-1),cqed(-2:-1)
      integer:: j1,j2,j3,j4
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+half*sign(one,x)
      
      s12=2d0
     &*(p(j1,4)*p(j2,4)-p(j1,1)*p(j2,1)-p(j1,2)*p(j2,2)-p(j1,3)*p(j2,3))
      s13=2d0
     &*(p(j1,4)*p(j3,4)-p(j1,1)*p(j3,1)-p(j1,2)*p(j3,2)-p(j1,3)*p(j3,3))
      s14=2d0
     &*(p(j1,4)*p(j4,4)-p(j1,1)*p(j4,1)-p(j1,2)*p(j4,2)-p(j1,3)*p(j4,3))
      s23=2d0
     &*(p(j2,4)*p(j3,4)-p(j2,1)*p(j3,1)-p(j2,2)*p(j3,2)-p(j2,3)*p(j3,3))
      s24=2d0
     &*(p(j2,4)*p(j4,4)-p(j2,1)*p(j4,1)-p(j2,2)*p(j4,2)-p(j2,3)*p(j4,3))
      s34=2d0
     &*(p(j3,4)*p(j4,4)-p(j3,1)*p(j4,1)-p(j3,2)*p(j4,2)-p(j3,3)*p(j4,3))
      v12=sqrt(1d0-4d0*mt**4/s12**2)
      xlog12=log((1d0-v12)/(1d0+v12))/v12
      xlm23=log(mt**2/abs(s23))
      xlm14=log(mt**2/abs(s14))
      xl14=log(musq/abs(s14))
      xl13=log(musq/abs(s13))
      xl23=log(musq/abs(s23))
      xl24=log(musq/abs(s24))
      xl34=log(musq/abs(s34))
      
      cnab(-2)= + xn**(-1) * ( 1.D0 )
      cnab(-2) = cnab(-2) + xn * (  - 1.D0 )

      cnab(-1)= + xn**(-1) * ( 5.D0/2.D0 + xlm24 - xlm23 - xlm14 + 
     &    xlm13 + xl34 + xl24 - xl23 - xl14 + xl13 + 1.D0/2.D0*xlog12 )
      cnab(-1) = cnab(-1) + impi * ( 2.D0*htheta(s12)*v12**(-1) + 2.D0*
     &    htheta(s13) + 2.D0*htheta(s14) + 2.D0*htheta(s23) + 2.D0*
     &    htheta(s24) + 2.D0*htheta(s34) )
      cnab(-1) = cnab(-1) + xn * ( 7.D0/6.D0 - 1.D0/2.D0*xlm24 - 1.D0/2.
     &    D0*xlm13 - 1.D0/2.D0*xl24 - 1.D0/2.D0*xl13 )
      cnab(-1) = cnab(-1) - 2.D0/3.D0*nf

      cqed(-2)= + xn**(-2) * (  - 1.D0 )
      cqed(-2) = cqed(-2) + 1.D0

      cqed(-1)= + xn**(-2) * (  - 5.D0/2.D0 - 1.D0/2.D0*xlm23 - 1.D0/2.D
     &    0*xlm14 - 1.D0/2.D0*xl23 - 1.D0/2.D0*xl14 )
      cqed(-1) = cqed(-1) + xn**(-1) * (  - 1.D0/2.D0*xlm24 + xlm23 + 
     &    xlm14 - 1.D0/2.D0*xlm13 - xl34 - 1.D0/2.D0*xl24 + xl23 + xl14
     &     - 1.D0/2.D0*xl13 + 2.D0/3.D0*nf - 1.D0/2.D0*xlog12 )
      cqed(-1) = cqed(-1) + xn**(-1)*impi * (  - 2.D0*htheta(s12)*
     &    v12**(-1) - 2.D0*htheta(s13) - 2.D0*htheta(s14) - 2.D0*
     &    htheta(s23) - 2.D0*htheta(s24) - 2.D0*htheta(s34) )
      cqed(-1) = cqed(-1) - 7.D0/6.D0 + 1.D0/2.D0*xlm23 + 1.D0/2.D0*
     &    xlm14 + 1.D0/2.D0*xl23 + 1.D0/2.D0*xl14

 
      return
      end
