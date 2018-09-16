      function Cgamgam(i1,i2,i3,i4,i5,i)
      implicit none
      include 'types.f'
      real(dp):: Cgamgam
      
CCCCCC Matrix element squared for
C     qbar(-p1)+q(-p2)=gamma(p3)+gamma(p4)+g(p5)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      integer:: i1,i2,i3,i4,i5,i
      real(dp),parameter::statfac=0.5_dp
      Cgamgam=s(i1,i2)
     .*((s(i2,i5)**2+s(i1,i5)**2)/(s(i1,i4)*s(i2,i4)*s(i1,i3)*s(i2,i3))
     & +(s(i2,i4)**2+s(i1,i4)**2)/(s(i1,i5)*s(i2,i5)*s(i1,i3)*s(i2,i3))
     & +(s(i2,i3)**2+s(i1,i3)**2)/(s(i1,i5)*s(i2,i5)*s(i1,i4)*s(i2,i4)))

      Cgamgam=16._dp*cf*xn*gsq*esq**2*Q(i)**4*Cgamgam*statfac

      return
      end
