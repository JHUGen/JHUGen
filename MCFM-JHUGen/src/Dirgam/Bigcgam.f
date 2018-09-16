      function BigCgam(i1,i2,i3,i4,i5,i)
      implicit none
      include 'types.f'
      real(dp):: BigCgam
      
CCCCCC Matrix element squared for
C     qbar(-p1)+q(-p2)=g(p3)+g(p4)+gamma(p5)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      integer:: i1,i2,i3,i4,i5,i
      BigCgam=
     & (xn*(s(i1,i3)*s(i2,i4)+s(i2,i3)*s(i1,i4))/s(i4,i3)-s(i1,i2)/xn)
     .*((s(i2,i5)**2+s(i1,i5)**2)/(s(i1,i4)*s(i2,i4)*s(i1,i3)*s(i2,i3))
     & +(s(i2,i4)**2+s(i4,i1)**2)/(s(i1,i5)*s(i2,i5)*s(i1,i3)*s(i2,i3))
     & +(s(i2,i3)**2+s(i1,i3)**2)/(s(i1,i5)*s(i2,i5)*s(i1,i4)*s(i2,i4)))

      BigCgam=32._dp*esq*gsq**2*Q(i)**2*BigCgam

      return
      end
