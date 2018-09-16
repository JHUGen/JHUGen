      function Bigagam(i1,i2,i3,i4,i5,i,j)
      implicit none
      include 'types.f'
      real(dp):: Bigagam
      
C     q_i(p1)+q_j(p2)->q_i(p3)+q_j(p4)+gamma(p5)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      integer:: i1,i2,i3,i4,i5,i,j
      real(dp):: ss,ssp,tt,ttp,uu,uup,e12,e34,e13,e14,e23,e24,eik
      e12=s(i1,i2)/(s(i1,i5)*s(i2,i5))
      e34=s(i3,i4)/(s(i3,i5)*s(i4,i5))
      e14=s(i1,i4)/(s(i1,i5)*s(i4,i5))
      e24=s(i2,i4)/(s(i2,i5)*s(i4,i5))
      e13=s(i1,i3)/(s(i1,i5)*s(i3,i5))
      e23=s(i2,i3)/(s(i2,i5)*s(i3,i5))
      eik=Q(i)*Q(j)*(-e12+e23+e14-e34)
     & +Q(i)**2*e13+Q(j)**2*e24
      ss=s(i1,i2)
      ssp=s(i3,i4)
      tt=s(i1,i3)
      ttp=s(i2,i4)
      uu=s(i1,i4)
      uup=s(i2,i3)
      Bigagam=8._dp*xn*CF*gsq**2*esq*eik*(ss**2+ssp**2+uu**2+uup**2)
     & /(tt*ttp)
      return
      end
