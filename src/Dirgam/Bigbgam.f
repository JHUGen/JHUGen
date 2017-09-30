      double precision function Bigbgam(i1,i2,i3,i4,i5,i)
      implicit none
CC    Matrix element for q(-p1)+q(-p2)->q(p3)+q(p4)+gamma(p5)
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      integer i1,i2,i3,i4,i5,i
      double precision eik,ss,sp,tt,tp,uu,up,e12,e34,e13,e14,e23,e24
      e12=s(i1,i2)/(s(i1,i5)*s(i2,i5))
      e34=s(i3,i4)/(s(i3,i5)*s(i4,i5))
      e14=s(i1,i4)/(s(i1,i5)*s(i4,i5))
      e24=s(i2,i4)/(s(i2,i5)*s(i4,i5))
      e13=s(i1,i3)/(s(i1,i5)*s(i3,i5))
      e23=s(i2,i3)/(s(i2,i5)*s(i3,i5))
      eik=-e12+e23+e14-e34+e13+e24
      ss=s(i1,i2)
      sp=s(i3,i4)
      tt=s(i1,i3)
      tp=s(i2,i4)
      uu=s(i1,i4)
      up=s(i2,i3)
      Bigbgam=8d0*xn*CF*gsq**2*esq*Q(i)**2*eik
      Bigbgam=Bigbgam*
     . ((ss**2+sp**2+uu**2+up**2)/(tt*tp)
     . +(ss**2+sp**2+tt**2+tp**2)/(uu*up)
     . -(ss**2+sp**2)*(ss*sp-tt*tp-uu*up)/(xn*tt*tp*uu*up))

      return
      end
