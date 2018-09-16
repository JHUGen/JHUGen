c--- abbreviated version of the routine in qqb_wp2jet_g.f
      subroutine addhel_wbj(i1,i2,i3,i4,i5,i6,i7,xmsq_us_ds)
      implicit none
      include 'types.f'

      integer:: i1,i2,i3,i4,i5,i6,i7
C     u(i1)+d~(i2)+b(i3)+b~(i4)+g(i5)+l(i6)+lbar(i7)
      real(dp):: xmsq_us_ds,
     .uLdR_dLdRp,uLsL_dLsLp,uLdR_dLdRm,uLsL_dLsLm

c--- basic set of +ve gluon helicity amplitudes for QQ -> QQ
      call testem_wbj(i1,i2,i3,i4,i5,i6,i7,uLdR_dLdRp,uLsL_dLsLp)

c--- generate -ve gluon helicity amplitudes this way for now
c--- under this symmetry, uL dR -> dR dL  <===>  uR uL -> dL uR
c---                and   uL dL -> dL dL  <===>  uL uL -> dL uL
      call testem_wbj(i2,i1,i4,i3,i5,i7,i6,uLdR_dLdRm,uLsL_dLsLm)

      xmsq_us_ds=(uLdR_dLdRp+uLsL_dLsLp
     &           +uLdR_dLdRm+uLsL_dLsLm)

      return
      end

      subroutine testem_wbj(j1,j2,j3,j4,j5,j6,j7,xmsqLR,xmsqLL)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
C     u(j1)+d~(j2)+b(j3)+b~(j4)+g(j5)+l(j6)+lbar(j7)

      real(dp):: ofac,xmsqLR,xmsqLL
      real(dp):: s167,s267,s134,s234,s345
      real(dp):: s67,s34
      complex(dp):: dLRa1,dLRa2,dLRna1,dLRna2,t2
      complex(dp):: dLLa1,dLLa2,dLLna1,dLLna2

      integer:: j1,j2,j3,j4,j5,j6,j7
      t2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      s67=s(j6,j7)
      s34=s(j3,j4)
      s134=s(j1,j3)+s(j1,j4)+s(j3,j4)
      s234=s(j2,j3)+s(j2,j4)+s(j3,j4)
      s167=s(j1,j6)+s(j1,j7)+s(j6,j7)
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)

      dLRna1=
     & -zb(j7,j1)
     & *(-za(j2,j3)*t2(j4,j1,j7,j6)*t2(j5,j3,j4,j2)
     & +(-za(j2,j5)*za(j2,j3))*t2(j5,j1,j7,j6)*zb(j4,j5))
     & /(za(j2,j5)*s67*s34*s167*s345)

     & +za(j2,j3)*za(j2,j3)*zb(j7,j1)*t2(j4,j1,j7,j6)
     & /(za(j2,j5)*za(j5,j3)*s67*s167*s345)
     & -za(j2,j6)*zb(j4,j1)*zb(j7,j6)*za(j6,j2)*t2(j5,j1,j4,j3)
     & /(za(j2,j5)*s67*s34*s134*s267)

     & -za(j2,j6)*(zb(j4,j1)*t2(j7,j2,j6,j3)*t2(j5,j3,j4,j2)
     & +zb(j4,j5)*zb(j5,j1)
     & *(-za(j5,j3)*zb(j7,j6)*za(j6,j2)
     &   +za(j2,j3)*t2(j7,j2,j6,j5)))
     & /(za(j2,j5)*s67*s34*s267*s345)

     & -za(j2,j6)*za(j2,j3)*zb(j4,j1)
     & *t2(j7,j2,j6,j3)/(za(j2,j5)*za(j5,j3)*s67*s267*s345)


      dLRa1=-(
     &  -za(j2,j3)*t2(j4,j2,j3,j6)
     & *t2(j7,j1,j5,j2)/(za(j2,j5)*za(j1,j5)*s67*s34*s234)

     & -za(j2,j3)*zb(j7,j1)*zb(j4,j3)*za(j3,j2)
     & *t2(j5,j1,j7,j6)/(za(j2,j5)*s67*s34*s167*s234)

     & -za(j2,j6)*t2(j7,j2,j6,j3)*t2(j4,j1,j5,j2)
     & /(za(j2,j5)*za(j1,j5)*s67*s34*s267)

     & -za(j2,j6)*zb(j4,j1)*zb(j7,j6)*za(j6,j2)*t2(j5,j1,j4,j3)
     & /(za(j2,j5)*s67*s34*s134*s267))/xn


      dLRa2=-(
     & -zb(j7,j1)*za(j2,j3)
     & *(za(j3,j4)*t2(j4,j1,j7,j6)+za(j3,j5)*t2(j5,j1,j7,j6))
     & /(za(j3,j5)*za(j5,j4)*s67*s167*s345)

     & +za(j2,j6)
     & *(za(j3,j4)*zb(j4,j1)+za(j3,j5)*zb(j5,j1))*t2(j7,j2,j6,j3)
     & /(za(j3,j5)*za(j5,j4)*s67*s267*s345))/xn

       dLRna2=-dLRna1-xn*(dLRa1+dLRa2)

c****************************************************************

      dLLna2=-za(j2,j4)*t2(j3,j2,j4,j6)
     & *t2(j7,j1,j5,j4)/(za(j4,j5)*za(j1,j5)*s67*s34*s234)
     & -za(j2,j4)*zb(j7,j1)*zb(j3,j2)*za(j2,j4)
     & *t2(j5,j1,j7,j6)/(za(j4,j5)*s67*s34*s167*s234)
     & +zb(j7,j1)
     & *(-za(j2,j4)*t2(j3,j1,j7,j6)*zb(j5,j3)*za(j3,j4)
     & +(za(j2,j4)*za(j5,j4))
     & *t2(j5,j1,j7,j6)*zb(j3,j5))
     & /(za(j4,j5)*s67*s34*s167*s345)
     & -za(j2,j6)*t2(j7,j2,j6,j4)*t2(j3,j1,j5,j4)
     & /(za(j4,j5)*za(j1,j5)*s67*s34*s267)
     & +za(j2,j6)*t2(j7,j2,j6,j4)
     & *(+zb(j3,j1)*zb(j5,j3)*za(j3,j4)
     & -zb(j3,j5)*zb(j5,j1)*za(j5,j4))
     & /(za(j4,j5)*s67*s34*s267*s345)

      dLLa1=-(-za(j2,j4)*t2(j3,j2,j4,j6)
     & *t2(j7,j1,j5,j2)/(za(j2,j5)*za(j1,j5)*s67*s34*s234)
     .-za(j2,j4)*zb(j7,j1)*zb(j3,j4)*za(j4,j2)
     & *t2(j5,j1,j7,j6)/(za(j2,j5)*s67*s34*s167*s234)
     & -za(j2,j6)*t2(j7,j2,j6,j4)*t2(j3,j1,j5,j2)
     & /(za(j2,j5)*za(j1,j5)*s67*s34*s267)
     & -za(j2,j6)*zb(j3,j1)*t2(j7,j2,j6,j2)*t2(j5,j1,j3,j4)
     & /(za(j2,j5)*s67*s34*s134*s267))/xn


      dLLa2=-(
     & +za(j2,j4)*zb(j7,j1)
     & *(za(j3,j4)*t2(j3,j1,j7,j6)+za(j5,j4)*t2(j5,j1,j7,j6))
     & /(za(j4,j5)*za(j3,j5)*s67*s167*s345)

     & -za(j2,j6)*t2(j7,j2,j6,j4)
     & *(za(j3,j4)*zb(j3,j1)+za(j5,j4)*zb(j5,j1))
     & /(za(j4,j5)*za(j3,j5)*s67*s267*s345))/xn

      dLLna1=-dLLna2-xn*(dLLa1+dLLa2)

      ofac=8._dp*gsq**3*gwsq**2*aveqq
      ofac=ofac*s67**2/((s67-wmass**2)**2+(wmass*wwidth)**2)

C---eg uL+uR
      xmsqLR=ofac*V*xn/8._dp*(abs(dLRna1)**2+abs(dLRna2)**2
     & +abs(dLRa1)**2+abs(dLRa2)**2
     & +2._dp/xn*real((dLRa1+dLRa2)*conjg(dLRna1+dLRna2)))

C--eg uL+sL
      xmsqLL=ofac*V*xn/8._dp*(abs(dLLna1)**2+abs(dLLna2)**2
     & +abs(dLLa1)**2+abs(dLLa2)**2
     & +2._dp/xn*real((dLLa1+dLLa2)*conjg(dLLna1+dLLna2)))


      return
      end


