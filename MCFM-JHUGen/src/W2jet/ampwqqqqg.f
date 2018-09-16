      subroutine ampwqqqqg(j1,j2,j3,j4,j5,j6,j7,
     & xmsqLR,xmsq1LR,xmsqRL,xmsqLL,xmsqiLL,xmsqiiLL)
      implicit none
      include 'types.f'
c--- amplitudes for the production of a W, four quarks and a gluon

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'lc.f'
      include 'first.f'

      real(dp):: ofac,xmsqLR,xmsq1LR,xmsqLL,xmsqiLL,xmsqiiLL,
     &  xmsqRL
      real(dp):: s167,s267,s123,s124,s134,s145,
     & s234,s345,s235,s467,s367
      real(dp):: s67,s23,s34,s14
c      complex(dp):: dLR1,dLR2,dLR3,dLR4,dLR5,dLR6,
c     &               dLR7,dLR8,dLR9,dLR10,dLR11,dLR12
      complex(dp):: dLRa1,dLRa2,dLRna1,dLRna2,t2
c      complex(dp):: eLR1,eLR2,eLR3,eLR4,eLR5,eLR6,
c     &               eLR7,eLR8,eLR9,eLR10,eLR11,eLR12
      complex(dp):: eLRa1,eLRa2,eLRna1,eLRna2
c      complex(dp):: fRL1,fRL2,fRL3,fRL4,fRL5,fRL6,
c     &               fRL7,fRL8,fRL9,fRL10,fRL11,fRL12
      complex(dp):: fRLa1,fRLa2,fRLna1,fRLna2

c      complex(dp):: dLL1,dLL2,dLL3,dLL4,dLL5,dLL6,
c     &               dLL7,dLL8,dLL9,dLL10,dLL11,dLL12
      complex(dp):: dLLa1,dLLa2,dLLna1,dLLna2
c      complex(dp):: eLL1,eLL2,eLL3,eLL4,eLL5,eLL6,
c     &               eLL7,eLL8,eLL9,eLL10,eLL11,eLL12
      complex(dp):: eLLa1,eLLa2,eLLna1,eLLna2
      complex(dp):: diLLa1,diLLa2,diLLna1,diLLna2

c      complex(dp):: fLL1,fLL2,fLL3,fLL4,fLL5,fLL6,
c     &               fLL7,fLL8,fLL9,fLL10,fLL11,fLL12
      complex(dp):: fLLa1,fLLa2,fLLna1,fLLna2
      complex(dp):: fiLLa1,fiLLa2,fiLLna1,fiLLna2

      integer:: j1,j2,j3,j4,j5,j6,j7

      t2(j1,j2,j3,j4)=zb(j1,j2)*za(j2,j4)+zb(j1,j3)*za(j3,j4)
      s67=s(j6,j7)
      s14=s(j1,j4)
      s23=s(j2,j3)
      s34=s(j3,j4)
      s123=s(j1,j2)+s(j1,j3)+s(j2,j3)
      s124=s(j1,j2)+s(j1,j4)+s(j2,j4)
      s134=s(j1,j3)+s(j1,j4)+s(j3,j4)
      s145=s(j1,j4)+s(j1,j5)+s(j4,j5)
      s234=s(j2,j3)+s(j2,j4)+s(j3,j4)
      s235=s(j2,j3)+s(j2,j5)+s(j3,j5)
      s167=s(j1,j6)+s(j1,j7)+s(j6,j7)
      s267=s(j2,j6)+s(j2,j7)+s(j6,j7)
      s367=s(j3,j6)+s(j3,j7)+s(j6,j7)
      s467=s(j4,j6)+s(j4,j7)+s(j6,j7)
      s345=s(j3,j4)+s(j3,j5)+s(j4,j5)

      if (first) then
c         write(*,*) 'Using new bit.f'
         first = .false.
      endif

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
C----obtained from above by exhange (2<-->4)

      eLRna1=-za(j4,j3)*t2(j2,j4,j3,j6)
     & *t2(j7,j1,j5,j2)/(za(j2,j5)*za(j1,j5)*s67*s23*s234)

     & -za(j4,j3)*zb(j7,j1)*t2(j2,j4,j3,j2)
     & *t2(j5,j1,j7,j6)/(za(j2,j5)*s67*s23*s167*s234)

     & +zb(j7,j1)
     & *(-za(j4,j3)*t2(j2,j1,j7,j6)*zb(j5,j3)*za(j3,j2)
     & +(za(j4,j2)*za(j5,j3)-za(j4,j5)*za(j2,j3))
     & *t2(j5,j1,j7,j6)*zb(j2,j5))
     & /(za(j2,j5)*s67*s23*s167*s235)

     & -zb(j7,j1)*za(j4,j3)
     & *t2(j5,j1,j7,j6)
     & /(za(j5,j2)*s67*s167*s235)

     & -za(j4,j6)*t2(j7,j4,j6,j3)*t2(j2,j1,j5,j2)
     & /(za(j2,j5)*za(j1,j5)*s67*s23*s467)

     & +za(j4,j6)
     & *(+zb(j2,j1)*t2(j7,j4,j6,j3)*zb(j5,j3)*za(j3,j2)
     & +zb(j2,j5)*zb(j5,j1)
     & *(-za(j5,j3)*t2(j7,j4,j6,j2)
     &   +za(j2,j3)*t2(j7,j4,j6,j5)))
     & /(za(j2,j5)*s67*s23*s467*s235)

     & +za(j4,j6)
     & *(zb(j5,j1))*t2(j7,j4,j6,j3)/(za(j5,j2)*s67*s467*s235)



      eLRa1=-(-zb(j7,j1)*za(j4,j3)
     & *(za(j3,j2)*t2(j2,j1,j7,j6)+za(j3,j5)*t2(j5,j1,j7,j6))
     & /(za(j3,j5)*za(j5,j2)*s67*s167*s235)

     & +za(j4,j6)
     & *(za(j3,j2)*zb(j2,j1)+za(j3,j5)*zb(j5,j1))*t2(j7,j4,j6,j3)
     & /(za(j3,j5)*za(j5,j2)*s67*s467*s235))/xn


      eLRa2=-(-za(j4,j3)*t2(j2,j4,j3,j6)
     & *t2(j7,j1,j5,j4)/(za(j4,j5)*za(j1,j5)*s67*s23*s234)

     &  -za(j4,j3)*zb(j7,j1)*t2(j2,j4,j3,j4)
     & *t2(j5,j1,j7,j6)/(za(j4,j5)*s67*s23*s167*s234)

     & -za(j4,j6)*t2(j7,j4,j6,j3)*t2(j2,j1,j5,j4)
     & /(za(j4,j5)*za(j1,j5)*s67*s23*s467)

     & -za(j4,j6)*zb(j2,j1)*t2(j7,j4,j6,j4)*t2(j5,j1,j2,j3)
     & /(za(j4,j5)*s67*s23*s123*s467))/xn


      eLRna2=-eLRna1-xn*(eLRa1+eLRa2)

c****************************************************************
C obtained from dLR by exchanging  (1<-->3)

      fRLa1=-(-za(j2,j1)*t2(j4,j2,j1,j6)
     & *t2(j7,j3,j5,j2)/(za(j2,j5)*za(j3,j5)*s67*s14*s124)

     & -za(j2,j1)*zb(j7,j3)*zb(j4,j1)*za(j1,j2)
     & *t2(j5,j3,j7,j6)/(za(j2,j5)*s67*s14*s367*s124)


     & -za(j2,j6)*t2(j7,j2,j6,j1)*t2(j4,j3,j5,j2)
     & /(za(j2,j5)*za(j3,j5)*s67*s14*s267)
     & -za(j2,j6)*zb(j4,j3)*zb(j7,j6)*za(j6,j2)*t2(j5,j3,j4,j1)
     & /(za(j2,j5)*s67*s14*s134*s267))/xn


      fRLa2=-(-zb(j7,j3)*za(j2,j1)
     & *(za(j1,j4)*t2(j4,j3,j7,j6)+za(j1,j5)*t2(j5,j3,j7,j6))
     & /(za(j1,j5)*za(j5,j4)*s67*s367*s145)
     &  +za(j2,j6)
     & *(za(j1,j4)*zb(j4,j3)+za(j1,j5)*zb(j5,j3))*t2(j7,j2,j6,j1)
     & /(za(j1,j5)*za(j5,j4)*s67*s267*s145))/xn


      fRLna2=-za(j2,j1)*t2(j4,j2,j1,j6)
     & *t2(j7,j3,j5,j1)/(za(j1,j5)*za(j3,j5)*s67*s14*s124)

     & -za(j2,j1)*zb(j7,j3)*zb(j4,j2)*za(j2,j1)
     & *t2(j5,j3,j7,j6)/(za(j1,j5)*s67*s14*s367*s124)

     & +zb(j7,j3)
     & *(-za(j2,j1)*t2(j4,j3,j7,j6)*zb(j5,j4)*za(j4,j1)
     & +(za(j2,j1)*za(j5,j1))
     & *t2(j5,j3,j7,j6)*zb(j4,j5))
     & /(za(j1,j5)*s67*s14*s367*s145)

     & -zb(j7,j3)*za(j2,j1)
     & *(za(j1,j4)*t2(j4,j3,j7,j6)+za(j1,j5)*t2(j5,j3,j7,j6))
     & /(za(j1,j5)*za(j5,j4)*s67*s367*s145)

     & -za(j2,j6)*t2(j7,j2,j6,j1)*t2(j4,j3,j5,j1)
     & /(za(j1,j5)*za(j3,j5)*s67*s14*s267)

     & +za(j2,j6)*t2(j7,j2,j6,j1)
     & *(+zb(j4,j3)*zb(j5,j4)*za(j4,j1)
     & -zb(j4,j5)*zb(j5,j3)*za(j5,j1))
     & /(za(j1,j5)*s67*s14*s267*s145)

     & +za(j2,j6)
     & *(za(j1,j4)*zb(j4,j3)+za(j1,j5)*zb(j5,j3))*t2(j7,j2,j6,j1)
     & /(za(j1,j5)*za(j5,j4)*s67*s267*s145)


      fRLna1=-fRLna2-xn*(fRLa1+fRLa2)

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

c****************************************************************
C  ell1 etc obtained from dLL1 etc by exchanging 2<-->4 and changing sign

c****************************************************************
      eLLa1=-(
     & -za(j4,j2)*zb(j7,j1)
     & *(za(j3,j2)*t2(j3,j1,j7,j6)+za(j5,j2)*t2(j5,j1,j7,j6))
     & /(za(j2,j5)*za(j3,j5)*s67*s167*s235)

     & +za(j4,j6)*t2(j7,j4,j6,j2)
     & *(za(j3,j2)*zb(j3,j1)+za(j5,j2)*zb(j5,j1))
     & /(za(j2,j5)*za(j3,j5)*s67*s467*s235))/xn


      eLLa2=-(+za(j4,j2)*t2(j3,j4,j2,j6)
     & *t2(j7,j1,j5,j4)/(za(j4,j5)*za(j1,j5)*s67*s23*s234)

     & +za(j4,j2)*zb(j7,j1)*zb(j3,j2)*za(j2,j4)
     & *t2(j5,j1,j7,j6)/(za(j4,j5)*s67*s23*s167*s234)

     & +za(j4,j6)*t2(j7,j4,j6,j2)*t2(j3,j1,j5,j4)
     & /(za(j4,j5)*za(j1,j5)*s67*s23*s467)

     & +za(j4,j6)*zb(j3,j1)*zb(j7,j6)*za(j6,j4)*t2(j5,j1,j3,j2)
     & /(za(j4,j5)*s67*s23*s123*s467)

     & )/xn



      eLLna1=+za(j4,j2)*t2(j3,j4,j2,j6)
     & *t2(j7,j1,j5,j2)/(za(j2,j5)*za(j1,j5)*s67*s23*s234)

     & +za(j4,j2)*zb(j7,j1)*zb(j3,j4)*za(j4,j2)
     & *t2(j5,j1,j7,j6)/(za(j2,j5)*s67*s23*s167*s234)

     & -zb(j7,j1)
     & *(-za(j4,j2)*t2(j3,j1,j7,j6)*zb(j5,j3)*za(j3,j2)
     & +(za(j4,j2)*za(j5,j2)-za(j4,j5)*za(j2,j2))
     & *t2(j5,j1,j7,j6)*zb(j3,j5))
     & /(za(j2,j5)*s67*s23*s167*s235)

     & +za(j4,j6)*t2(j7,j4,j6,j2)*t2(j3,j1,j5,j2)
     & /(za(j2,j5)*za(j1,j5)*s67*s23*s467)

     & -za(j4,j6)
     & *(+zb(j3,j1)*t2(j7,j4,j6,j2)*zb(j5,j3)*za(j3,j2)
     & +zb(j3,j5)*zb(j5,j1)
     & *(-za(j5,j2)*t2(j7,j4,j6,j2)
     &   +za(j2,j2)*t2(j7,j4,j6,j5)))
     & /(za(j2,j5)*s67*s23*s467*s235)


      eLLna2=-eLLna1-xn*(eLLa1+eLLa2)

c****************************************************************
C  fLL1 etc obtained from dLL1 etc by exchanging (1<-->3) and changing sign

c****************************************************************

      fLLa1=-(+za(j2,j4)*t2(j1,j2,j4,j6)
     & *t2(j7,j3,j5,j2)/(za(j2,j5)*za(j3,j5)*s67*s14*s124)
     & +za(j2,j4)*zb(j7,j3)*zb(j1,j4)*za(j4,j2)
     & *t2(j5,j3,j7,j6)/(za(j2,j5)*s67*s14*s367*s124)
     & +za(j2,j6)*t2(j7,j2,j6,j4)*t2(j1,j3,j5,j2)
     & /(za(j2,j5)*za(j3,j5)*s67*s14*s267)
     & +za(j2,j6)*zb(j1,j3)*zb(j7,j6)*za(j6,j2)*t2(j5,j3,j1,j4)
     & /(za(j2,j5)*s67*s14*s134*s267))/xn


      fLLa2=-(
     & -za(j2,j4)*zb(j7,j3)
     & *(za(j1,j4)*t2(j1,j3,j7,j6)+za(j5,j4)*t2(j5,j3,j7,j6))
     & /(za(j4,j5)*za(j1,j5)*s67*s367*s145)
     & +za(j2,j6)*t2(j7,j2,j6,j4)
     & *(za(j1,j4)*zb(j1,j3)+za(j5,j4)*zb(j5,j3))
     & /(za(j4,j5)*za(j1,j5)*s67*s267*s145))/xn

      fLLna2=+za(j2,j4)*t2(j1,j2,j4,j6)
     & *t2(j7,j3,j5,j4)/(za(j4,j5)*za(j3,j5)*s67*s14*s124)

     & +za(j2,j4)*zb(j7,j3)*zb(j1,j2)*za(j2,j4)
     & *t2(j5,j3,j7,j6)/(za(j4,j5)*s67*s14*s367*s124)

     & -zb(j7,j3)
     & *(-za(j2,j4)*t2(j1,j3,j7,j6)*zb(j5,j1)*za(j1,j4)
     & +za(j2,j4)*za(j5,j4)*t2(j5,j3,j7,j6)*zb(j1,j5))
     & /(za(j4,j5)*s67*s14*s367*s145)

     & +za(j2,j6)*t2(j7,j2,j6,j4)*t2(j1,j3,j5,j4)
     & /(za(j4,j5)*za(j3,j5)*s67*s14*s267)

     & -za(j2,j6)
     & *(zb(j1,j3)*zb(j5,j1)*za(j1,j4)
     & -zb(j1,j5)*zb(j5,j3)*za(j5,j4))*t2(j7,j2,j6,j4)
     & /(za(j4,j5)*s67*s14*s267*s145)

       fLLna1=-fLLna2-xn*(fLLa1+fLLa2)

c****************************************************************

C    ans=
c    - t(cb,i2,i1)*D(i3,i4)/2/n*(dLR1+dLR2+dLR3+dLR7+dLR8+dLR9)
c    - t(cb,i4,i3)*D(i1,i2)/2/n*(dLR5+dLR6+dLR11+dLR12)
c + t(cb,i2,i3)*D(i1,i4)/2*(dLR3-dLR4+dLR6+dLR8+dLR9-dLR10+dLR12)
c + t(cb,i4,i1)*D(i2,i3)/2*(dLR1+dLR2+dLR4+dLR5+dLR7+dLR10+dLR11)

      if (colourchoice ==1) then
      diLLa1=eLLna1
      diLLa2=eLLna2
      diLLna1=dLLna1
      diLLna2=dLLna2

      fiLLa1=fLLna1
      fiLLa2=fLLna2
      fiLLna1=dLLna1
      fiLLna2=dLLna2
      else
      diLLa1=dLLa1+eLLna1
      diLLa2=dLLa2+eLLna2
      diLLna1=dLLna1+eLLa1
      diLLna2=dLLna2+eLLa2

      fiLLa1=dLLa1+fLLna1
      fiLLa2=dLLa2+fLLna2
      fiLLna1=dLLna1+fLLa1
      fiLLna2=dLLna2+fLLa2
      endif


      ofac=8._dp*gsq**3*gwsq**2*aveqq
      ofac=ofac*s67**2/((s67-wmass**2)**2+(wmass*wwidth)**2)

      if (colourchoice == 1) then

C---eg uL+uR
      xmsqLR=ofac*V*xn/8._dp*(abs(dLRna1)**2+abs(dLRna2)**2)

C---eg uL+dR (second term)
      xmsq1LR=ofac*V*xn/8._dp*(abs(eLRna1)**2+abs(eLRna2)**2)

C--eg uL+sL
      xmsqLL=ofac*V*xn/8._dp*(abs(dLLna1)**2+abs(dLLna2)**2)

C--eg uL+dL
      xmsqiLL=ofac*V*xn/8._dp*(abs(diLLa1)**2+abs(diLLa2)**2
     & +abs(diLLna1)**2+abs(diLLna2)**2)

C--eg uL+uL
      xmsqiiLL=ofac*V*xn/8._dp*(abs(fiLLna1)**2+abs(fiLLna2)**2
     & +abs(fiLLa1)**2+abs(fiLLa2)**2)

C--eg uR+uL
      xmsqRL=ofac*V*xn/8._dp*(abs(fRLna1)**2+abs(fRLna2)**2)

      else

C---eg uL+uR
      xmsqLR=ofac*V*xn/8._dp*(abs(dLRna1)**2+abs(dLRna2)**2
     & +abs(dLRa1)**2+abs(dLRa2)**2
     & +2._dp/xn*real((dLRa1+dLRa2)*conjg(dLRna1+dLRna2)))
C---eg uL+dR (second term)
      xmsq1LR=ofac*V*xn/8._dp*(abs(eLRna1)**2+abs(eLRna2)**2
     & +abs(eLRa1)**2+abs(eLRa2)**2
     & +2._dp/xn*real((eLRa1+eLRa2)*conjg(eLRna1+eLRna2)))

C--eg uL+sL
      xmsqLL=ofac*V*xn/8._dp*(abs(dLLna1)**2+abs(dLLna2)**2
     & +abs(dLLa1)**2+abs(dLLa2)**2
     & +2._dp/xn*real((dLLa1+dLLa2)*conjg(dLLna1+dLLna2)))

C--eg uL+dL
      xmsqiLL=ofac*V*xn/8._dp*(abs(diLLa1)**2+abs(diLLa2)**2
     & +abs(diLLna1)**2+abs(diLLna2)**2
     & +2._dp/xn*real((diLLa1+diLLa2)*conjg(diLLna1+diLLna2)))
C--eg uL+uL
      xmsqiiLL=ofac*V*xn/8._dp*(abs(fiLLna1)**2+abs(fiLLna2)**2
     & +abs(fiLLa1)**2+abs(fiLLa2)**2
     & +2._dp/xn*real((fiLLa1+fiLLa2)*conjg(fiLLna1+fiLLna2)))

C--eg uR+uL
      xmsqRL=ofac*V*xn/8._dp*(abs(fRLna1)**2+abs(fRLna2)**2
     & +abs(fRLa1)**2+abs(fRLa2)**2
     & +2._dp/xn*real((fRLa1+fRLa2)*conjg(fRLna1+fRLna2)))
      endif
      return
      end


