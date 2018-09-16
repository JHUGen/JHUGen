      subroutine dkWWamps(j1,j2,p3,p4,p5,p6,p7,Fa,Fb,Fc,Fd)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: j,j1,j2,p1,p2,p3,p4,p5,p6,p7
      complex(dp):: zab,iza,izb,Fa(2,2),Fb(2,2),Fc(2,2),Fd(2,2)
      real(dp):: s12,s34,s134,s234,s567
      zab(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      iza(p1,p2)=cone/za(p1,p2)
      izb(p1,p2)=cone/zb(p1,p2)
      do j=1,2
      if (j==1) then
      p1=j1
      p2=j2
      else
      p1=j2
      p2=j1
      endif
      s134=s(p1,p3)+s(p1,p4)+s(p3,p4)
      s234=s(p2,p3)+s(p2,p4)+s(p3,p4)
      s567=s(p5,p6)+s(p5,p7)+s(p6,p7)
      s34=s(p3,p4)
      s12=s(p1,p2)
c--- First label is the quark helicity
c--- Second label is the gluon helicity
      Fa(j,1)= + izb(p5,p7)*s567**(-1)*s34**(-1) * ( za(p2,p3)*za(p2,p7
     &    )*zb(p1,p6)*zb(p2,p4)*s234**(-1) + za(p2,p3)*za(p3,p7)*zb(p1,
     &    p6)*zb(p3,p4)*s234**(-1) )
      Fa(j,1) = Fa(j,1) + izb(p5,p7)*izb(p6,p7)*s567**(-1)*s34**(-1)
     &  * (  - za(p2,p3)*za(p2,p5)*zb(p1,p6)*zb(p2,p4)*zb(p5,p6)*
     &    s234**(-1) - za(p2,p3)*za(p3,p5)*zb(p1,p6)*zb(p3,p4)*zb(p5,p6
     &    )*s234**(-1) )

      Fa(j,2)= + iza(p5,p7)*iza(p6,p7)*s567**(-1)*s34**(-1) * ( za(p2,
     &    p3)*za(p2,p5)*za(p5,p6)*zb(p1,p6)*zb(p2,p4)*s234**(-1) + za(
     &    p2,p3)*za(p3,p5)*za(p5,p6)*zb(p1,p6)*zb(p3,p4)*s234**(-1) )
      Fa(j,2) = Fa(j,2) + iza(p6,p7)*s567**(-1)*s34**(-1) * ( za(p2,p3)
     &    *za(p2,p5)*zb(p1,p7)*zb(p2,p4)*s234**(-1) + za(p2,p3)*za(p3,
     &    p5)*zb(p1,p7)*zb(p3,p4)*s234**(-1) )

      Fb(j,1)= + izb(p5,p7)*s567**(-1)*s34**(-1) * (  - za(p1,p3)*za(p2
     &    ,p7)*zb(p1,p4)*zb(p1,p6)*s134**(-1) + za(p2,p7)*za(p3,p4)*zb(
     &    p1,p4)*zb(p4,p6)*s134**(-1) )
      Fb(j,1) = Fb(j,1) + izb(p5,p7)*izb(p6,p7)*s567**(-1)*s34**(-1)
     &  * ( za(p1,p3)*za(p2,p5)*zb(p1,p4)*zb(p1,p6)*zb(p5,p6)*
     &    s134**(-1) - za(p2,p5)*za(p3,p4)*zb(p1,p4)*zb(p4,p6)*zb(p5,p6
     &    )*s134**(-1) )

      Fb(j,2)= + iza(p5,p7)*iza(p6,p7)*s567**(-1)*s34**(-1) * (  - za(
     &    p1,p3)*za(p2,p5)*za(p5,p6)*zb(p1,p4)*zb(p1,p6)*s134**(-1) + 
     &    za(p2,p5)*za(p3,p4)*za(p5,p6)*zb(p1,p4)*zb(p4,p6)*s134**(-1)
     &     )
      Fb(j,2) = Fb(j,2) + iza(p6,p7)*s567**(-1)*s34**(-1) * (  - za(p1,
     &    p3)*za(p2,p5)*zb(p1,p4)*zb(p1,p7)*s134**(-1) + za(p2,p5)*za(
     &    p3,p4)*zb(p1,p4)*zb(p4,p7)*s134**(-1) )

      Fc(j,1)= + izb(p5,p7)*s12**(-1)*s567**(-1)*s34**(-1) * (  - za(p2
     &    ,p3)*zb(p1,p4)*zab(p7,p1,p2,p6) + za(p2,p7)*zb(p1,p6)*zab(p3,
     &    p1,p2,p4) - za(p3,p7)*zb(p4,p6)*zab(p2,p3,p4,p1) )
      Fc(j,1) = Fc(j,1) + izb(p5,p7)*izb(p6,p7)*s12**(-1)*s567**(-1)*
     & s34**(-1) * ( za(p2,p3)*zb(p1,p4)*zb(p5,p6)*zab(p5,p1,p2,p6) - 
     &    za(p2,p5)*zb(p1,p6)*zb(p5,p6)*zab(p3,p1,p2,p4) + za(p3,p5)*
     &    zb(p4,p6)*zb(p5,p6)*zab(p2,p3,p4,p1) )

      Fc(j,2)= + iza(p5,p7)*iza(p6,p7)*s12**(-1)*s567**(-1)*s34**(-1)
     &  * (  - za(p2,p3)*za(p5,p6)*zb(p1,p4)*zab(p5,p1,p2,p6) + za(p2,
     &    p5)*za(p5,p6)*zb(p1,p6)*zab(p3,p1,p2,p4) - za(p3,p5)*za(p5,p6
     &    )*zb(p4,p6)*zab(p2,p3,p4,p1) )
      Fc(j,2) = Fc(j,2) + iza(p6,p7)*s12**(-1)*s567**(-1)*s34**(-1)
     &  * (  - za(p2,p3)*zb(p1,p4)*zab(p5,p1,p2,p7) + za(p2,p5)*zb(p1,
     &    p7)*zab(p3,p1,p2,p4) - za(p3,p5)*zb(p4,p7)*zab(p2,p3,p4,p1) )

      enddo
      Fd(:,:)=-Fc(:,:)
      return
      end
