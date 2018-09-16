      subroutine susana(p1,p2,p3,p4,p5,p6,za,zb,a,b)
      implicit none
      include 'types.f'
      
C     Author: R.K.Ellis
C     December 2003
C     This routine calculates the abelian (a)
C     and the non-abelian (b) base processes
C     for W-pair production. 
C     The indices represent the various polarizations
C     for the W's. 
C     Thus -1 is like left-handed
C     Thus +1 is like right-handed
C     Thus  0 is like longitudinal
C     These functions allow us to calculate W-pair production 
C     without decay correlations. 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      real(dp):: s12,s34,s56,t134
      complex(dp):: a(-1:+1,-1:+1),b(-1:+1,-1:+1)
      complex(dp):: A6treea,a6treeb
      integer:: p1,p2,p3,p4,p5,p6
      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      t134=s(p1,p3)+s(p1,p4)+s(p3,p4)

      a(-1,-1)=a6treea(p1,p2,p3,p4,p5,p6,za,zb)
      b(-1,-1)=a6treeb(p1,p2,p3,p4,p5,p6,za,zb)
      a(-1,+1)=a6treea(p1,p2,p3,p4,p6,p5,za,zb)
      b(-1,+1)=a6treeb(p1,p2,p3,p4,p6,p5,za,zb)
      a(+1,-1)=a6treea(p1,p2,p4,p3,p5,p6,za,zb)
      b(+1,-1)=a6treeb(p1,p2,p4,p3,p5,p6,za,zb)
      a(+1,+1)=a6treea(p1,p2,p4,p3,p6,p5,za,zb)
      b(+1,+1)=a6treeb(p1,p2,p4,p3,p6,p5,za,zb)

      a(-1,0) =za(p1,p3)*( 
     &    +zb(p2,p4)*(za(p2,p5)*zb(p2,p5)-za(p2,p6)*zb(p2,p6)) 
     &    +za(p5,p6)*(zb(p2,p5)*zb(p4,p6)+zb(p2,p6)*zb(p4,p5)))
     &   /(rt2*t134*s34*s56)

      a(+1,0) =za(p1,p4)*( 
     &    +zb(p2,p3)*(za(p2,p5)*zb(p2,p5)-za(p2,p6)*zb(p2,p6)) 
     &    +za(p5,p6)*(zb(p2,p5)*zb(p3,p6)+zb(p2,p6)*zb(p3,p5)))
     &   /(rt2*t134*s34*s56)

      a(0,-1)=(za(p1,p3)*(za(p6,p2)*zb(p2,p3)+za(p6,p5)*zb(p5,p3))
     &        -za(p1,p4)*(za(p6,p2)*zb(p2,p4)+za(p6,p5)*zb(p5,p4))
     & )*zb(p2,p5)/(rt2*t134*s34*s56)

      a(0,+1)=(za(p1,p3)*(za(p5,p2)*zb(p2,p3)+za(p5,p6)*zb(p6,p3))
     &        -za(p1,p4)*(za(p5,p2)*zb(p2,p4)+za(p5,p6)*zb(p6,p4))
     & )*zb(p2,p6)/(rt2*t134*s34*s56)


      a(0,0)=((za(p1,p3)*zb(p2,p3)-za(p1,p4)*zb(p2,p4))
     & *(za(p2,p5)*zb(p2,p5)-za(p2,p6)*zb(p2,p6)+za(p5,p6)*zb(p5,p6))
     & +2._dp*za(p5,p6)*zb(p2,p6)
     & *(za(p1,p3)*zb(p3,p5)-za(p1,p4)*zb(p4,p5))
     & )/(2._dp*t134*s34*s56)

      b(-1,0)=
     &  +(za(p1,p3)*zb(p1,p2)*(za(p1,p6)*zb(p4,p6)-za(p1,p5)*zb(p4,p5))
     &  + za(p1,p3)*zb(p3,p4)*(za(p3,p6)*zb(p2,p6)-za(p3,p5)*zb(p2,p5))
     &  + za(p1,p2)*zb(p2,p4)*(za(p3,p5)*zb(p2,p5)-za(p3,p6)*zb(p2,p6))
     &  + za(p3,p4)*zb(p2,p4)*(za(p1,p5)*zb(p4,p5)-za(p1,p6)*zb(p4,p6))
     &  + 2._dp*za(p1,p5)*za(p3,p6)*zb(p2,p5)*zb(p4,p6)
     &  - 2._dp*za(p1,p6)*za(p3,p5)*zb(p2,p6)*zb(p4,p5)
     & )/(2._dp*rt2*s12*s34*s56)


      b(+1,0)=
     &  +(za(p1,p4)*zb(p1,p2)*(za(p1,p6)*zb(p3,p6)-za(p1,p5)*zb(p3,p5))
     &  + za(p1,p4)*zb(p4,p3)*(za(p4,p6)*zb(p2,p6)-za(p4,p5)*zb(p2,p5))
     &  + za(p1,p2)*zb(p2,p3)*(za(p4,p5)*zb(p2,p5)-za(p4,p6)*zb(p2,p6))
     &  + za(p4,p3)*zb(p2,p3)*(za(p1,p5)*zb(p3,p5)-za(p1,p6)*zb(p3,p6))
     &  + 2._dp*za(p1,p5)*za(p4,p6)*zb(p2,p5)*zb(p3,p6)
     &  - 2._dp*za(p1,p6)*za(p4,p5)*zb(p2,p6)*zb(p3,p5)
     & )/(2._dp*rt2*s12*s34*s56)


      b(0,-1)=(
     &  +(za(p1,p2)*zb(p2,p5)+za(p1,p6)*zb(p5,p6))
     &  *(za(p4,p6)*zb(p2,p4)-za(p3,p6)*zb(p2,p3))
     &  +(za(p1,p3)*zb(p3,p5)-za(p1,p4)*zb(p4,p5))
     &  *(za(p1,p6)*zb(p1,p2)+za(p5,p6)*zb(p2,p5))
     &  -2._dp*za(p1,p4)*za(p3,p6)*zb(p2,p5)*zb(p3,p4)
     &  -2._dp*za(p1,p6)*za(p3,p4)*zb(p2,p3)*zb(p4,p5)
     & )/(2._dp*rt2*s12*s34*s56)


      b(0,1)=(
     &  +(za(p1,p2)*zb(p2,p6)+za(p1,p5)*zb(p6,p5))
     &  *(za(p4,p5)*zb(p2,p4)-za(p3,p5)*zb(p2,p3))
     &  +(za(p1,p3)*zb(p3,p6)-za(p1,p4)*zb(p4,p6))
     &  *(za(p1,p5)*zb(p1,p2)+za(p6,p5)*zb(p2,p6))
     &  -2._dp*za(p1,p4)*za(p3,p5)*zb(p2,p6)*zb(p3,p4)
     &  -2._dp*za(p1,p5)*za(p3,p4)*zb(p2,p3)*zb(p4,p6)
     & )/(2._dp*rt2*s12*s34*s56)


      b(0,0)=
     & ((za(p1,p5)*(za(p1,p4)*zb(p4,p5)-za(p1,p3)*zb(p3,p5))
     &  +za(p1,p6)*(za(p1,p3)*zb(p3,p6)-za(p1,p4)*zb(p4,p6)))*zb(p1,p2)
     & +za(p1,p2)*((za(p3,p5)*zb(p2,p5)-za(p3,p6)*zb(p2,p6))*zb(p2,p3)
     &            +(za(p4,p6)*zb(p2,p6)-za(p4,p5)*zb(p2,p5))*zb(p2,p4))
     & )/(4._dp*s12*s34*s56)

      b(0,0)=b(0,0)
     & +(za(p1,p4)*zb(p2,p4)*(za(p3,p5)*zb(p3,p5)-za(p3,p6)*zb(p3,p6))
     &  +za(p1,p5)*zb(p2,p5)*(za(p3,p6)*zb(p3,p6)-za(p4,p6)*zb(p4,p6))
     &  +za(p1,p3)*zb(p2,p3)*(za(p4,p6)*zb(p4,p6)-za(p4,p5)*zb(p4,p5))
     &  +za(p1,p6)*zb(p2,p6)*(za(p4,p5)*zb(p4,p5)-za(p3,p5)*zb(p3,p5))
     &  )/(2._dp*s12*s34*s56)

      return
      end
