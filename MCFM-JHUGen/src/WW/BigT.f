      function BigT(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: BigT

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6,i1,i2,i3,i4
      complex(dp):: z2
      real(dp):: t134,t234,del3,del12,del34,del56
C---statement function
      z2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

      t134=s(j1,j3)+s(j1,j4)+s(j3,j4)
      t234=s(j2,j3)+s(j2,j4)+s(j3,j4)
      del12=s(j1,j2)-s(j3,j4)-s(j5,j6)
      del34=s(j3,j4)-s(j1,j2)-s(j5,j6)
      del56=s(j5,j6)-s(j1,j2)-s(j3,j4)
      Del3=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2._dp*(s(j1,j2)*s(j3,j4)+s(j1,j2)*s(j5,j6)+s(j3,j4)*s(j5,j6))
c----first line
      BigT=1.5_dp*s(j1,j2)*del12*(t134-t234)
     & *z2(j6,j1,j2,j5)*z2(j3,j1,j2,j4)/(z2(j2,j5,j6,j1)*Del3**2)
c----second line
      BigT=BigT
     & -0.5_dp*(3._dp*s(j1,j2)+2._dp*t134)*z2(j6,j1,j2,j5)*z2(j3,j1,j2,j4)
     & /(z2(j2,j5,j6,j1)*Del3)

c----third line
      BigT=BigT+t134/(z2(j2,j5,j6,j1)**2*Del3)
     & *(zb(j1,j4)*za(j2,j6)
     & *(za(j3,j6)*zb(j6,j5)*del56-za(j3,j4)*zb(j4,j5)*del34)
     &  -zb(j1,j5)*za(j2,j3)
     & *(za(j6,j5)*zb(j5,j4)*del56-za(j6,j3)*zb(j3,j4)*del34))

c----fourth line
      BigT=BigT+za(j3,j6)*zb(j4,j5)*s(j1,j2)*t134/(z2(j2,j5,j6,j1)*Del3)
     & -za(j3,j4)*zb(j5,j6)*z2(j6,j1,j2,j4)**2/(z2(j2,j5,j6,j1)*Del3)
     & +2._dp*za(j1,j6)*zb(j2,j4)
     & *(za(j6,j5)*zb(j5,j4)*del56-za(j6,j3)*zb(j3,j4)*del34)
     & /(zb(j3,j4)*za(j5,j6)*Del3)

c---"fourth" line from the bottom
      BigT=BigT+2._dp*z2(j6,j2,j5,j4)/(z2(j2,j5,j6,j1)*Del3)
     & *((za(j6,j5)*zb(j5,j2)*za(j2,j1)*zb(j1,j4)*del56
     &   -za(j6,j2)*zb(j2,j1)*za(j1,j3)*zb(j3,j4)*del34
     &   +z2(j6,j2,j5,j4)*s(j1,j2)*del12)/(zb(j3,j4)*za(j5,j6))
     &  +2._dp*z2(j3,j2,j6,j5)*s(j1,j2))

c---"third" line from the bottom
      BigT=BigT-zb(j1,j4)*za(j2,j6)*z2(j3,j2,j6,j5)/z2(j2,j5,j6,j1)**2
     & +2._dp*zb(j1,j5)*za(j2,j3)*z2(j6,j2,j5,j4)/z2(j2,j5,j6,j1)**2
     &     -zb(j1,j4)*za(j2,j6)*z2(j6,j2,j5,j4)*del12
     & /(zb(j3,j4)*za(j5,j6)*z2(j2,j5,j6,j1)**2)

c---penultimate two lines
      BigT=BigT+0.5_dp*(
     .+3._dp*za(j6,j2)*zb(j2,j4)*za(j6,j1)*zb(j1,j4)/(zb(j3,j4)*za(j5,j6))
     &    +za(j3,j2)*zb(j2,j5)*za(j3,j1)*zb(j1,j5)/(za(j3,j4)*zb(j5,j6))
     &    +zb(j1,j4)*za(j1,j6)*zb(j4,j5)/zb(j3,j4)
     &    -zb(j2,j4)*za(j2,j6)*za(j3,j6)/za(j5,j6)
     &    +za(j2,j3)*zb(j2,j5)*za(j3,j6)/za(j3,j4)
     &    -za(j1,j3)*zb(j1,j5)*zb(j4,j5)/zb(j5,j6)
     .+4._dp*za(j3,j6)*zb(j4,j5))/z2(j2,j5,j6,j1)

c---last line
      BigT=BigT+0.5_dp*(
     &  (za(j1,j6)*zb(j2,j4))**2/(zb(j3,j4)*za(j5,j6))
     & -(za(j1,j3)*zb(j2,j5))**2/(za(j3,j4)*zb(j5,j6)))/z2(j1,j5,j6,j2)
      BigT=BigT-0.5_dp*(zb(j1,j4)*za(j2,j6))**2
     & *(t134*del12+2._dp*s(j3,j4)*s(j5,j6))
     & /(zb(j3,j4)*za(j5,j6)*z2(j2,j5,j6,j1)**3)


      return
      end

