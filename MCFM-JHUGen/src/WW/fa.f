      function fa(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: fa
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      integer:: i1,i2,i3,i4
      complex(dp):: Lsm1_2mht,z2,BigT,I3m,flipbit

      real(dp):: t134,t234,del3,del12
  
C---statement function  
      z2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

c     Taken from Eq.(2.14)
      t134=s(j1,j3)+s(j1,j4)+s(j3,j4)
      t234=s(j2,j3)+s(j2,j4)+s(j3,j4)
      del12=s(j1,j2)-s(j3,j4)-s(j5,j6)
      del3=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2._dp*(s(j1,j2)*s(j3,j4)+s(j1,j2)*s(j5,j6)+s(j3,j4)*s(j5,j6))

      fa=(za(j1,j3)**2*zb(j2,j5)**2
     & /(za(j3,j4)*zb(j5,j6)*t134*(z2(j1,j5,j6,j2)))
     & -(z2(j2,j5,j6,j4)*z2(j6,j2,j5,j1))**2
     & /(zb(j3,j4)*za(j5,j6)*t134*(z2(j2,j5,j6,j1)**3)))
     & *Lsm1_2mht(s(j1,j2),t134,s(j3,j4),s(j5,j6))

      fa=fa+flipbit(j1,j2,j3,j4,j5,j6,za,zb)
     &     -flipbit(j2,j1,j5,j6,j3,j4,zb,za)

      fa=fa
     & +BigT(j1,j2,j3,j4,j5,j6,za,zb)*I3m(s(j1,j2),s(j3,j4),s(j5,j6))
     
      fa=fa+0.5_dp*(t234*del12+2._dp*s(j3,j4)*s(j5,j6))
     & /(z2(j2,j5,j6,j1)*del3)
     & *(zb(j4,j5)**2/(zb(j3,j4)*zb(j5,j6))
     &  +za(j3,j6)**2/(za(j3,j4)*za(j5,j6)))
      fa=fa+za(j3,j6)*zb(j4,j5)*(t134-t234)/(z2(j2,j5,j6,j1)*del3)
     & -0.5_dp*z2(j6,j2,j5,j4)**2
     & /(zb(j3,j4)*za(j5,j6)*t134*z2(j2,j5,j6,j1))
      return
      end



      function flipbit(j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: flipbit
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      integer:: i1,i2,i3,i4
      complex(dp):: L0,L1,Lnrat,z2,L34_12,fb

      real(dp):: t134
  
      z2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      t134=s(j1,j3)+s(j1,j4)+s(j3,j4)

      fb=0.5_dp*(za(j6,j1)*zb(j1,j4))**2*t134
     ./(zb(j3,j4)*za(j5,j6)*z2(j2,j5,j6,j1)*t134**2)*L1(-s(j3,j4),-t134)

      fb=fb+2._dp*za(j6,j1)*zb(j1,j4)*z2(j6,j2,j5,j4)*L0(-t134,-s(j3,j4))
     ./(zb(j3,j4)*za(j5,j6)*z2(j2,j5,j6,j1)*s(j3,j4))

      fb=fb-za(j1,j6)*za(j2,j6)*zb(j1,j4)**2*t134*L0(-t134,-s(j3,j4))
     ./(zb(j3,j4)*za(j5,j6)*z2(j2,j5,j6,j1)**2*s(j3,j4))

      fb=fb-0.5_dp*za(j2,j6)*zb(j1,j4)*z2(j6,j2,j5,j4)
     .*(Lnrat(-t134,-s(j3,j4))+Lnrat(-s(j1,j2),-s(j3,j4)))
     ./(zb(j3,j4)*za(j5,j6)*z2(j2,j5,j6,j1)**2)
      
      fb=fb-0.75_dp*z2(j6,j2,j5,j4)**2
     ./(zb(j3,j4)*za(j5,j6)*t134*z2(j2,j5,j6,j1))
     .*(Lnrat(-t134,-s(j3,j4))+Lnrat(-s(j1,j2),-s(j3,j4)))

      fb=fb+L34_12(j1,j2,j3,j4,j5,j6,za,zb)*lnrat(-s(j3,j4),-s(j1,j2))
      flipbit=fb
      return
      end

