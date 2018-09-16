      subroutine mbc(st,j1,j2,j3,j4,j5,j6,za,zb,bcoeff)
      implicit none
      include 'types.f'

c----Bubble coefficients extracted from BDK 11.5, 11.8
c----after performing the transformation
c    (1-->4)
c    (2-->3)
c    (3-->1)
c    (4-->2)
c    (5-->5)
c    (6-->6)
      integer:: j1,j2,j3,j4,j5,j6,j,k
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'blabels.f'
      include 'kprocess.f'
      character*9 st
      complex(dp):: bcoeff(8),zab2,izab2,bub6
      complex(dp):: iza(6,6),izb(6,6)
      real(dp):: is(6,6),t,t134ms34,t234ms34,t134ms56,t234ms56
      real(dp):: IDelta
      zab2(j1,j2,j3,j4)=+za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)
      izab2(j1,j2,j3,j4)=cone/zab2(j1,j2,j3,j4)
      t(j1,j2,j3)=s(j1,j2)+s(j2,j3)+s(j1,j3)
      do j=1,6
      do k=j+1,6
      is(j,k)=1._dp/s(j,k)
      iza(j,k)=cone/za(j,k)
      izb(j,k)=cone/zb(j,k)
      is(k,j)=is(j,k)
      iza(k,j)=-iza(j,k)
      izb(k,j)=-izb(j,k)
      enddo
      enddo

      IDelta=s(j1,j2)**2+s(j3,j4)**2+s(j5,j6)**2
     & -2._dp*s(j1,j2)*s(j3,j4)
     & -2._dp*s(j3,j4)*s(j5,j6)
     & -2._dp*s(j1,j2)*s(j5,j6)
      IDelta=1._dp/IDelta

      t134ms34=t(j1,j3,j4)-s(j3,j4)
      t234ms34=t(j2,j3,j4)-s(j3,j4)
      t134ms56=t(j1,j3,j4)-s(j5,j6)
      t234ms56=t(j2,j3,j4)-s(j5,j6)

      do j=1,8
      bcoeff(j)=czip
      enddo

      if     (st=='q+qb-g+g+') then

      bcoeff(b34)= + t134ms34**(-2) * (  - za(j4,j3)*za(j1,j5)**2*zb(j4
     &    ,j1)**2*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b34) = bcoeff(b34) + t134ms34**(-1) * ( za(j3,j1)*za(j1,j5
     &    )*za(j2,j5)*zb(j4,j1)*iza(j1,j2)**3*iza(j5,j6) + za(j3,j2)*
     &    za(j1,j5)**2*zb(j4,j1)*iza(j1,j2)**3*iza(j5,j6) - za(j3,j5)*
     &    za(j1,j5)*zb(j4,j1)*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b34) = bcoeff(b34) + t234ms34**(-2) * (  - za(j4,j3)*za(j2
     &    ,j5)**2*zb(j4,j2)**2*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b34) = bcoeff(b34) + t234ms34**(-1) * (  - za(j3,j1)*za(j2
     &    ,j5)**2*zb(j4,j2)*iza(j1,j2)**3*iza(j5,j6) - za(j3,j2)*za(j1,
     &    j5)*za(j2,j5)*zb(j4,j2)*iza(j1,j2)**3*iza(j5,j6) - za(j3,j5)*
     &    za(j2,j5)*zb(j4,j2)*iza(j1,j2)**2*iza(j5,j6) )

      bcoeff(b56)= + t234ms56**(-2) * (  - za(j3,j1)**2*za(j5,j6)*zb(j1
     &    ,j6)**2*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b56) = bcoeff(b56) + t234ms56**(-1) * (  - za(j3,j1)**2*
     &    za(j2,j5)*zb(j1,j6)*iza(j4,j3)*iza(j1,j2)**3 - za(j3,j1)*za(
     &    j3,j2)*za(j1,j5)*zb(j1,j6)*iza(j4,j3)*iza(j1,j2)**3 - za(j3,
     &    j1)*za(j3,j5)*zb(j1,j6)*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b56) = bcoeff(b56) + t134ms56**(-2) * (  - za(j3,j2)**2*
     &    za(j5,j6)*zb(j2,j6)**2*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b56) = bcoeff(b56) + t134ms56**(-1) * ( za(j3,j1)*za(j3,j2
     &    )*za(j2,j5)*zb(j2,j6)*iza(j4,j3)*iza(j1,j2)**3 + za(j3,j2)**2
     &    *za(j1,j5)*zb(j2,j6)*iza(j4,j3)*iza(j1,j2)**3 - za(j3,j2)*za(
     &    j3,j5)*zb(j2,j6)*iza(j4,j3)*iza(j1,j2)**2 )

      bcoeff(b134)= + t134ms34**(-2) * ( za(j4,j3)*za(j1,j5)**2*zb(j4,
     &    j1)**2*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b134) = bcoeff(b134) + t134ms34**(-1) * (  - za(j3,j1)*za(
     &    j1,j5)*za(j2,j5)*zb(j4,j1)*iza(j1,j2)**3*iza(j5,j6) - za(j3,
     &    j2)*za(j1,j5)**2*zb(j4,j1)*iza(j1,j2)**3*iza(j5,j6) + za(j3,
     &    j5)*za(j1,j5)*zb(j4,j1)*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b134) = bcoeff(b134) + t134ms56**(-2) * ( za(j3,j2)**2*za(
     &    j5,j6)*zb(j2,j6)**2*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b134) = bcoeff(b134) + t134ms56**(-1) * (  - za(j3,j1)*za(
     &    j3,j2)*za(j2,j5)*zb(j2,j6)*iza(j4,j3)*iza(j1,j2)**3 - za(j3,
     &    j2)**2*za(j1,j5)*zb(j2,j6)*iza(j4,j3)*iza(j1,j2)**3 + za(j3,
     &    j2)*za(j3,j5)*zb(j2,j6)*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b134) = bcoeff(b134) - za(j3,j1)*za(j3,j5)*za(j2,j5)*iza(
     & j4,j3)*iza(j1,j2)**3*iza(j5,j6) - za(j3,j2)*za(j3,j5)*za(j1,j5)*
     &    iza(j4,j3)*iza(j1,j2)**3*iza(j5,j6)

      bcoeff(b234)= + t234ms34**(-2) * ( za(j4,j3)*za(j2,j5)**2*zb(j4,
     &    j2)**2*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b234) = bcoeff(b234) + t234ms34**(-1) * ( za(j3,j1)*za(j2,
     &    j5)**2*zb(j4,j2)*iza(j1,j2)**3*iza(j5,j6) + za(j3,j2)*za(j1,
     &    j5)*za(j2,j5)*zb(j4,j2)*iza(j1,j2)**3*iza(j5,j6) + za(j3,j5)*
     &    za(j2,j5)*zb(j4,j2)*iza(j1,j2)**2*iza(j5,j6) )
      bcoeff(b234) = bcoeff(b234) + t234ms56**(-2) * ( za(j3,j1)**2*za(
     &    j5,j6)*zb(j1,j6)**2*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b234) = bcoeff(b234) + t234ms56**(-1) * ( za(j3,j1)**2*za(
     &    j2,j5)*zb(j1,j6)*iza(j4,j3)*iza(j1,j2)**3 + za(j3,j1)*za(j3,
     &    j2)*za(j1,j5)*zb(j1,j6)*iza(j4,j3)*iza(j1,j2)**3 + za(j3,j1)*
     &    za(j3,j5)*zb(j1,j6)*iza(j4,j3)*iza(j1,j2)**2 )
      bcoeff(b234) = bcoeff(b234) + za(j3,j1)*za(j3,j5)*za(j2,j5)*iza(
     & j4,j3)*iza(j1,j2)**3*iza(j5,j6) + za(j3,j2)*za(j3,j5)*za(j1,j5)*
     &    iza(j4,j3)*iza(j1,j2)**3*iza(j5,j6)

      bcoeff(rat)= + t134ms34**(-1) * ( za(j4,j3)*za(j1,j5)**2*zb(j4,j1
     &    )**2*iza(j1,j2)**2*iza(j5,j6)*is(j4,j3) )
      bcoeff(rat) = bcoeff(rat) + t234ms34**(-1) * ( za(j4,j3)*za(j2,j5
     &    )**2*zb(j4,j2)**2*iza(j1,j2)**2*iza(j5,j6)*is(j4,j3) )
      bcoeff(rat) = bcoeff(rat) + t234ms56**(-1) * ( za(j3,j1)**2*za(j5
     &    ,j6)*zb(j1,j6)**2*iza(j4,j3)*iza(j1,j2)**2*is(j5,j6) )
      bcoeff(rat) = bcoeff(rat) + t134ms56**(-1) * ( za(j3,j2)**2*za(j5
     &    ,j6)*zb(j2,j6)**2*iza(j4,j3)*iza(j1,j2)**2*is(j5,j6) )
      bcoeff(rat) = bcoeff(rat) + za(j3,j5)**2*iza(j4,j3)*iza(j1,j2)**2
     & *iza(j5,j6) - zb(j4,j6)**2*iza(j1,j2)**2*izb(j4,j3)*izb(j5,j6)

      if     ((kcase==kHWWint) .or. (kcase==kWWqqbr)
     &   .or. (kcase==kHWWHpi) .or. (kcase==kggWW4l)
     &   .or. (kcase==kggWWbx)
     &   .or. (kcase==kggVV4l) .or. (kcase==kggVVbx)) then
      bcoeff(b12)=czip
      bcoeff(b12zm)=czip
      bcoeff(b12m)=czip
      elseif ((kcase==kHZZint) .or. (kcase==kZZlept)
     &   .or. (kcase==kHZZHpi) .or. (kcase==kggZZ4l)
     &   .or. (kcase==kggZZbx)) then
      bcoeff(b12zm)=czip
      bcoeff(b12)=bcoeff(b12zm)
      bcoeff(b1)=-bcoeff(b12)-bcoeff(b34)-bcoeff(b56)
     & -bcoeff(b134)-bcoeff(b234)
      else
      write(6,*) 'Unimplemented case in mbc.f'
      stop
      endif

      elseif (st=='q+qb-g+g-') then

      bcoeff(b34)= + t134ms34**(-2) * ( za(j3,j1)**2*zb(j4,j3)*zb(j1,j6
     &    )**2*izb(j5,j6)*izab2(j1,j4,j3,j2)**2*s(j4,j3) )
      bcoeff(b34) = bcoeff(b34) + t134ms34**(-1) * ( 2._dp*za(j3,j1)**2*
     &    zb(j4,j3)*zb(j1,j6)**2*izb(j5,j6)*izab2(j1,j4,j3,j2)**2 - 2._dp
     &    *za(j3,j1)*zb(j4,j3)*zb(j1,j6)*izb(j5,j6)*zab2(j3,j4,j1,j2)*
     &    zab2(j1,j4,j3,j6)*izab2(j1,j4,j3,j2)**3 )
      bcoeff(b34) = bcoeff(b34) + t234ms34**(-2) * ( za(j4,j3)*za(j2,j5
     &    )**2*zb(j4,j2)**2*iza(j5,j6)*izab2(j1,j3,j4,j2)**2*s(j4,j3) )
      bcoeff(b34) = bcoeff(b34) + t234ms34**(-1) * ( 2._dp*za(j4,j3)*za(
     &    j2,j5)**2*zb(j4,j2)**2*iza(j5,j6)*izab2(j1,j3,j4,j2)**2 - 2._dp
     &    *za(j4,j3)*za(j2,j5)*zb(j4,j2)*iza(j5,j6)*zab2(j1,j3,j2,j4)*
     &    zab2(j5,j3,j4,j2)*izab2(j1,j3,j4,j2)**3 )
      bcoeff(b34) = bcoeff(b34) + IDelta**2 * (  - 3._dp*zab2(j3,j1,j2,
     &    j4)*zab2(j2,j4,j3,j1)*zab2(j5,j1,j2,j6)*izab2(j1,j4,j3,j2)*s(
     &    j4,j3) - 3._dp*zab2(j3,j1,j2,j4)*zab2(j2,j4,j3,j1)*zab2(j5,j1,
     &    j2,j6)*izab2(j1,j4,j3,j2)*s(j1,j2) + 3._dp*zab2(j3,j1,j2,j4)*
     &    zab2(j2,j4,j3,j1)*zab2(j5,j1,j2,j6)*izab2(j1,j4,j3,j2)*s(j5,
     &    j6) - 3._dp*zab2(j3,j2,j1,j4)*zab2(j2,j3,j4,j1)*zab2(j5,j2,j1,
     &    j6)*izab2(j1,j3,j4,j2)*s(j4,j3) - 3._dp*zab2(j3,j2,j1,j4)*
     &    zab2(j2,j3,j4,j1)*zab2(j5,j2,j1,j6)*izab2(j1,j3,j4,j2)*s(j1,
     &    j2) + 3._dp*zab2(j3,j2,j1,j4)*zab2(j2,j3,j4,j1)*zab2(j5,j2,j1,
     &    j6)*izab2(j1,j3,j4,j2)*s(j5,j6) )
      bcoeff(b34) = bcoeff(b34) + IDelta * ( 2._dp*za(j4,j3)*za(j3,j1)*
     &    za(j5,j6)*zb(j4,j3)*zb(j4,j6)*zb(j2,j6)*izab2(j1,j4,j3,j2)**3
     &    *t(j4,j3,j1) - 2._dp*za(j4,j3)*za(j3,j1)*za(j5,j6)*zb(j4,j3)*
     &    zb(j4,j6)*zb(j2,j6)*izab2(j1,j4,j3,j2)**3*t(j4,j3,j2) + 2._dp*
     &    za(j4,j3)*za(j3,j5)*za(j1,j5)*zb(j4,j3)*zb(j4,j2)*zb(j5,j6)*
     &    izab2(j1,j3,j4,j2)**3*t(j3,j4,j1) - 2._dp*za(j4,j3)*za(j3,j5)*
     &    za(j1,j5)*zb(j4,j3)*zb(j4,j2)*zb(j5,j6)*izab2(j1,j3,j4,j2)**3
     &    *t(j3,j4,j2) - za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(j4,j6)*izab2(
     &    j1,j4,j3,j2)**2*s(j4,j3) + za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(
     &    j4,j6)*izab2(j1,j4,j3,j2)**2*s(j1,j2) + za(j4,j3)*za(j3,j5)*
     &    zb(j4,j3)*zb(j4,j6)*izab2(j1,j4,j3,j2)**2*s(j5,j6) - za(j4,j3
     &    )*za(j3,j5)*zb(j4,j3)*zb(j4,j6)*izab2(j1,j3,j4,j2)**2*s(j4,j3
     &    ) + za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(j4,j6)*izab2(j1,j3,j4,j2
     &    )**2*s(j1,j2) + za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(j4,j6)*
     &    izab2(j1,j3,j4,j2)**2*s(j5,j6) + 2._dp*za(j4,j3)*za(j1,j5)*zb(
     &    j4,j2)*zb(j4,j6)*izab2(j1,j3,j4,j2)**3*t(j3,j4,j1)*t(j3,j4,j2
     &    ) )
      bcoeff(b34) = bcoeff(b34) + IDelta * (  - 2._dp*za(j4,j3)*za(j1,j5
     &    )*zb(j4,j2)*zb(j4,j6)*izab2(j1,j3,j4,j2)**3*t(j3,j4,j2)**2 -
     &    2._dp*za(j4,j3)*za(j2,j5)*zb(j4,j2)*zb(j4,j6)*izab2(j1,j3,j4,
     &    j2)**2*t(j3,j4,j1) + 2._dp*za(j4,j3)*za(j2,j5)*zb(j4,j2)*zb(j4
     &    ,j6)*izab2(j1,j3,j4,j2)**2*t(j3,j4,j2) + 1._dp/2._dp*za(j4,j3)*
     &    zb(j4,j6)**2*izb(j5,j6)*zab2(j2,j4,j3,j1)*izab2(j1,j4,j3,j2)
     &     + 1._dp/2._dp*za(j4,j3)*zb(j4,j6)**2*izb(j5,j6)*zab2(j2,j3,j4,
     &    j1)*izab2(j1,j3,j4,j2) - za(j4,j3)*zb(j4,j6)**2*izb(j5,j6)*
     &    izab2(j1,j3,j4,j2)**2*s(j4,j3)*t(j3,j4,j2) + za(j4,j3)*zb(j4,
     &    j6)**2*izb(j5,j6)*izab2(j1,j3,j4,j2)**2*s(j1,j2)*t(j3,j4,j2)
     &     + za(j4,j3)*zb(j4,j6)**2*izb(j5,j6)*izab2(j1,j3,j4,j2)**2*s(
     &    j5,j6)*t(j3,j4,j2) - 2._dp*za(j3,j1)*za(j3,j5)*zb(j4,j3)*zb(j1
     &    ,j6)*izab2(j1,j4,j3,j2)**2*t(j4,j3,j1) + 2._dp*za(j3,j1)*za(j3
     &    ,j5)*zb(j4,j3)*zb(j1,j6)*izab2(j1,j4,j3,j2)**2*t(j4,j3,j2) +
     &    2._dp*za(j3,j1)*za(j3,j5)*zb(j4,j3)*zb(j2,j6)*izab2(j1,j4,j3,
     &    j2)**3*t(j4,j3,j1)**2 )
      bcoeff(b34) = bcoeff(b34) + IDelta * (  - 2._dp*za(j3,j1)*za(j3,j5
     &    )*zb(j4,j3)*zb(j2,j6)*izab2(j1,j4,j3,j2)**3*t(j4,j3,j1)*t(j4,
     &    j3,j2) + 1._dp/2._dp*za(j3,j5)**2*zb(j4,j3)*iza(j5,j6)*zab2(j2,
     &    j4,j3,j1)*izab2(j1,j4,j3,j2) + 1._dp/2._dp*za(j3,j5)**2*zb(j4,
     &    j3)*iza(j5,j6)*zab2(j2,j3,j4,j1)*izab2(j1,j3,j4,j2) - za(j3,
     &    j5)**2*zb(j4,j3)*iza(j5,j6)*izab2(j1,j4,j3,j2)**2*s(j4,j3)*t(
     &    j4,j3,j1) + za(j3,j5)**2*zb(j4,j3)*iza(j5,j6)*izab2(j1,j4,j3,
     &    j2)**2*s(j1,j2)*t(j4,j3,j1) + za(j3,j5)**2*zb(j4,j3)*iza(j5,
     &    j6)*izab2(j1,j4,j3,j2)**2*s(j5,j6)*t(j4,j3,j1) - za(j3,j5)*
     &    zb(j4,j6)*zab2(j2,j4,j3,j1)*izab2(j1,j4,j3,j2) - za(j3,j5)*
     &    zb(j4,j6)*zab2(j2,j3,j4,j1)*izab2(j1,j3,j4,j2) )

      bcoeff(b12zm)= + IDelta**2 * ( 3._dp*zab2(j3,j1,j2,j4)*zab2(j2,j4,
     &    j3,j1)*zab2(j5,j1,j2,j6)*izab2(j1,j4,j3,j2)*s(j4,j3) + 3._dp*
     &    zab2(j3,j1,j2,j4)*zab2(j2,j4,j3,j1)*zab2(j5,j1,j2,j6)*izab2(
     &    j1,j4,j3,j2)*s(j1,j2) - 3._dp*zab2(j3,j1,j2,j4)*zab2(j2,j4,j3,
     &    j1)*zab2(j5,j1,j2,j6)*izab2(j1,j4,j3,j2)*s(j5,j6) - 3._dp*
     &    zab2(j3,j1,j2,j4)*zab2(j2,j6,j5,j1)*zab2(j5,j1,j2,j6)*izab2(
     &    j1,j6,j5,j2)*s(j4,j3) + 3._dp*zab2(j3,j1,j2,j4)*zab2(j2,j6,j5,
     &    j1)*zab2(j5,j1,j2,j6)*izab2(j1,j6,j5,j2)*s(j1,j2) + 3._dp*
     &    zab2(j3,j1,j2,j4)*zab2(j2,j6,j5,j1)*zab2(j5,j1,j2,j6)*izab2(
     &    j1,j6,j5,j2)*s(j5,j6) + 3._dp*zab2(j3,j2,j1,j4)*zab2(j2,j3,j4,
     &    j1)*zab2(j5,j2,j1,j6)*izab2(j1,j3,j4,j2)*s(j4,j3) + 3._dp*
     &    zab2(j3,j2,j1,j4)*zab2(j2,j3,j4,j1)*zab2(j5,j2,j1,j6)*izab2(
     &    j1,j3,j4,j2)*s(j1,j2) - 3._dp*zab2(j3,j2,j1,j4)*zab2(j2,j3,j4,
     &    j1)*zab2(j5,j2,j1,j6)*izab2(j1,j3,j4,j2)*s(j5,j6) - 3._dp*
     &    zab2(j3,j2,j1,j4)*zab2(j2,j5,j6,j1)*zab2(j5,j2,j1,j6)*izab2(
     &    j1,j5,j6,j2)*s(j4,j3) )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta**2 * ( 3._dp*zab2(j3,j2,j1,
     &    j4)*zab2(j2,j5,j6,j1)*zab2(j5,j2,j1,j6)*izab2(j1,j5,j6,j2)*s(
     &    j1,j2) + 3._dp*zab2(j3,j2,j1,j4)*zab2(j2,j5,j6,j1)*zab2(j5,j2,
     &    j1,j6)*izab2(j1,j5,j6,j2)*s(j5,j6) )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta * (  - 2._dp*za(j4,j3)*za(
     &    j3,j1)*za(j5,j6)*zb(j4,j3)*zb(j4,j6)*zb(j2,j6)*izab2(j1,j4,j3
     &    ,j2)**3*t(j4,j3,j1) + 2._dp*za(j4,j3)*za(j3,j1)*za(j5,j6)*zb(
     &    j4,j3)*zb(j4,j6)*zb(j2,j6)*izab2(j1,j4,j3,j2)**3*t(j4,j3,j2)
     &     - 2._dp*za(j4,j3)*za(j3,j5)*za(j1,j5)*zb(j4,j3)*zb(j4,j2)*zb(
     &    j5,j6)*izab2(j1,j3,j4,j2)**3*t(j3,j4,j1) + 2._dp*za(j4,j3)*za(
     &    j3,j5)*za(j1,j5)*zb(j4,j3)*zb(j4,j2)*zb(j5,j6)*izab2(j1,j3,j4
     &    ,j2)**3*t(j3,j4,j2) + za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(j4,j6)
     &    *izab2(j1,j4,j3,j2)**2*s(j4,j3) - za(j4,j3)*za(j3,j5)*zb(j4,
     &    j3)*zb(j4,j6)*izab2(j1,j4,j3,j2)**2*s(j1,j2) - za(j4,j3)*za(
     &    j3,j5)*zb(j4,j3)*zb(j4,j6)*izab2(j1,j4,j3,j2)**2*s(j5,j6) +
     &    za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(j4,j6)*izab2(j1,j3,j4,j2)**2
     &    *s(j4,j3) - za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(j4,j6)*izab2(j1,
     &    j3,j4,j2)**2*s(j1,j2) - za(j4,j3)*za(j3,j5)*zb(j4,j3)*zb(j4,
     &    j6)*izab2(j1,j3,j4,j2)**2*s(j5,j6) - 2._dp*za(j4,j3)*za(j1,j5)
     &    *za(j5,j6)*zb(j4,j2)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)**
     &    3*t(j6,j5,j1) )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta * ( 2._dp*za(j4,j3)*za(j1,
     &    j5)*za(j5,j6)*zb(j4,j2)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2
     &    )**3*t(j6,j5,j2) - 2._dp*za(j4,j3)*za(j1,j5)*zb(j4,j2)*zb(j4,
     &    j6)*izab2(j1,j3,j4,j2)**3*t(j3,j4,j1)*t(j3,j4,j2) + 2._dp*za(
     &    j4,j3)*za(j1,j5)*zb(j4,j2)*zb(j4,j6)*izab2(j1,j3,j4,j2)**3*t(
     &    j3,j4,j2)**2 + 2._dp*za(j4,j3)*za(j2,j5)*zb(j4,j2)*zb(j4,j6)*
     &    izab2(j1,j3,j4,j2)**2*t(j3,j4,j1) - 2._dp*za(j4,j3)*za(j2,j5)*
     &    zb(j4,j2)*zb(j4,j6)*izab2(j1,j3,j4,j2)**2*t(j3,j4,j2) - 1._dp/
     &    2._dp*za(j4,j3)*zb(j4,j6)**2*izb(j5,j6)*zab2(j2,j4,j3,j1)*
     &    izab2(j1,j4,j3,j2) - 1._dp/2._dp*za(j4,j3)*zb(j4,j6)**2*izb(j5,
     &    j6)*zab2(j2,j3,j4,j1)*izab2(j1,j3,j4,j2) + za(j4,j3)*zb(j4,j6
     &    )**2*izb(j5,j6)*izab2(j1,j3,j4,j2)**2*s(j4,j3)*t(j3,j4,j2) -
     &    za(j4,j3)*zb(j4,j6)**2*izb(j5,j6)*izab2(j1,j3,j4,j2)**2*s(j1,
     &    j2)*t(j3,j4,j2) - za(j4,j3)*zb(j4,j6)**2*izb(j5,j6)*izab2(j1,
     &    j3,j4,j2)**2*s(j5,j6)*t(j3,j4,j2) - 2._dp*za(j3,j1)*za(j3,j5)*
     &    za(j5,j6)*zb(j4,j3)*zb(j2,j6)*zb(j5,j6)*izab2(j1,j5,j6,j2)**3
     &    *t(j5,j6,j1) )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta * ( 2._dp*za(j3,j1)*za(j3,
     &    j5)*za(j5,j6)*zb(j4,j3)*zb(j2,j6)*zb(j5,j6)*izab2(j1,j5,j6,j2
     &    )**3*t(j5,j6,j2) + 2._dp*za(j3,j1)*za(j3,j5)*zb(j4,j3)*zb(j1,
     &    j6)*izab2(j1,j4,j3,j2)**2*t(j4,j3,j1) - 2._dp*za(j3,j1)*za(j3,
     &    j5)*zb(j4,j3)*zb(j1,j6)*izab2(j1,j4,j3,j2)**2*t(j4,j3,j2) - 2.
     &    _dp*za(j3,j1)*za(j3,j5)*zb(j4,j3)*zb(j2,j6)*izab2(j1,j4,j3,j2)
     &    **3*t(j4,j3,j1)**2 + 2._dp*za(j3,j1)*za(j3,j5)*zb(j4,j3)*zb(j2
     &    ,j6)*izab2(j1,j4,j3,j2)**3*t(j4,j3,j1)*t(j4,j3,j2) - 2._dp*za(
     &    j3,j1)*za(j5,j6)*zb(j4,j6)*zb(j2,j6)*izab2(j1,j5,j6,j2)**3*t(
     &    j5,j6,j1)*t(j5,j6,j2) + 2._dp*za(j3,j1)*za(j5,j6)*zb(j4,j6)*
     &    zb(j2,j6)*izab2(j1,j5,j6,j2)**3*t(j5,j6,j2)**2 + 2._dp*za(j3,
     &    j2)*za(j5,j6)*zb(j4,j6)*zb(j2,j6)*izab2(j1,j5,j6,j2)**2*t(j5,
     &    j6,j1) - 2._dp*za(j3,j2)*za(j5,j6)*zb(j4,j6)*zb(j2,j6)*izab2(
     &    j1,j5,j6,j2)**2*t(j5,j6,j2) - 1._dp/2._dp*za(j3,j5)**2*zb(j4,j3
     &    )*iza(j5,j6)*zab2(j2,j4,j3,j1)*izab2(j1,j4,j3,j2) - 1._dp/2._dp
     &    *za(j3,j5)**2*zb(j4,j3)*iza(j5,j6)*zab2(j2,j3,j4,j1)*izab2(j1
     &    ,j3,j4,j2) )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta * ( za(j3,j5)**2*zb(j4,j3)
     &    *iza(j5,j6)*izab2(j1,j4,j3,j2)**2*s(j4,j3)*t(j4,j3,j1) - za(
     &    j3,j5)**2*zb(j4,j3)*iza(j5,j6)*izab2(j1,j4,j3,j2)**2*s(j1,j2)
     &    *t(j4,j3,j1) - za(j3,j5)**2*zb(j4,j3)*iza(j5,j6)*izab2(j1,j4,
     &    j3,j2)**2*s(j5,j6)*t(j4,j3,j1) - 1._dp/2._dp*za(j3,j5)**2*zb(j5
     &    ,j6)*iza(j4,j3)*zab2(j2,j5,j6,j1)*izab2(j1,j5,j6,j2) - 1._dp/2.
     &    _dp*za(j3,j5)**2*zb(j5,j6)*iza(j4,j3)*zab2(j2,j6,j5,j1)*izab2(
     &    j1,j6,j5,j2) - za(j3,j5)**2*zb(j5,j6)*iza(j4,j3)*izab2(j1,j6,
     &    j5,j2)**2*s(j4,j3)*t(j6,j5,j1) - za(j3,j5)**2*zb(j5,j6)*iza(
     &    j4,j3)*izab2(j1,j6,j5,j2)**2*s(j1,j2)*t(j6,j5,j1) + za(j3,j5)
     &    **2*zb(j5,j6)*iza(j4,j3)*izab2(j1,j6,j5,j2)**2*s(j5,j6)*t(j6,
     &    j5,j1) + 2._dp*za(j3,j5)*za(j1,j5)*zb(j4,j1)*zb(j5,j6)*izab2(
     &    j1,j6,j5,j2)**2*t(j6,j5,j1) - 2._dp*za(j3,j5)*za(j1,j5)*zb(j4,
     &    j1)*zb(j5,j6)*izab2(j1,j6,j5,j2)**2*t(j6,j5,j2) - 2._dp*za(j3,
     &    j5)*za(j1,j5)*zb(j4,j2)*zb(j5,j6)*izab2(j1,j6,j5,j2)**3*t(j6,
     &    j5,j1)**2 )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta * ( 2._dp*za(j3,j5)*za(j1,
     &    j5)*zb(j4,j2)*zb(j5,j6)*izab2(j1,j6,j5,j2)**3*t(j6,j5,j1)*t(
     &    j6,j5,j2) - za(j3,j5)*za(j5,j6)*zb(j4,j6)*zb(j5,j6)*izab2(j1,
     &    j5,j6,j2)**2*s(j4,j3) - za(j3,j5)*za(j5,j6)*zb(j4,j6)*zb(j5,
     &    j6)*izab2(j1,j5,j6,j2)**2*s(j1,j2) + za(j3,j5)*za(j5,j6)*zb(
     &    j4,j6)*zb(j5,j6)*izab2(j1,j5,j6,j2)**2*s(j5,j6) - za(j3,j5)*
     &    za(j5,j6)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)**2*s(j4,j3)
     &     - za(j3,j5)*za(j5,j6)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)
     &    **2*s(j1,j2) + za(j3,j5)*za(j5,j6)*zb(j4,j6)*zb(j5,j6)*izab2(
     &    j1,j6,j5,j2)**2*s(j5,j6) + za(j3,j5)*zb(j4,j6)*zab2(j2,j4,j3,
     &    j1)*izab2(j1,j4,j3,j2) + za(j3,j5)*zb(j4,j6)*zab2(j2,j3,j4,j1
     &    )*izab2(j1,j3,j4,j2) + za(j3,j5)*zb(j4,j6)*zab2(j2,j5,j6,j1)*
     &    izab2(j1,j5,j6,j2) + za(j3,j5)*zb(j4,j6)*zab2(j2,j6,j5,j1)*
     &    izab2(j1,j6,j5,j2) - 1._dp/2._dp*za(j5,j6)*zb(j4,j6)**2*izb(j4,
     &    j3)*zab2(j2,j5,j6,j1)*izab2(j1,j5,j6,j2) - 1._dp/2._dp*za(j5,j6
     &    )*zb(j4,j6)**2*izb(j4,j3)*zab2(j2,j6,j5,j1)*izab2(j1,j6,j5,j2
     &    ) )
      bcoeff(b12zm) = bcoeff(b12zm) + IDelta * (  - za(j5,j6)*zb(j4,j6)
     &    **2*izb(j4,j3)*izab2(j1,j5,j6,j2)**2*s(j4,j3)*t(j5,j6,j2) -
     &    za(j5,j6)*zb(j4,j6)**2*izb(j4,j3)*izab2(j1,j5,j6,j2)**2*s(j1,
     &    j2)*t(j5,j6,j2) + za(j5,j6)*zb(j4,j6)**2*izb(j4,j3)*izab2(j1,
     &    j5,j6,j2)**2*s(j5,j6)*t(j5,j6,j2) )
      bcoeff(b12zm) = bcoeff(b12zm) - 2._dp*za(j3,j1)*za(j3,j2)*zb(j1,j6
     & )*zb(j2,j6)*iza(j4,j3)*izb(j5,j6)*izab2(j1,j4,j3,j2)**2 + 2._dp*
     &    za(j3,j1)*zb(j2,j6)*iza(j4,j3)*izb(j5,j6)*zab2(j3,j4,j2,j6)*
     &    izab2(j1,j4,j3,j2)**3*t(j4,j3,j1) - 2._dp*za(j1,j5)*za(j2,j5)*
     &    zb(j4,j1)*zb(j4,j2)*iza(j5,j6)*izb(j4,j3)*izab2(j1,j3,j4,j2)
     &    **2 + 2._dp*za(j1,j5)*zb(j4,j2)*iza(j5,j6)*izb(j4,j3)*zab2(j5,
     &    j3,j1,j4)*izab2(j1,j3,j4,j2)**3*t(j3,j4,j2) - iza(j4,j3)*izb(
     &    j5,j6)*zab2(j3,j4,j2,j6)**2*izab2(j1,j4,j3,j2)**2 - iza(j5,j6
     &    )*izb(j4,j3)*zab2(j5,j3,j1,j4)**2*izab2(j1,j3,j4,j2)**2

      bcoeff(b56)= + t234ms56**(-2) * ( za(j1,j5)**2*zb(j4,j1)**2*zb(j5
     &    ,j6)*izb(j4,j3)*izab2(j1,j6,j5,j2)**2*s(j5,j6) )
      bcoeff(b56) = bcoeff(b56) + t234ms56**(-1) * ( 2._dp*za(j1,j5)**2*
     &    zb(j4,j1)**2*zb(j5,j6)*izb(j4,j3)*izab2(j1,j6,j5,j2)**2 - 2._dp
     &    *za(j1,j5)*zb(j4,j1)*zb(j5,j6)*izb(j4,j3)*zab2(j1,j6,j5,j4)*
     &    zab2(j5,j6,j1,j2)*izab2(j1,j6,j5,j2)**3 )
      bcoeff(b56) = bcoeff(b56) + t134ms56**(-2) * ( za(j3,j2)**2*za(j5
     &    ,j6)*zb(j2,j6)**2*iza(j4,j3)*izab2(j1,j5,j6,j2)**2*s(j5,j6) )
      bcoeff(b56) = bcoeff(b56) + t134ms56**(-1) * ( 2._dp*za(j3,j2)**2*
     &    za(j5,j6)*zb(j2,j6)**2*iza(j4,j3)*izab2(j1,j5,j6,j2)**2 - 2._dp
     &    *za(j3,j2)*za(j5,j6)*zb(j2,j6)*iza(j4,j3)*zab2(j3,j5,j6,j2)*
     &    zab2(j1,j5,j2,j6)*izab2(j1,j5,j6,j2)**3 )
      bcoeff(b56) = bcoeff(b56) + IDelta**2 * ( 3._dp*zab2(j3,j1,j2,j4)*
     &    zab2(j2,j6,j5,j1)*zab2(j5,j1,j2,j6)*izab2(j1,j6,j5,j2)*s(j4,
     &    j3) - 3._dp*zab2(j3,j1,j2,j4)*zab2(j2,j6,j5,j1)*zab2(j5,j1,j2,
     &    j6)*izab2(j1,j6,j5,j2)*s(j1,j2) - 3._dp*zab2(j3,j1,j2,j4)*
     &    zab2(j2,j6,j5,j1)*zab2(j5,j1,j2,j6)*izab2(j1,j6,j5,j2)*s(j5,
     &    j6) + 3._dp*zab2(j3,j2,j1,j4)*zab2(j2,j5,j6,j1)*zab2(j5,j2,j1,
     &    j6)*izab2(j1,j5,j6,j2)*s(j4,j3) - 3._dp*zab2(j3,j2,j1,j4)*
     &    zab2(j2,j5,j6,j1)*zab2(j5,j2,j1,j6)*izab2(j1,j5,j6,j2)*s(j1,
     &    j2) - 3._dp*zab2(j3,j2,j1,j4)*zab2(j2,j5,j6,j1)*zab2(j5,j2,j1,
     &    j6)*izab2(j1,j5,j6,j2)*s(j5,j6) )
      bcoeff(b56) = bcoeff(b56) + IDelta * ( 2._dp*za(j4,j3)*za(j1,j5)*
     &    za(j5,j6)*zb(j4,j2)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)**3
     &    *t(j6,j5,j1) - 2._dp*za(j4,j3)*za(j1,j5)*za(j5,j6)*zb(j4,j2)*
     &    zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)**3*t(j6,j5,j2) + 2._dp*
     &    za(j3,j1)*za(j3,j5)*za(j5,j6)*zb(j4,j3)*zb(j2,j6)*zb(j5,j6)*
     &    izab2(j1,j5,j6,j2)**3*t(j5,j6,j1) - 2._dp*za(j3,j1)*za(j3,j5)*
     &    za(j5,j6)*zb(j4,j3)*zb(j2,j6)*zb(j5,j6)*izab2(j1,j5,j6,j2)**3
     &    *t(j5,j6,j2) + 2._dp*za(j3,j1)*za(j5,j6)*zb(j4,j6)*zb(j2,j6)*
     &    izab2(j1,j5,j6,j2)**3*t(j5,j6,j1)*t(j5,j6,j2) - 2._dp*za(j3,j1
     &    )*za(j5,j6)*zb(j4,j6)*zb(j2,j6)*izab2(j1,j5,j6,j2)**3*t(j5,j6
     &    ,j2)**2 - 2._dp*za(j3,j2)*za(j5,j6)*zb(j4,j6)*zb(j2,j6)*izab2(
     &    j1,j5,j6,j2)**2*t(j5,j6,j1) + 2._dp*za(j3,j2)*za(j5,j6)*zb(j4,
     &    j6)*zb(j2,j6)*izab2(j1,j5,j6,j2)**2*t(j5,j6,j2) + 1._dp/2._dp*
     &    za(j3,j5)**2*zb(j5,j6)*iza(j4,j3)*zab2(j2,j5,j6,j1)*izab2(j1,
     &    j5,j6,j2) + 1._dp/2._dp*za(j3,j5)**2*zb(j5,j6)*iza(j4,j3)*zab2(
     &    j2,j6,j5,j1)*izab2(j1,j6,j5,j2) )
      bcoeff(b56) = bcoeff(b56) + IDelta * ( za(j3,j5)**2*zb(j5,j6)*
     &    iza(j4,j3)*izab2(j1,j6,j5,j2)**2*s(j4,j3)*t(j6,j5,j1) + za(j3
     &    ,j5)**2*zb(j5,j6)*iza(j4,j3)*izab2(j1,j6,j5,j2)**2*s(j1,j2)*
     &    t(j6,j5,j1) - za(j3,j5)**2*zb(j5,j6)*iza(j4,j3)*izab2(j1,j6,
     &    j5,j2)**2*s(j5,j6)*t(j6,j5,j1) - 2._dp*za(j3,j5)*za(j1,j5)*zb(
     &    j4,j1)*zb(j5,j6)*izab2(j1,j6,j5,j2)**2*t(j6,j5,j1) + 2._dp*za(
     &    j3,j5)*za(j1,j5)*zb(j4,j1)*zb(j5,j6)*izab2(j1,j6,j5,j2)**2*t(
     &    j6,j5,j2) + 2._dp*za(j3,j5)*za(j1,j5)*zb(j4,j2)*zb(j5,j6)*
     &    izab2(j1,j6,j5,j2)**3*t(j6,j5,j1)**2 - 2._dp*za(j3,j5)*za(j1,
     &    j5)*zb(j4,j2)*zb(j5,j6)*izab2(j1,j6,j5,j2)**3*t(j6,j5,j1)*t(
     &    j6,j5,j2) + za(j3,j5)*za(j5,j6)*zb(j4,j6)*zb(j5,j6)*izab2(j1,
     &    j5,j6,j2)**2*s(j4,j3) + za(j3,j5)*za(j5,j6)*zb(j4,j6)*zb(j5,
     &    j6)*izab2(j1,j5,j6,j2)**2*s(j1,j2) - za(j3,j5)*za(j5,j6)*zb(
     &    j4,j6)*zb(j5,j6)*izab2(j1,j5,j6,j2)**2*s(j5,j6) + za(j3,j5)*
     &    za(j5,j6)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)**2*s(j4,j3)
     &     + za(j3,j5)*za(j5,j6)*zb(j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)
     &    **2*s(j1,j2) )
      bcoeff(b56) = bcoeff(b56) + IDelta * (  - za(j3,j5)*za(j5,j6)*zb(
     &    j4,j6)*zb(j5,j6)*izab2(j1,j6,j5,j2)**2*s(j5,j6) - za(j3,j5)*
     &    zb(j4,j6)*zab2(j2,j5,j6,j1)*izab2(j1,j5,j6,j2) - za(j3,j5)*
     &    zb(j4,j6)*zab2(j2,j6,j5,j1)*izab2(j1,j6,j5,j2) + 1._dp/2._dp*
     &    za(j5,j6)*zb(j4,j6)**2*izb(j4,j3)*zab2(j2,j5,j6,j1)*izab2(j1,
     &    j5,j6,j2) + 1._dp/2._dp*za(j5,j6)*zb(j4,j6)**2*izb(j4,j3)*zab2(
     &    j2,j6,j5,j1)*izab2(j1,j6,j5,j2) + za(j5,j6)*zb(j4,j6)**2*izb(
     &    j4,j3)*izab2(j1,j5,j6,j2)**2*s(j4,j3)*t(j5,j6,j2) + za(j5,j6)
     &    *zb(j4,j6)**2*izb(j4,j3)*izab2(j1,j5,j6,j2)**2*s(j1,j2)*t(j5,
     &    j6,j2) - za(j5,j6)*zb(j4,j6)**2*izb(j4,j3)*izab2(j1,j5,j6,j2)
     &    **2*s(j5,j6)*t(j5,j6,j2) )

      bcoeff(b134)= + t134ms34**(-2) * (  - za(j3,j1)**2*zb(j4,j3)*zb(
     &    j1,j6)**2*izb(j5,j6)*izab2(j1,j4,j3,j2)**2*s(j4,j3) )
      bcoeff(b134) = bcoeff(b134) + t134ms34**(-1) * (  - 2._dp*za(j3,j1
     &    )**2*zb(j4,j3)*zb(j1,j6)**2*izb(j5,j6)*izab2(j1,j4,j3,j2)**2
     &     + 2._dp*za(j3,j1)*zb(j4,j3)*zb(j1,j6)*izb(j5,j6)*zab2(j3,j4,
     &    j1,j2)*zab2(j1,j4,j3,j6)*izab2(j1,j4,j3,j2)**3 )
      bcoeff(b134) = bcoeff(b134) + t134ms56**(-2) * (  - za(j3,j2)**2*
     &    za(j5,j6)*zb(j2,j6)**2*iza(j4,j3)*izab2(j1,j5,j6,j2)**2*s(j5,
     &    j6) )
      bcoeff(b134) = bcoeff(b134) + t134ms56**(-1) * (  - 2._dp*za(j3,j2
     &    )**2*za(j5,j6)*zb(j2,j6)**2*iza(j4,j3)*izab2(j1,j5,j6,j2)**2
     &     + 2._dp*za(j3,j2)*za(j5,j6)*zb(j2,j6)*iza(j4,j3)*zab2(j3,j5,
     &    j6,j2)*zab2(j1,j5,j2,j6)*izab2(j1,j5,j6,j2)**3 )
      bcoeff(b134) = bcoeff(b134) + 2._dp*za(j3,j1)*za(j3,j2)*zb(j1,j6)*
     & zb(j2,j6)*iza(j4,j3)*izb(j5,j6)*izab2(j1,j4,j3,j2)**2 - 2._dp*za(
     &    j3,j1)*zb(j2,j6)*iza(j4,j3)*izb(j5,j6)*zab2(j3,j4,j2,j6)*
     &    izab2(j1,j4,j3,j2)**3*t(j4,j3,j1) + iza(j4,j3)*izb(j5,j6)*
     &    zab2(j3,j4,j2,j6)**2*izab2(j1,j4,j3,j2)**2

      bcoeff(b234)= + t234ms34**(-2) * (  - za(j4,j3)*za(j2,j5)**2*zb(
     &    j4,j2)**2*iza(j5,j6)*izab2(j1,j3,j4,j2)**2*s(j4,j3) )
      bcoeff(b234) = bcoeff(b234) + t234ms34**(-1) * (  - 2._dp*za(j4,j3
     &    )*za(j2,j5)**2*zb(j4,j2)**2*iza(j5,j6)*izab2(j1,j3,j4,j2)**2
     &     + 2._dp*za(j4,j3)*za(j2,j5)*zb(j4,j2)*iza(j5,j6)*zab2(j1,j3,
     &    j2,j4)*zab2(j5,j3,j4,j2)*izab2(j1,j3,j4,j2)**3 )
      bcoeff(b234) = bcoeff(b234) + t234ms56**(-2) * (  - za(j1,j5)**2*
     &    zb(j4,j1)**2*zb(j5,j6)*izb(j4,j3)*izab2(j1,j6,j5,j2)**2*s(j5,
     &    j6) )
      bcoeff(b234) = bcoeff(b234) + t234ms56**(-1) * (  - 2._dp*za(j1,j5
     &    )**2*zb(j4,j1)**2*zb(j5,j6)*izb(j4,j3)*izab2(j1,j6,j5,j2)**2
     &     + 2._dp*za(j1,j5)*zb(j4,j1)*zb(j5,j6)*izb(j4,j3)*zab2(j1,j6,
     &    j5,j4)*zab2(j5,j6,j1,j2)*izab2(j1,j6,j5,j2)**3 )
      bcoeff(b234) = bcoeff(b234) + 2._dp*za(j1,j5)*za(j2,j5)*zb(j4,j1)*
     & zb(j4,j2)*iza(j5,j6)*izb(j4,j3)*izab2(j1,j3,j4,j2)**2 - 2._dp*za(
     &    j1,j5)*zb(j4,j2)*iza(j5,j6)*izb(j4,j3)*zab2(j5,j3,j1,j4)*
     &    izab2(j1,j3,j4,j2)**3*t(j3,j4,j2) + iza(j5,j6)*izb(j4,j3)*
     &    zab2(j5,j3,j1,j4)**2*izab2(j1,j3,j4,j2)**2

      bcoeff(rat)= + t134ms34**(-1) * (  - za(j3,j1)**2*zb(j4,j3)*zb(j1
     &    ,j6)**2*izb(j5,j6)*izab2(j1,j4,j3,j2)**2 )
      bcoeff(rat) = bcoeff(rat) + t234ms34**(-1) * (  - za(j4,j3)*za(j2
     &    ,j5)**2*zb(j4,j2)**2*iza(j5,j6)*izab2(j1,j3,j4,j2)**2 )
      bcoeff(rat) = bcoeff(rat) + t234ms56**(-1) * (  - za(j1,j5)**2*
     &    zb(j4,j1)**2*zb(j5,j6)*izb(j4,j3)*izab2(j1,j6,j5,j2)**2 )
      bcoeff(rat) = bcoeff(rat) + t134ms56**(-1) * (  - za(j3,j2)**2*
     &    za(j5,j6)*zb(j2,j6)**2*iza(j4,j3)*izab2(j1,j5,j6,j2)**2 )
      bcoeff(rat) = bcoeff(rat) + IDelta * (  - za(j3,j5)**2*iza(j4,j3)
     &    *iza(j5,j6)*zab2(j2,j3,j4,j1)*izab2(j1,j3,j4,j2)*s(j4,j3) +
     &    za(j3,j5)**2*iza(j4,j3)*iza(j5,j6)*zab2(j2,j3,j4,j1)*izab2(j1
     &    ,j3,j4,j2)*s(j1,j2) - za(j3,j5)**2*iza(j4,j3)*iza(j5,j6)*
     &    zab2(j2,j3,j4,j1)*izab2(j1,j3,j4,j2)*s(j5,j6) - 2._dp*za(j3,j5
     &    )*zb(j4,j6)*zab2(j2,j4,j3,j1)*izab2(j1,j4,j3,j2) - 2._dp*za(j3
     &    ,j5)*zb(j4,j6)*zab2(j2,j3,j4,j1)*izab2(j1,j3,j4,j2) - zb(j4,
     &    j6)**2*izb(j4,j3)*izb(j5,j6)*zab2(j2,j4,j3,j1)*izab2(j1,j4,j3
     &    ,j2)*s(j4,j3) + zb(j4,j6)**2*izb(j4,j3)*izb(j5,j6)*zab2(j2,j4
     &    ,j3,j1)*izab2(j1,j4,j3,j2)*s(j1,j2) - zb(j4,j6)**2*izb(j4,j3)
     &    *izb(j5,j6)*zab2(j2,j4,j3,j1)*izab2(j1,j4,j3,j2)*s(j5,j6) )
      bcoeff(rat) = bcoeff(rat) - za(j3,j5)*zb(j4,j6)*izab2(j1,j4,j3,j2
     & )**2 - za(j3,j5)*zb(j4,j6)*izab2(j1,j3,j4,j2)**2 + iza(j4,j3)*
     &    izb(j5,j6)*zab2(j3,j4,j2,j6)**2*izab2(j1,j4,j3,j2)**2 + iza(
     &    j5,j6)*izb(j4,j3)*zab2(j5,j3,j1,j4)**2*izab2(j1,j3,j4,j2)**2

      if     ((kcase==kHWWint) .or. (kcase==kWWqqbr)
     &   .or. (kcase==kHWWHpi) .or. (kcase==kggWW4l)
     &   .or. (kcase==kggWWbx)
     &   .or. (kcase==kggVV4l) .or. (kcase==kggVVbx)) then
      bcoeff(b12)=bub6(j2,j1,j4,j3,j6,j5,zb,za)
      bcoeff(b12m)=bcoeff(b12zm)-bcoeff(b12)
      elseif ((kcase==kHZZint) .or. (kcase==kZZlept)
     &   .or. (kcase==kHZZHpi) .or. (kcase==kggZZ4l)
     &   .or. (kcase==kggZZbx)) then
      bcoeff(b12)=bcoeff(b12zm)
      bcoeff(b1)=-bcoeff(b12)-bcoeff(b34)-bcoeff(b56)
     & -bcoeff(b134)-bcoeff(b234)
      else
      write(6,*) 'Unimplemented case in mbc.f'
      stop
      endif

      endif

      do j=1,8
      bcoeff(j)=im*bcoeff(j)
      enddo

      return
      end




