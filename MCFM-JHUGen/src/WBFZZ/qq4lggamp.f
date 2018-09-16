      subroutine qq4lggamp(i1,i2,i3,i4,i5,i6,i7,i8,
     & b7,b8,za,zb,amp78,amp87)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      integer:: jdu,j78,h34,h56,i1,i2,i3,i4,i5,i6,i7,i8,b7,b8,
     & p1,p2,p3,p4,p5,p6,p7,p8
      complex(dp):: zab2,zba2,
     & amp78(2,2,2,2,2,2),amp87(2,2,2,2,2,2),
     & bmp78(2,2,2,2,2,2),bmp87(2,2,2,2,2,2),
     & propz34,propz56,propz3456,gamz3456(2,2,2),gamll56(2,2),
     & zab3,zba3,zab4,zba4,iza,izb,gamz34(2,2,2),gamz56(2,2,2)
      real(dp):: t3,t4,s356,s456,
     & s34,s56,s18,s27,s78,s278,s178,s156,s234,s2347,s3456,
     & xq1,xq2,xl1,xl2,xr1,xr2
C-----Begin statement functions
      iza(i1,i2)=cone/za(i1,i2)
      izb(i1,i2)=cone/zb(i1,i2)
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      zba2(i1,i2,i3,i4)=zb(i1,i2)*za(i2,i4)+zb(i1,i3)*za(i3,i4)
      zab3(i1,i2,i3,i4,i5)=
     & za(i1,i2)*zb(i2,i5)+za(i1,i3)*zb(i3,i5)+za(i1,i4)*zb(i4,i5)
      zba3(i1,i2,i3,i4,i5)=
     & zb(i1,i2)*za(i2,i5)+zb(i1,i3)*za(i3,i5)+zb(i1,i4)*za(i4,i5)
      zab4(i1,i2,i3,i4,i5,i6)=
     & +za(i1,i2)*zb(i2,i6)+za(i1,i3)*zb(i3,i6)
     & +za(i1,i4)*zb(i4,i6)+za(i1,i5)*zb(i5,i6)
      zba4(i1,i2,i3,i4,i5,i6)=
     & +zb(i1,i2)*za(i2,i6)+zb(i1,i3)*za(i3,i6)
     & +zb(i1,i4)*za(i4,i6)+zb(i1,i5)*za(i5,i6)
      t3(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i3,i1)
      t4(i1,i2,i3,i4)=s(i1,i2)+s(i1,i3)+s(i1,i4)
     &             +s(i2,i3)+s(i2,i4)+s(i3,i4)
C-----end statement functions

      p1=i1
      p2=i2
      s34=s(i3,i4)
      s56=s(i5,i6)
      s278=t3(i2,i7,i8)
      s178=t3(i1,i7,i8)
      s156=t3(i1,i5,i6)
      s234=t3(i2,i3,i4)
      s356=t3(i3,i5,i6)
      s456=t3(i4,i5,i6)
      s3456=t4(i3,i4,i5,i6)

      propz3456=s3456-czmass2
      propz34=s34-czmass2
      propz56=s56-czmass2

C---set up couplings dependent on whether we are doing 34- or 56- line
      if (i3+i4 == 7) then
      xl1=l1
      xr1=r1
      xq1=q1
      xl2=l2
      xr2=r2
      xq2=q2
      elseif (i3+i4 == 11) then
      xl1=l2
      xr1=r2
      xq1=q2
      xl2=l1
      xr2=r1
      xq2=q1
      else
      write(6,*) 'Unexpected case qq4lggamp.f'
      stop
      endif

      do jdu=1,2
      gamz34(jdu,1,1)=Q(jdu)*xq1/s34+L(jdu)*xl1/propz34
      gamz34(jdu,1,2)=Q(jdu)*xq1/s34+L(jdu)*xr1/propz34
      gamz34(jdu,2,1)=Q(jdu)*xq1/s34+R(jdu)*xl1/propz34
      gamz34(jdu,2,2)=Q(jdu)*xq1/s34+R(jdu)*xr1/propz34
      gamz3456(jdu,1,1)=Q(jdu)*xq1/s3456+L(jdu)*xl1/propz3456
      gamz3456(jdu,1,2)=Q(jdu)*xq1/s3456+L(jdu)*xr1/propz3456
      gamz3456(jdu,2,1)=Q(jdu)*xq1/s3456+R(jdu)*xl1/propz3456
      gamz3456(jdu,2,2)=Q(jdu)*xq1/s3456+R(jdu)*xr1/propz3456
      gamz56(jdu,1,1)=Q(jdu)*xq2/s56+L(jdu)*xl2/propz56
      gamz56(jdu,1,2)=Q(jdu)*xq2/s56+L(jdu)*xr2/propz56
      gamz56(jdu,2,1)=Q(jdu)*xq2/s56+R(jdu)*xl2/propz56
      gamz56(jdu,2,2)=Q(jdu)*xq2/s56+R(jdu)*xr2/propz56
      gamll56(1,1)=xq1*xq2/s56+xl1*xl2/propz56
      gamll56(1,2)=xq1*xq2/s56+xl1*xr2/propz56
      gamll56(2,1)=xq1*xq2/s56+xr1*xl2/propz56
      gamll56(2,2)=xq2*xq2/s56+xr1*xr2/propz56
      enddo

      do j78=1,2
      if (j78 == 1) then
      p7=i7
      p8=i8
      elseif (j78 == 2) then
      p7=i8
      p8=i7
      endif

      s78=s(p7,p8)
      s27=s(p2,p7)
      s18=s(p1,p8)
      s2347=t4(i2,i3,i4,p7)

      do jdu=1,2
      do h56=1,2
      if (h56 == 1) then
      p5=i5
      p6=i6
      elseif (h56 == 2) then
      p5=i6
      p6=i5
      endif
      p3=i3
      p4=i4
      if (j78==1) then
C---  bmp78(jdu,h12,h34,h56,h7,h8)
      bmp78(jdu,1,1,h56,1,1)= + gamz3456(jdu,1,1)*gamll56(1,h56) * ( 4.D
     &    0*za(p2,p7)*za(p3,p5)*za(p7,p8)*zb(p1,p4)*zb(p3,p6)*izb(p1,p8
     &    )*zba3(p1,p2,p7,p8,p3)*s78**(-1)*s278**(-1)*s356**(-1) + 4.D0
     &    *za(p2,p7)*za(p3,p5)*za(p7,p8)*zb(p1,p4)*zb(p5,p6)*izb(p1,p8)
     &    *zba3(p1,p2,p7,p8,p5)*s78**(-1)*s278**(-1)*s356**(-1) + 4.D0*
     &    za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p3,p6)*zba2(p1,p2,p7,p8)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*s27**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p2,p7)*za(p3,p5)*zb(p1,p4)*
     &    zb(p5,p6)*zba2(p1,p2,p7,p8)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,
     &    p7,p8,p5)*s27**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p2,p7)*
     &    za(p4,p5)*za(p7,p8)*zb(p1,p4)*zb(p4,p6)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3)*s78**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*
     &    za(p4,p5)*zb(p1,p4)*zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p1,p7)*
     &    izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*s27**(-1)*s278**(-1)*
     &    s456**(-1) + 4.D0*za(p2,p7)*za(p5,p6)*za(p7,p8)*zb(p1,p6)*zb(
     &    p4,p6)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*s78**(-1)*s278**(-1)*
     &    s456**(-1) )
      bmp78(jdu,1,1,h56,1,1) = bmp78(jdu,1,1,h56,1,1) + gamz3456(jdu,1,
     & 1)*gamll56(1,h56) * ( 4.D0*za(p2,p7)*za(p5,p6)*zb(p1,p6)*zb(p4,
     &    p6)*zba2(p1,p2,p7,p8)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,
     &    p3)*s27**(-1)*s278**(-1)*s456**(-1) + 4.D0*za(p2,p8)*za(p3,p5
     &    )*za(p7,p8)*zb(p1,p4)*zb(p3,p6)*izb(p1,p7)*zba3(p1,p2,p7,p8,
     &    p3)*s78**(-1)*s278**(-1)*s356**(-1) + 4.D0*za(p2,p8)*za(p3,p5
     &    )*za(p7,p8)*zb(p1,p4)*zb(p5,p6)*izb(p1,p7)*zba3(p1,p2,p7,p8,
     &    p5)*s78**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p2,p8)*za(p4,p5
     &    )*za(p7,p8)*zb(p1,p4)*zb(p4,p6)*izb(p1,p7)*zba3(p1,p2,p7,p8,
     &    p3)*s78**(-1)*s278**(-1)*s456**(-1) + 4.D0*za(p2,p8)*za(p5,p6
     &    )*za(p7,p8)*zb(p1,p6)*zb(p4,p6)*izb(p1,p7)*zba3(p1,p2,p7,p8,
     &    p3)*s78**(-1)*s278**(-1)*s456**(-1) )
      bmp78(jdu,1,1,h56,2,1)= + gamz3456(jdu,1,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p3)*za(p2,p8)*za(p1,p8)*za(p3,p5)*zb(p2,p7)*zb(
     &    p1,p4)*zb(p1,p7)*zb(p3,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s27**(-1)*s356**(-1) + 4.D0*za(p2,p3)*za(p2,p8)*za(p1,p8)*za(
     &    p4,p5)*zb(p2,p7)*zb(p1,p4)*zb(p1,p7)*zb(p4,p6)*iza(p7,p8)*
     &    izb(p7,p8)*s18**(-1)*s27**(-1)*s456**(-1) - 4.D0*za(p2,p3)*
     &    za(p2,p8)*za(p1,p8)*za(p5,p6)*zb(p2,p7)*zb(p1,p6)*zb(p1,p7)*
     &    zb(p4,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s456**(-1) - 4.D0*za(p2,p3)*za(p1,p8)*za(p3,p5)*zb(p1,p7)**2*
     &    zb(p3,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p3)*za(p1,p8)*za(p4,p5)*
     &    zb(p1,p7)**2*zb(p4,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p7,p8
     &    )*s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p3)*za(p1,p8)*
     &    za(p5,p6)*zb(p1,p7)**2*zb(p4,p6)*zba2(p6,p1,p7,p8)*iza(p7,p8)
     &    *izb(p7,p8)*s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p5)*
     &    za(p2,p8)*za(p1,p8)*za(p3,p5)*zb(p2,p7)*zb(p1,p4)*zb(p1,p7)*
     &    zb(p5,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s356**(-1) )
      bmp78(jdu,1,1,h56,2,1) = bmp78(jdu,1,1,h56,2,1) + gamz3456(jdu,1,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p2,p5)*za(p1,p8)*za(p3,p5)*zb(
     &    p1,p7)**2*zb(p5,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s178**(-1)*s356**(-1) - 4.D0*za(p2,p8)**2*za(p3,p5)
     &    *zb(p2,p7)*zb(p1,p4)*zb(p3,p6)*zba2(p7,p2,p8,p3)*iza(p7,p8)*
     &    izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p2,p8)**
     &    2*za(p3,p5)*zb(p2,p7)*zb(p1,p4)*zb(p5,p6)*zba2(p7,p2,p8,p5)*
     &    iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0*
     &    za(p2,p8)**2*za(p4,p5)*zb(p2,p7)*zb(p1,p4)*zb(p4,p6)*zba2(p7,
     &    p2,p8,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s456**(-1) - 4.D0*za(p2,p8)**2*za(p5,p6)*zb(p2,p7)*zb(p1,p6)*
     &    zb(p4,p6)*zba2(p7,p2,p8,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) )
      bmp78(jdu,1,1,h56,1,2)= + gamz3456(jdu,1,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p3)*za(p1,p7)*za(p3,p5)*zb(p1,p8)**2*zb(p3,p6)*
     &    zba2(p4,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p3)*za(p1,p7)*za(p4,p5)*zb(p1,p8)**2*
     &    zb(p4,p6)*zba2(p4,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p2,p3)*za(p1,p7)*za(p5,p6)*
     &    zb(p1,p8)**2*zb(p4,p6)*zba2(p6,p1,p8,p7)*iza(p7,p8)*izb(p7,p8
     &    )*s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p5)*za(p1,p7)*
     &    za(p3,p5)*zb(p1,p8)**2*zb(p5,p6)*zba2(p4,p1,p8,p7)*iza(p7,p8)
     &    *izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) - 4.D0*za(p2,p7)
     &    **2*za(p3,p5)*zb(p2,p8)*zb(p1,p4)*zb(p3,p6)*zba2(p8,p2,p7,p3)
     &    *iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) - 4.D0
     &    *za(p2,p7)**2*za(p3,p5)*zb(p2,p8)*zb(p1,p4)*zb(p5,p6)*zba2(p8
     &    ,p2,p7,p5)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p7)**2*za(p4,p5)*zb(p2,p8)*zb(p1,p4)*
     &    zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) )
      bmp78(jdu,1,1,h56,1,2) = bmp78(jdu,1,1,h56,1,2) + gamz3456(jdu,1,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p2,p7)**2*za(p5,p6)*zb(p2,p8)*
     &    zb(p1,p6)*zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*
     &    s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p3,p5)*
     &    zb(p1,p8)*zb(p3,p6)*zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(
     &    p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s356**(-1) - 4.D0*za(p2
     &    ,p7)*za(p3,p5)*zb(p1,p8)*zb(p5,p6)*zba2(p4,p1,p8,p7)*zba2(p8,
     &    p2,p7,p5)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p7)*za(p4,p5)*zb(p1,p8)*zb(p4,p6)*
     &    zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s27**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p5,p6)*zb(
     &    p1,p8)*zb(p4,p6)*zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p7,
     &    p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s456**(-1) )
      bmp78(jdu,1,1,h56,2,2)= + gamz3456(jdu,1,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(p3,p6)*zb(p7,p8)*
     &    iza(p2,p8)*zba3(p4,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*
     &    s356**(-1) - 4.D0*za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*zb(
     &    p7,p8)*iza(p2,p7)*zba3(p4,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*
     &    zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*zba3(p4,p1,p7,p8,p2)*
     &    s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p3)*za(p4,p5)*
     &    zb(p1,p7)*zb(p4,p6)*zb(p7,p8)*iza(p2,p8)*zba3(p4,p1,p7,p8,p2)
     &    *s78**(-1)*s178**(-1)*s456**(-1) + 4.D0*za(p2,p3)*za(p4,p5)*
     &    zb(p1,p8)*zb(p4,p6)*zb(p7,p8)*iza(p2,p7)*zba3(p4,p1,p7,p8,p2)
     &    *s78**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p3)*za(p4,p5)*
     &    zb(p1,p8)*zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*
     &    zba3(p4,p1,p7,p8,p2)*s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*
     &    za(p2,p3)*za(p5,p6)*zb(p1,p7)*zb(p4,p6)*zb(p7,p8)*iza(p2,p8)*
     &    zba3(p6,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*s456**(-1) )
      bmp78(jdu,1,1,h56,2,2) = bmp78(jdu,1,1,h56,2,2) + gamz3456(jdu,1,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p2,p3)*za(p5,p6)*zb(p1,p8)*zb(
     &    p4,p6)*zb(p7,p8)*iza(p2,p7)*zba3(p6,p1,p7,p8,p2)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p2,p3)*za(p5,p6)*zb(p1,p8)*
     &    zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*zba3(p6,p1,
     &    p7,p8,p2)*s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p5)*
     &    za(p3,p5)*zb(p1,p7)*zb(p5,p6)*zb(p7,p8)*iza(p2,p8)*zba3(p4,p1
     &    ,p7,p8,p2)*s78**(-1)*s178**(-1)*s356**(-1) - 4.D0*za(p2,p5)*
     &    za(p3,p5)*zb(p1,p8)*zb(p5,p6)*zb(p7,p8)*iza(p2,p7)*zba3(p4,p1
     &    ,p7,p8,p2)*s78**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p5)*
     &    za(p3,p5)*zb(p1,p8)*zb(p5,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*
     &    iza(p2,p8)*zba3(p4,p1,p7,p8,p2)*s18**(-1)*s178**(-1)*
     &    s356**(-1) )
      bmp78(jdu,2,1,h56,1,1)= + gamz3456(jdu,2,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p7)*za(p1,p3)*za(p3,p5)*zb(p2,p1)**2*zb(p3,p6)*
     &    zab2(p8,p2,p7,p4)*izb(p1,p7)*izb(p1,p8)*s27**(-1)*s278**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p7)*za(p1,p3)*za(p4,p5)*zb(p2,p1)**2*
     &    zb(p4,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*izb(p1,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p1,p3)*za(p5,p6)*
     &    zb(p2,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p6)*izb(p1,p7)*izb(p1,p8
     &    )*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p1,p5)*
     &    za(p3,p5)*zb(p2,p1)**2*zb(p5,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)
     &    *izb(p1,p8)*s27**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p2,p7)*
     &    za(p1,p8)*za(p3,p5)*za(p3,p8)*zb(p2,p1)*zb(p2,p4)*zb(p3,p6)*
     &    izb(p1,p7)*s18**(-1)*s27**(-1)*s356**(-1) - 4.D0*za(p2,p7)*
     &    za(p1,p8)*za(p3,p5)*za(p5,p8)*zb(p2,p1)*zb(p2,p4)*zb(p5,p6)*
     &    izb(p1,p7)*s18**(-1)*s27**(-1)*s356**(-1) + 4.D0*za(p2,p7)*
     &    za(p1,p8)*za(p3,p8)*za(p4,p5)*zb(p2,p1)*zb(p2,p4)*zb(p4,p6)*
     &    izb(p1,p7)*s18**(-1)*s27**(-1)*s456**(-1) )
      bmp78(jdu,2,1,h56,1,1) = bmp78(jdu,2,1,h56,1,1) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p2,p7)*za(p1,p8)*za(p3,p8)*za(
     &    p5,p6)*zb(p2,p1)*zb(p2,p6)*zb(p4,p6)*izb(p1,p7)*s18**(-1)*
     &    s27**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p3,p5)*za(p7,p8)*zb(
     &    p2,p1)*zb(p3,p6)*zab2(p7,p2,p8,p4)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p3)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p3,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p3)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p4,p6)*zab2(p7,p2,p8,p4)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p3)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p1)*zb(p4,p6)*zab2(p7,p2,p8,p6)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p6)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) )
      bmp78(jdu,2,1,h56,1,1) = bmp78(jdu,2,1,h56,1,1) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * ( 4.D0*za(p1,p5)*za(p3,p5)*za(p7,p8)*zb(p2,
     &    p1)*zb(p5,p6)*zab2(p7,p2,p8,p4)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p5)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p5,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p7)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p3,p6)*zab2(p3,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p1,p7)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p1,p7)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p1,p7)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p6)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p1,p8)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p3,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) )
      bmp78(jdu,2,1,h56,1,1) = bmp78(jdu,2,1,h56,1,1) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p1,p8)*za(p3,p5)*za(p7,p8)*zb(
     &    p2,p4)*zb(p3,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p1,p8)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p1,p8)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p1,p8)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p6)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p1,p8)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p6)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s456**(-1) )
      bmp78(jdu,2,1,h56,2,1)= + gamz3456(jdu,2,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p8)*za(p1,p3)*za(p3,p5)*zb(p2,p7)**2*zb(p3,p6)*
     &    zab2(p8,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p8)*za(p1,p3)*za(p4,p5)*zb(p2,p7)**2*
     &    zb(p4,p6)*zab2(p8,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p2,p8)*za(p1,p3)*za(p5,p6)*
     &    zb(p2,p7)**2*zb(p4,p6)*zab2(p8,p2,p7,p6)*iza(p7,p8)*izb(p7,p8
     &    )*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p8)*za(p1,p5)*
     &    za(p3,p5)*zb(p2,p7)**2*zb(p5,p6)*zab2(p8,p2,p7,p4)*iza(p7,p8)
     &    *izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p1,p8)
     &    **2*za(p3,p5)*zb(p2,p4)*zb(p1,p7)*zb(p3,p6)*zab2(p3,p1,p8,p7)
     &    *iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) - 4.D0
     &    *za(p1,p8)**2*za(p3,p5)*zb(p2,p4)*zb(p1,p7)*zb(p5,p6)*zab2(p5
     &    ,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s356**(-1) + 4.D0*za(p1,p8)**2*za(p4,p5)*zb(p2,p4)*zb(p1,p7)*
     &    zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s456**(-1) )
      bmp78(jdu,2,1,h56,2,1) = bmp78(jdu,2,1,h56,2,1) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p1,p8)**2*za(p5,p6)*zb(p2,p6)*
     &    zb(p1,p7)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p1,p8)*za(p3,p5)*
     &    zb(p2,p7)*zb(p3,p6)*zab2(p3,p1,p8,p7)*zab2(p8,p2,p7,p4)*iza(
     &    p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s356**(-1) - 4.D0*za(p1
     &    ,p8)*za(p3,p5)*zb(p2,p7)*zb(p5,p6)*zab2(p5,p1,p8,p7)*zab2(p8,
     &    p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s356**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*zb(p2,p7)*zb(p4,p6)*
     &    zab2(p3,p1,p8,p7)*zab2(p8,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s27**(-1)*s456**(-1) - 4.D0*za(p1,p8)*za(p5,p6)*zb(
     &    p2,p7)*zb(p4,p6)*zab2(p3,p1,p8,p7)*zab2(p8,p2,p7,p6)*iza(p7,
     &    p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s456**(-1) )
      bmp78(jdu,2,1,h56,1,2)= + gamz3456(jdu,2,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p7)*za(p1,p3)*za(p1,p7)*za(p3,p5)*zb(p2,p4)*zb(
     &    p2,p8)*zb(p1,p8)*zb(p3,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s27**(-1)*s356**(-1) + 4.D0*za(p2,p7)*za(p1,p3)*za(p1,p7)*za(
     &    p4,p5)*zb(p2,p4)*zb(p2,p8)*zb(p1,p8)*zb(p4,p6)*iza(p7,p8)*
     &    izb(p7,p8)*s18**(-1)*s27**(-1)*s456**(-1) - 4.D0*za(p2,p7)*
     &    za(p1,p3)*za(p1,p7)*za(p5,p6)*zb(p2,p6)*zb(p2,p8)*zb(p1,p8)*
     &    zb(p4,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s456**(-1) - 4.D0*za(p2,p7)*za(p1,p3)*za(p3,p5)*zb(p2,p8)**2*
     &    zb(p3,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p2,p7)*za(p1,p3)*za(p4,p5)*
     &    zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p7,p8
     &    )*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p1,p3)*
     &    za(p5,p6)*zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,p6)*iza(p7,p8)
     &    *izb(p7,p8)*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*
     &    za(p1,p5)*za(p1,p7)*za(p3,p5)*zb(p2,p4)*zb(p2,p8)*zb(p1,p8)*
     &    zb(p5,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s356**(-1) )
      bmp78(jdu,2,1,h56,1,2) = bmp78(jdu,2,1,h56,1,2) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p2,p7)*za(p1,p5)*za(p3,p5)*zb(
     &    p2,p8)**2*zb(p5,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p7,p8)*
     &    s27**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p1,p7)**2*za(p3,p5)
     &    *zb(p2,p4)*zb(p1,p8)*zb(p3,p6)*zab2(p3,p1,p7,p8)*iza(p7,p8)*
     &    izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) - 4.D0*za(p1,p7)**
     &    2*za(p3,p5)*zb(p2,p4)*zb(p1,p8)*zb(p5,p6)*zab2(p5,p1,p7,p8)*
     &    iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*
     &    za(p1,p7)**2*za(p4,p5)*zb(p2,p4)*zb(p1,p8)*zb(p4,p6)*zab2(p3,
     &    p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s456**(-1) - 4.D0*za(p1,p7)**2*za(p5,p6)*zb(p2,p6)*zb(p1,p8)*
     &    zb(p4,p6)*zab2(p3,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s456**(-1) )
      bmp78(jdu,2,1,h56,2,2)= + gamz3456(jdu,2,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,p4)*zb(p1,p8)*zb(p3,p6)*
     &    zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*s18**(-1)*s178**(-1)*
     &    s356**(-1) - 4.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,p4)*zb(p1,p8)*
     &    zb(p5,p6)*zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*s18**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p1)**2*za(p4,p5)*zb(p2,p4)
     &    *zb(p1,p8)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*
     &    s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p1)**2*za(p5,p6)
     &    *zb(p2,p6)*zb(p1,p8)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p2,p7)*
     &    iza(p2,p8)*s18**(-1)*s178**(-1)*s456**(-1) + 4.D0*za(p2,p1)*
     &    za(p1,p3)*za(p3,p5)*zb(p2,p7)*zb(p1,p8)*zb(p3,p6)*zb(p4,p7)*
     &    iza(p2,p8)*s18**(-1)*s27**(-1)*s356**(-1) - 4.D0*za(p2,p1)*
     &    za(p1,p3)*za(p4,p5)*zb(p2,p7)*zb(p1,p8)*zb(p4,p6)*zb(p4,p7)*
     &    iza(p2,p8)*s18**(-1)*s27**(-1)*s456**(-1) + 4.D0*za(p2,p1)*
     &    za(p1,p3)*za(p5,p6)*zb(p2,p7)*zb(p1,p8)*zb(p4,p6)*zb(p6,p7)*
     &    iza(p2,p8)*s18**(-1)*s27**(-1)*s456**(-1) )
      bmp78(jdu,2,1,h56,2,2) = bmp78(jdu,2,1,h56,2,2) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * ( 4.D0*za(p2,p1)*za(p1,p5)*za(p3,p5)*zb(p2,
     &    p7)*zb(p1,p8)*zb(p4,p7)*zb(p5,p6)*iza(p2,p8)*s18**(-1)*
     &    s27**(-1)*s356**(-1) + 4.D0*za(p2,p1)*za(p3,p5)*zb(p2,p4)*zb(
     &    p3,p6)*zb(p7,p8)*zab2(p3,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p1)*za(p3,p5)*zb(p2,p4)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p3,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p1)*za(p3,p5)*zb(p2,p4)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p5,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p1)*za(p3,p5)*zb(p2,p4)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p5,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p2,p1)*za(p4,p5)*zb(p2,p4)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p3,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p2,p1)*za(p4,p5)*zb(p2,p4)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p3,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) )
      bmp78(jdu,2,1,h56,2,2) = bmp78(jdu,2,1,h56,2,2) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * ( 4.D0*za(p2,p1)*za(p5,p6)*zb(p2,p6)*zb(p4,
     &    p6)*zb(p7,p8)*zab2(p3,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p2,p1)*za(p5,p6)*zb(p2,p6)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p3,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p3,p5)*zb(p2,p7)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p3)*za(p3,p5)*zb(p2,p7)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p3)*za(p3,p5)*zb(p2,p8)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p3)*za(p4,p5)*zb(p2,p7)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p3)*za(p4,p5)*zb(p2,p7)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) )
      bmp78(jdu,2,1,h56,2,2) = bmp78(jdu,2,1,h56,2,2) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p1,p3)*za(p4,p5)*zb(p2,p8)*zb(
     &    p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p6)*zb(p2,p7)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p6)*zb(p2,p7)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p6)*zb(p2,p8)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p1,p5)*za(p3,p5)*zb(p2,p7)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p5)*za(p3,p5)*zb(p2,p7)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p5)*za(p3,p5)*zb(p2,p8)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) )

      bmp78(jdu,1,2,h56,1,1)= + gamz3456(jdu,1,2)*gamll56(2,h56) * ( 4.D
     &    0*za(p2,p7)*za(p3,p5)*za(p7,p8)*zb(p1,p3)*zb(p3,p6)*izb(p1,p8
     &    )*zba3(p1,p2,p7,p8,p4)*s78**(-1)*s278**(-1)*s356**(-1) + 4.D0
     &    *za(p2,p7)*za(p3,p5)*zb(p1,p3)*zb(p3,p6)*zba2(p1,p2,p7,p8)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p4)*s27**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p2,p7)*za(p4,p5)*za(p7,p8)*
     &    zb(p1,p3)*zb(p4,p6)*izb(p1,p8)*zba3(p1,p2,p7,p8,p4)*s78**(-1)
     &    *s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p4,p5)*za(p7,p8)*
     &    zb(p1,p3)*zb(p5,p6)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*s78**(-1)
     &    *s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p4,p5)*zb(p1,p3)*
     &    zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,
     &    p7,p8,p4)*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*
     &    za(p4,p5)*zb(p1,p3)*zb(p5,p6)*zba2(p1,p2,p7,p8)*izb(p1,p7)*
     &    izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*s27**(-1)*s278**(-1)*
     &    s456**(-1) - 4.D0*za(p2,p7)*za(p5,p6)*za(p7,p8)*zb(p1,p6)*zb(
     &    p3,p6)*izb(p1,p8)*zba3(p1,p2,p7,p8,p4)*s78**(-1)*s278**(-1)*
     &    s356**(-1) )
      bmp78(jdu,1,2,h56,1,1) = bmp78(jdu,1,2,h56,1,1) + gamz3456(jdu,1,
     & 2)*gamll56(2,h56) * (  - 4.D0*za(p2,p7)*za(p5,p6)*zb(p1,p6)*zb(
     &    p3,p6)*zba2(p1,p2,p7,p8)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,
     &    p8,p4)*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0*za(p2,p8)*za(p3
     &    ,p5)*za(p7,p8)*zb(p1,p3)*zb(p3,p6)*izb(p1,p7)*zba3(p1,p2,p7,
     &    p8,p4)*s78**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p2,p8)*za(p4
     &    ,p5)*za(p7,p8)*zb(p1,p3)*zb(p4,p6)*izb(p1,p7)*zba3(p1,p2,p7,
     &    p8,p4)*s78**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p8)*za(p4
     &    ,p5)*za(p7,p8)*zb(p1,p3)*zb(p5,p6)*izb(p1,p7)*zba3(p1,p2,p7,
     &    p8,p5)*s78**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p8)*za(p5
     &    ,p6)*za(p7,p8)*zb(p1,p6)*zb(p3,p6)*izb(p1,p7)*zba3(p1,p2,p7,
     &    p8,p4)*s78**(-1)*s278**(-1)*s356**(-1) )
      bmp78(jdu,1,2,h56,2,1)= + gamz3456(jdu,1,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p4)*za(p2,p8)*za(p1,p8)*za(p3,p5)*zb(p2,p7)*zb(
     &    p1,p3)*zb(p1,p7)*zb(p3,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s27**(-1)*s356**(-1) + 4.D0*za(p2,p4)*za(p2,p8)*za(p1,p8)*za(
     &    p4,p5)*zb(p2,p7)*zb(p1,p3)*zb(p1,p7)*zb(p4,p6)*iza(p7,p8)*
     &    izb(p7,p8)*s18**(-1)*s27**(-1)*s456**(-1) + 4.D0*za(p2,p4)*
     &    za(p2,p8)*za(p1,p8)*za(p5,p6)*zb(p2,p7)*zb(p1,p6)*zb(p1,p7)*
     &    zb(p3,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s356**(-1) - 4.D0*za(p2,p4)*za(p1,p8)*za(p3,p5)*zb(p1,p7)**2*
     &    zb(p3,p6)*zba2(p3,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p4)*za(p1,p8)*za(p4,p5)*
     &    zb(p1,p7)**2*zb(p4,p6)*zba2(p3,p1,p7,p8)*iza(p7,p8)*izb(p7,p8
     &    )*s18**(-1)*s178**(-1)*s456**(-1) + 4.D0*za(p2,p4)*za(p1,p8)*
     &    za(p5,p6)*zb(p1,p7)**2*zb(p3,p6)*zba2(p6,p1,p7,p8)*iza(p7,p8)
     &    *izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p5)*
     &    za(p2,p8)*za(p1,p8)*za(p4,p5)*zb(p2,p7)*zb(p1,p3)*zb(p1,p7)*
     &    zb(p5,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s456**(-1) )
      bmp78(jdu,1,2,h56,2,1) = bmp78(jdu,1,2,h56,2,1) + gamz3456(jdu,1,
     & 2)*gamll56(2,h56) * ( 4.D0*za(p2,p5)*za(p1,p8)*za(p4,p5)*zb(p1,
     &    p7)**2*zb(p5,p6)*zba2(p3,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p8)**2*za(p3,p5)
     &    *zb(p2,p7)*zb(p1,p3)*zb(p3,p6)*zba2(p7,p2,p8,p4)*iza(p7,p8)*
     &    izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0*za(p2,p8)**
     &    2*za(p4,p5)*zb(p2,p7)*zb(p1,p3)*zb(p4,p6)*zba2(p7,p2,p8,p4)*
     &    iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*s456**(-1) + 4.D0*
     &    za(p2,p8)**2*za(p4,p5)*zb(p2,p7)*zb(p1,p3)*zb(p5,p6)*zba2(p7,
     &    p2,p8,p5)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s456**(-1) + 4.D0*za(p2,p8)**2*za(p5,p6)*zb(p2,p7)*zb(p1,p6)*
     &    zb(p3,p6)*zba2(p7,p2,p8,p4)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) )
      bmp78(jdu,1,2,h56,1,2)= + gamz3456(jdu,1,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p4)*za(p1,p7)*za(p3,p5)*zb(p1,p8)**2*zb(p3,p6)*
     &    zba2(p3,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p4)*za(p1,p7)*za(p4,p5)*zb(p1,p8)**2*
     &    zb(p4,p6)*zba2(p3,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p2,p4)*za(p1,p7)*za(p5,p6)*
     &    zb(p1,p8)**2*zb(p3,p6)*zba2(p6,p1,p8,p7)*iza(p7,p8)*izb(p7,p8
     &    )*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p5)*za(p1,p7)*
     &    za(p4,p5)*zb(p1,p8)**2*zb(p5,p6)*zba2(p3,p1,p8,p7)*iza(p7,p8)
     &    *izb(p7,p8)*s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p7)
     &    **2*za(p3,p5)*zb(p2,p8)*zb(p1,p3)*zb(p3,p6)*zba2(p8,p2,p7,p4)
     &    *iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0
     &    *za(p2,p7)**2*za(p4,p5)*zb(p2,p8)*zb(p1,p3)*zb(p4,p6)*zba2(p8
     &    ,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s456**(-1) + 4.D0*za(p2,p7)**2*za(p4,p5)*zb(p2,p8)*zb(p1,p3)*
     &    zb(p5,p6)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) )
      bmp78(jdu,1,2,h56,1,2) = bmp78(jdu,1,2,h56,1,2) + gamz3456(jdu,1,
     & 2)*gamll56(2,h56) * ( 4.D0*za(p2,p7)**2*za(p5,p6)*zb(p2,p8)*zb(
     &    p1,p6)*zb(p3,p6)*zba2(p8,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*
     &    s27**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p2,p7)*za(p3,p5)*
     &    zb(p1,p8)*zb(p3,p6)*zba2(p3,p1,p8,p7)*zba2(p8,p2,p7,p4)*iza(
     &    p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s356**(-1) + 4.D0*za(p2
     &    ,p7)*za(p4,p5)*zb(p1,p8)*zb(p4,p6)*zba2(p3,p1,p8,p7)*zba2(p8,
     &    p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s456**(-1) + 4.D0*za(p2,p7)*za(p4,p5)*zb(p1,p8)*zb(p5,p6)*
     &    zba2(p3,p1,p8,p7)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s27**(-1)*s456**(-1) + 4.D0*za(p2,p7)*za(p5,p6)*zb(
     &    p1,p8)*zb(p3,p6)*zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p4)*iza(p7,
     &    p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s356**(-1) )
      bmp78(jdu,1,2,h56,2,2)= + gamz3456(jdu,1,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p4)*za(p3,p5)*zb(p1,p7)*zb(p3,p6)*zb(p7,p8)*
     &    iza(p2,p8)*zba3(p3,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*
     &    s356**(-1) - 4.D0*za(p2,p4)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*zb(
     &    p7,p8)*iza(p2,p7)*zba3(p3,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p4)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*
     &    zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*zba3(p3,p1,p7,p8,p2)*
     &    s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p4)*za(p4,p5)*
     &    zb(p1,p7)*zb(p4,p6)*zb(p7,p8)*iza(p2,p8)*zba3(p3,p1,p7,p8,p2)
     &    *s78**(-1)*s178**(-1)*s456**(-1) + 4.D0*za(p2,p4)*za(p4,p5)*
     &    zb(p1,p8)*zb(p4,p6)*zb(p7,p8)*iza(p2,p7)*zba3(p3,p1,p7,p8,p2)
     &    *s78**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p4)*za(p4,p5)*
     &    zb(p1,p8)*zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*
     &    zba3(p3,p1,p7,p8,p2)*s18**(-1)*s178**(-1)*s456**(-1) + 4.D0*
     &    za(p2,p4)*za(p5,p6)*zb(p1,p7)*zb(p3,p6)*zb(p7,p8)*iza(p2,p8)*
     &    zba3(p6,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*s356**(-1) )
      bmp78(jdu,1,2,h56,2,2) = bmp78(jdu,1,2,h56,2,2) + gamz3456(jdu,1,
     & 2)*gamll56(2,h56) * ( 4.D0*za(p2,p4)*za(p5,p6)*zb(p1,p8)*zb(p3,
     &    p6)*zb(p7,p8)*iza(p2,p7)*zba3(p6,p1,p7,p8,p2)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p2,p4)*za(p5,p6)*zb(p1,p8)*
     &    zb(p3,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*zba3(p6,p1,
     &    p7,p8,p2)*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p5)*
     &    za(p4,p5)*zb(p1,p7)*zb(p5,p6)*zb(p7,p8)*iza(p2,p8)*zba3(p3,p1
     &    ,p7,p8,p2)*s78**(-1)*s178**(-1)*s456**(-1) + 4.D0*za(p2,p5)*
     &    za(p4,p5)*zb(p1,p8)*zb(p5,p6)*zb(p7,p8)*iza(p2,p7)*zba3(p3,p1
     &    ,p7,p8,p2)*s78**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p5)*
     &    za(p4,p5)*zb(p1,p8)*zb(p5,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*
     &    iza(p2,p8)*zba3(p3,p1,p7,p8,p2)*s18**(-1)*s178**(-1)*
     &    s456**(-1) )
      bmp78(jdu,2,2,h56,1,1)= + gamz3456(jdu,2,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p7)*za(p1,p4)*za(p3,p5)*zb(p2,p1)**2*zb(p3,p6)*
     &    zab2(p8,p2,p7,p3)*izb(p1,p7)*izb(p1,p8)*s27**(-1)*s278**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p7)*za(p1,p4)*za(p4,p5)*zb(p2,p1)**2*
     &    zb(p4,p6)*zab2(p8,p2,p7,p3)*izb(p1,p7)*izb(p1,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p2,p7)*za(p1,p4)*za(p5,p6)*
     &    zb(p2,p1)**2*zb(p3,p6)*zab2(p8,p2,p7,p6)*izb(p1,p7)*izb(p1,p8
     &    )*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0*za(p2,p7)*za(p1,p5)*
     &    za(p4,p5)*zb(p2,p1)**2*zb(p5,p6)*zab2(p8,p2,p7,p3)*izb(p1,p7)
     &    *izb(p1,p8)*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*
     &    za(p1,p8)*za(p3,p5)*za(p4,p8)*zb(p2,p1)*zb(p2,p3)*zb(p3,p6)*
     &    izb(p1,p7)*s18**(-1)*s27**(-1)*s356**(-1) + 4.D0*za(p2,p7)*
     &    za(p1,p8)*za(p4,p5)*za(p4,p8)*zb(p2,p1)*zb(p2,p3)*zb(p4,p6)*
     &    izb(p1,p7)*s18**(-1)*s27**(-1)*s456**(-1) + 4.D0*za(p2,p7)*
     &    za(p1,p8)*za(p4,p5)*za(p5,p8)*zb(p2,p1)*zb(p2,p3)*zb(p5,p6)*
     &    izb(p1,p7)*s18**(-1)*s27**(-1)*s456**(-1) )
      bmp78(jdu,2,2,h56,1,1) = bmp78(jdu,2,2,h56,1,1) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * ( 4.D0*za(p2,p7)*za(p1,p8)*za(p4,p8)*za(p5,
     &    p6)*zb(p2,p1)*zb(p2,p6)*zb(p3,p6)*izb(p1,p7)*s18**(-1)*
     &    s27**(-1)*s356**(-1) + 4.D0*za(p1,p4)*za(p3,p5)*za(p7,p8)*zb(
     &    p2,p1)*zb(p3,p6)*zab2(p7,p2,p8,p3)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p4)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p3,p6)*zab2(p8,p2,p7,p3)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p4)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p4,p6)*zab2(p7,p2,p8,p3)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p4)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p3)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p4)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p1)*zb(p3,p6)*zab2(p7,p2,p8,p6)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p4)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p1)*zb(p3,p6)*zab2(p8,p2,p7,p6)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) )
      bmp78(jdu,2,2,h56,1,1) = bmp78(jdu,2,2,h56,1,1) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * (  - 4.D0*za(p1,p5)*za(p4,p5)*za(p7,p8)*zb(
     &    p2,p1)*zb(p5,p6)*zab2(p7,p2,p8,p3)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p5)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p5,p6)*zab2(p8,p2,p7,p3)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p7)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p3,p6)*zab2(p4,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p1,p7)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p4,p6)*zab2(p4,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p7)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p7)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p6)*zb(p3,p6)*zab2(p4,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p1,p8)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p3,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) )
      bmp78(jdu,2,2,h56,1,1) = bmp78(jdu,2,2,h56,1,1) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * (  - 4.D0*za(p1,p8)*za(p3,p5)*za(p7,p8)*zb(
     &    p2,p3)*zb(p3,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p4,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p4,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p8)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p6)*zb(p3,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p1,p8)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p6)*zb(p3,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s356**(-1) )
      bmp78(jdu,2,2,h56,2,1)= + gamz3456(jdu,2,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p8)*za(p1,p4)*za(p3,p5)*zb(p2,p7)**2*zb(p3,p6)*
     &    zab2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p8)*za(p1,p4)*za(p4,p5)*zb(p2,p7)**2*
     &    zb(p4,p6)*zab2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p2,p8)*za(p1,p4)*za(p5,p6)*
     &    zb(p2,p7)**2*zb(p3,p6)*zab2(p8,p2,p7,p6)*iza(p7,p8)*izb(p7,p8
     &    )*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0*za(p2,p8)*za(p1,p5)*
     &    za(p4,p5)*zb(p2,p7)**2*zb(p5,p6)*zab2(p8,p2,p7,p3)*iza(p7,p8)
     &    *izb(p7,p8)*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p1,p8)
     &    **2*za(p3,p5)*zb(p2,p3)*zb(p1,p7)*zb(p3,p6)*zab2(p4,p1,p8,p7)
     &    *iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0
     &    *za(p1,p8)**2*za(p4,p5)*zb(p2,p3)*zb(p1,p7)*zb(p4,p6)*zab2(p4
     &    ,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s456**(-1) + 4.D0*za(p1,p8)**2*za(p4,p5)*zb(p2,p3)*zb(p1,p7)*
     &    zb(p5,p6)*zab2(p5,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s456**(-1) )
      bmp78(jdu,2,2,h56,2,1) = bmp78(jdu,2,2,h56,2,1) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * ( 4.D0*za(p1,p8)**2*za(p5,p6)*zb(p2,p6)*zb(
     &    p1,p7)*zb(p3,p6)*zab2(p4,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s178**(-1)*s356**(-1) - 4.D0*za(p1,p8)*za(p3,p5)*
     &    zb(p2,p7)*zb(p3,p6)*zab2(p4,p1,p8,p7)*zab2(p8,p2,p7,p3)*iza(
     &    p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s356**(-1) + 4.D0*za(p1
     &    ,p8)*za(p4,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p4,p1,p8,p7)*zab2(p8,
     &    p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s456**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*zb(p2,p7)*zb(p5,p6)*
     &    zab2(p5,p1,p8,p7)*zab2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s27**(-1)*s456**(-1) + 4.D0*za(p1,p8)*za(p5,p6)*zb(
     &    p2,p7)*zb(p3,p6)*zab2(p4,p1,p8,p7)*zab2(p8,p2,p7,p6)*iza(p7,
     &    p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s356**(-1) )
      bmp78(jdu,2,2,h56,1,2)= + gamz3456(jdu,2,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p7)*za(p1,p4)*za(p1,p7)*za(p3,p5)*zb(p2,p3)*zb(
     &    p2,p8)*zb(p1,p8)*zb(p3,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s27**(-1)*s356**(-1) + 4.D0*za(p2,p7)*za(p1,p4)*za(p1,p7)*za(
     &    p4,p5)*zb(p2,p3)*zb(p2,p8)*zb(p1,p8)*zb(p4,p6)*iza(p7,p8)*
     &    izb(p7,p8)*s18**(-1)*s27**(-1)*s456**(-1) + 4.D0*za(p2,p7)*
     &    za(p1,p4)*za(p1,p7)*za(p5,p6)*zb(p2,p6)*zb(p2,p8)*zb(p1,p8)*
     &    zb(p3,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s356**(-1) - 4.D0*za(p2,p7)*za(p1,p4)*za(p3,p5)*zb(p2,p8)**2*
     &    zb(p3,p6)*zab2(p7,p2,p8,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p2,p7)*za(p1,p4)*za(p4,p5)*
     &    zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,p3)*iza(p7,p8)*izb(p7,p8
     &    )*s27**(-1)*s278**(-1)*s456**(-1) + 4.D0*za(p2,p7)*za(p1,p4)*
     &    za(p5,p6)*zb(p2,p8)**2*zb(p3,p6)*zab2(p7,p2,p8,p6)*iza(p7,p8)
     &    *izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0*za(p2,p7)*
     &    za(p1,p5)*za(p1,p7)*za(p4,p5)*zb(p2,p3)*zb(p2,p8)*zb(p1,p8)*
     &    zb(p5,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s456**(-1) )
      bmp78(jdu,2,2,h56,1,2) = bmp78(jdu,2,2,h56,1,2) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * ( 4.D0*za(p2,p7)*za(p1,p5)*za(p4,p5)*zb(p2,
     &    p8)**2*zb(p5,p6)*zab2(p7,p2,p8,p3)*iza(p7,p8)*izb(p7,p8)*
     &    s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p1,p7)**2*za(p3,p5)
     &    *zb(p2,p3)*zb(p1,p8)*zb(p3,p6)*zab2(p4,p1,p7,p8)*iza(p7,p8)*
     &    izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p1,p7)**
     &    2*za(p4,p5)*zb(p2,p3)*zb(p1,p8)*zb(p4,p6)*zab2(p4,p1,p7,p8)*
     &    iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*s456**(-1) + 4.D0*
     &    za(p1,p7)**2*za(p4,p5)*zb(p2,p3)*zb(p1,p8)*zb(p5,p6)*zab2(p5,
     &    p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s456**(-1) + 4.D0*za(p1,p7)**2*za(p5,p6)*zb(p2,p6)*zb(p1,p8)*
     &    zb(p3,p6)*zab2(p4,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s356**(-1) )
      bmp78(jdu,2,2,h56,2,2)= + gamz3456(jdu,2,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,p3)*zb(p1,p8)*zb(p3,p6)*
     &    zab2(p4,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*s18**(-1)*s178**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p1)**2*za(p4,p5)*zb(p2,p3)*zb(p1,p8)*
     &    zb(p4,p6)*zab2(p4,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*s18**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p2,p1)**2*za(p4,p5)*zb(p2,p3)
     &    *zb(p1,p8)*zb(p5,p6)*zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*
     &    s18**(-1)*s178**(-1)*s456**(-1) + 4.D0*za(p2,p1)**2*za(p5,p6)
     &    *zb(p2,p6)*zb(p1,p8)*zb(p3,p6)*zab2(p4,p1,p8,p7)*iza(p2,p7)*
     &    iza(p2,p8)*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p1)*
     &    za(p1,p4)*za(p3,p5)*zb(p2,p7)*zb(p1,p8)*zb(p3,p6)*zb(p3,p7)*
     &    iza(p2,p8)*s18**(-1)*s27**(-1)*s356**(-1) - 4.D0*za(p2,p1)*
     &    za(p1,p4)*za(p4,p5)*zb(p2,p7)*zb(p1,p8)*zb(p3,p7)*zb(p4,p6)*
     &    iza(p2,p8)*s18**(-1)*s27**(-1)*s456**(-1) - 4.D0*za(p2,p1)*
     &    za(p1,p4)*za(p5,p6)*zb(p2,p7)*zb(p1,p8)*zb(p3,p6)*zb(p6,p7)*
     &    iza(p2,p8)*s18**(-1)*s27**(-1)*s356**(-1) )
      bmp78(jdu,2,2,h56,2,2) = bmp78(jdu,2,2,h56,2,2) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * (  - 4.D0*za(p2,p1)*za(p1,p5)*za(p4,p5)*zb(
     &    p2,p7)*zb(p1,p8)*zb(p3,p7)*zb(p5,p6)*iza(p2,p8)*s18**(-1)*
     &    s27**(-1)*s456**(-1) + 4.D0*za(p2,p1)*za(p3,p5)*zb(p2,p3)*zb(
     &    p3,p6)*zb(p7,p8)*zab2(p4,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p1)*za(p3,p5)*zb(p2,p3)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p4,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p2,p1)*za(p4,p5)*zb(p2,p3)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p4,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p2,p1)*za(p4,p5)*zb(p2,p3)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p4,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p2,p1)*za(p4,p5)*zb(p2,p3)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p5,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p2,p1)*za(p4,p5)*zb(p2,p3)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p5,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) )
      bmp78(jdu,2,2,h56,2,2) = bmp78(jdu,2,2,h56,2,2) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * (  - 4.D0*za(p2,p1)*za(p5,p6)*zb(p2,p6)*zb(
     &    p3,p6)*zb(p7,p8)*zab2(p4,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p2,p1)*za(p5,p6)*zb(p2,p6)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p4,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p1,p4)*za(p3,p5)*zb(p2,p7)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p4)*za(p3,p5)*zb(p2,p7)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p4)*za(p3,p5)*zb(p2,p8)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p4)*za(p4,p5)*zb(p2,p7)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p4)*za(p4,p5)*zb(p2,p7)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) )
      bmp78(jdu,2,2,h56,2,2) = bmp78(jdu,2,2,h56,2,2) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * (  - 4.D0*za(p1,p4)*za(p4,p5)*zb(p2,p8)*zb(
     &    p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p4)*za(p5,p6)*zb(p2,p7)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p4)*za(p5,p6)*zb(p2,p7)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p4)*za(p5,p6)*zb(p2,p8)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p5)*za(p4,p5)*zb(p2,p7)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p5)*za(p4,p5)*zb(p2,p7)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p5)*za(p4,p5)*zb(p2,p8)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) )
      elseif (j78==2) then
      bmp87(jdu,1,1,h56,1,1)= + gamz3456(jdu,1,1)*gamll56(1,h56) * ( 4.D
     &    0*za(p2,p7)*za(p3,p5)*za(p7,p8)*zb(p1,p4)*zb(p3,p6)*izb(p1,p8
     &    )*zba3(p1,p2,p7,p8,p3)*s78**(-1)*s278**(-1)*s356**(-1) + 4.D0
     &    *za(p2,p7)*za(p3,p5)*za(p7,p8)*zb(p1,p4)*zb(p5,p6)*izb(p1,p8)
     &    *zba3(p1,p2,p7,p8,p5)*s78**(-1)*s278**(-1)*s356**(-1) + 4.D0*
     &    za(p2,p7)*za(p3,p5)*zb(p1,p4)*zb(p3,p6)*zba2(p1,p2,p7,p8)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*s27**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p2,p7)*za(p3,p5)*zb(p1,p4)*
     &    zb(p5,p6)*zba2(p1,p2,p7,p8)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,
     &    p7,p8,p5)*s27**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p2,p7)*
     &    za(p4,p5)*za(p7,p8)*zb(p1,p4)*zb(p4,p6)*izb(p1,p8)*zba3(p1,p2
     &    ,p7,p8,p3)*s78**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*
     &    za(p4,p5)*zb(p1,p4)*zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p1,p7)*
     &    izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*s27**(-1)*s278**(-1)*
     &    s456**(-1) + 4.D0*za(p2,p7)*za(p5,p6)*za(p7,p8)*zb(p1,p6)*zb(
     &    p4,p6)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*s78**(-1)*s278**(-1)*
     &    s456**(-1) )
      bmp87(jdu,1,1,h56,1,1) = bmp87(jdu,1,1,h56,1,1) + gamz3456(jdu,1,
     & 1)*gamll56(1,h56) * ( 4.D0*za(p2,p7)*za(p5,p6)*zb(p1,p6)*zb(p4,
     &    p6)*zba2(p1,p2,p7,p8)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,
     &    p3)*s27**(-1)*s278**(-1)*s456**(-1) + 4.D0*za(p2,p8)*za(p3,p5
     &    )*za(p7,p8)*zb(p1,p4)*zb(p3,p6)*izb(p1,p7)*zba3(p1,p2,p7,p8,
     &    p3)*s78**(-1)*s278**(-1)*s356**(-1) + 4.D0*za(p2,p8)*za(p3,p5
     &    )*za(p7,p8)*zb(p1,p4)*zb(p5,p6)*izb(p1,p7)*zba3(p1,p2,p7,p8,
     &    p5)*s78**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p2,p8)*za(p4,p5
     &    )*za(p7,p8)*zb(p1,p4)*zb(p4,p6)*izb(p1,p7)*zba3(p1,p2,p7,p8,
     &    p3)*s78**(-1)*s278**(-1)*s456**(-1) + 4.D0*za(p2,p8)*za(p5,p6
     &    )*za(p7,p8)*zb(p1,p6)*zb(p4,p6)*izb(p1,p7)*zba3(p1,p2,p7,p8,
     &    p3)*s78**(-1)*s278**(-1)*s456**(-1) )
      bmp87(jdu,1,1,h56,1,2)= + gamz3456(jdu,1,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p3)*za(p2,p8)*za(p1,p8)*za(p3,p5)*zb(p2,p7)*zb(
     &    p1,p4)*zb(p1,p7)*zb(p3,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s27**(-1)*s356**(-1) + 4.D0*za(p2,p3)*za(p2,p8)*za(p1,p8)*za(
     &    p4,p5)*zb(p2,p7)*zb(p1,p4)*zb(p1,p7)*zb(p4,p6)*iza(p7,p8)*
     &    izb(p7,p8)*s18**(-1)*s27**(-1)*s456**(-1) - 4.D0*za(p2,p3)*
     &    za(p2,p8)*za(p1,p8)*za(p5,p6)*zb(p2,p7)*zb(p1,p6)*zb(p1,p7)*
     &    zb(p4,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s456**(-1) - 4.D0*za(p2,p3)*za(p1,p8)*za(p3,p5)*zb(p1,p7)**2*
     &    zb(p3,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p3)*za(p1,p8)*za(p4,p5)*
     &    zb(p1,p7)**2*zb(p4,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p7,p8
     &    )*s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p3)*za(p1,p8)*
     &    za(p5,p6)*zb(p1,p7)**2*zb(p4,p6)*zba2(p6,p1,p7,p8)*iza(p7,p8)
     &    *izb(p7,p8)*s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p5)*
     &    za(p2,p8)*za(p1,p8)*za(p3,p5)*zb(p2,p7)*zb(p1,p4)*zb(p1,p7)*
     &    zb(p5,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s356**(-1) )
      bmp87(jdu,1,1,h56,1,2) = bmp87(jdu,1,1,h56,1,2) + gamz3456(jdu,1,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p2,p5)*za(p1,p8)*za(p3,p5)*zb(
     &    p1,p7)**2*zb(p5,p6)*zba2(p4,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s178**(-1)*s356**(-1) - 4.D0*za(p2,p8)**2*za(p3,p5)
     &    *zb(p2,p7)*zb(p1,p4)*zb(p3,p6)*zba2(p7,p2,p8,p3)*iza(p7,p8)*
     &    izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p2,p8)**
     &    2*za(p3,p5)*zb(p2,p7)*zb(p1,p4)*zb(p5,p6)*zba2(p7,p2,p8,p5)*
     &    iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0*
     &    za(p2,p8)**2*za(p4,p5)*zb(p2,p7)*zb(p1,p4)*zb(p4,p6)*zba2(p7,
     &    p2,p8,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s456**(-1) - 4.D0*za(p2,p8)**2*za(p5,p6)*zb(p2,p7)*zb(p1,p6)*
     &    zb(p4,p6)*zba2(p7,p2,p8,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) )
      bmp87(jdu,1,1,h56,2,1)= + gamz3456(jdu,1,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p3)*za(p1,p7)*za(p3,p5)*zb(p1,p8)**2*zb(p3,p6)*
     &    zba2(p4,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p3)*za(p1,p7)*za(p4,p5)*zb(p1,p8)**2*
     &    zb(p4,p6)*zba2(p4,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p2,p3)*za(p1,p7)*za(p5,p6)*
     &    zb(p1,p8)**2*zb(p4,p6)*zba2(p6,p1,p8,p7)*iza(p7,p8)*izb(p7,p8
     &    )*s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p5)*za(p1,p7)*
     &    za(p3,p5)*zb(p1,p8)**2*zb(p5,p6)*zba2(p4,p1,p8,p7)*iza(p7,p8)
     &    *izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) - 4.D0*za(p2,p7)
     &    **2*za(p3,p5)*zb(p2,p8)*zb(p1,p4)*zb(p3,p6)*zba2(p8,p2,p7,p3)
     &    *iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) - 4.D0
     &    *za(p2,p7)**2*za(p3,p5)*zb(p2,p8)*zb(p1,p4)*zb(p5,p6)*zba2(p8
     &    ,p2,p7,p5)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p7)**2*za(p4,p5)*zb(p2,p8)*zb(p1,p4)*
     &    zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) )
      bmp87(jdu,1,1,h56,2,1) = bmp87(jdu,1,1,h56,2,1) + gamz3456(jdu,1,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p2,p7)**2*za(p5,p6)*zb(p2,p8)*
     &    zb(p1,p6)*zb(p4,p6)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*
     &    s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p3,p5)*
     &    zb(p1,p8)*zb(p3,p6)*zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(
     &    p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s356**(-1) - 4.D0*za(p2
     &    ,p7)*za(p3,p5)*zb(p1,p8)*zb(p5,p6)*zba2(p4,p1,p8,p7)*zba2(p8,
     &    p2,p7,p5)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p7)*za(p4,p5)*zb(p1,p8)*zb(p4,p6)*
     &    zba2(p4,p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s27**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p5,p6)*zb(
     &    p1,p8)*zb(p4,p6)*zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p7,
     &    p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s456**(-1) )
      bmp87(jdu,1,1,h56,2,2)= + gamz3456(jdu,1,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p3)*za(p3,p5)*zb(p1,p7)*zb(p3,p6)*zb(p7,p8)*
     &    iza(p2,p8)*zba3(p4,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*
     &    s356**(-1) - 4.D0*za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*zb(
     &    p7,p8)*iza(p2,p7)*zba3(p4,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p3)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*
     &    zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*zba3(p4,p1,p7,p8,p2)*
     &    s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p3)*za(p4,p5)*
     &    zb(p1,p7)*zb(p4,p6)*zb(p7,p8)*iza(p2,p8)*zba3(p4,p1,p7,p8,p2)
     &    *s78**(-1)*s178**(-1)*s456**(-1) + 4.D0*za(p2,p3)*za(p4,p5)*
     &    zb(p1,p8)*zb(p4,p6)*zb(p7,p8)*iza(p2,p7)*zba3(p4,p1,p7,p8,p2)
     &    *s78**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p3)*za(p4,p5)*
     &    zb(p1,p8)*zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*
     &    zba3(p4,p1,p7,p8,p2)*s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*
     &    za(p2,p3)*za(p5,p6)*zb(p1,p7)*zb(p4,p6)*zb(p7,p8)*iza(p2,p8)*
     &    zba3(p6,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*s456**(-1) )
      bmp87(jdu,1,1,h56,2,2) = bmp87(jdu,1,1,h56,2,2) + gamz3456(jdu,1,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p2,p3)*za(p5,p6)*zb(p1,p8)*zb(
     &    p4,p6)*zb(p7,p8)*iza(p2,p7)*zba3(p6,p1,p7,p8,p2)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p2,p3)*za(p5,p6)*zb(p1,p8)*
     &    zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*zba3(p6,p1,
     &    p7,p8,p2)*s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p5)*
     &    za(p3,p5)*zb(p1,p7)*zb(p5,p6)*zb(p7,p8)*iza(p2,p8)*zba3(p4,p1
     &    ,p7,p8,p2)*s78**(-1)*s178**(-1)*s356**(-1) - 4.D0*za(p2,p5)*
     &    za(p3,p5)*zb(p1,p8)*zb(p5,p6)*zb(p7,p8)*iza(p2,p7)*zba3(p4,p1
     &    ,p7,p8,p2)*s78**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p5)*
     &    za(p3,p5)*zb(p1,p8)*zb(p5,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*
     &    iza(p2,p8)*zba3(p4,p1,p7,p8,p2)*s18**(-1)*s178**(-1)*
     &    s356**(-1) )
      bmp87(jdu,2,1,h56,1,1)= + gamz3456(jdu,2,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p7)*za(p1,p3)*za(p3,p5)*zb(p2,p1)**2*zb(p3,p6)*
     &    zab2(p8,p2,p7,p4)*izb(p1,p7)*izb(p1,p8)*s27**(-1)*s278**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p7)*za(p1,p3)*za(p4,p5)*zb(p2,p1)**2*
     &    zb(p4,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*izb(p1,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p1,p3)*za(p5,p6)*
     &    zb(p2,p1)**2*zb(p4,p6)*zab2(p8,p2,p7,p6)*izb(p1,p7)*izb(p1,p8
     &    )*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p1,p5)*
     &    za(p3,p5)*zb(p2,p1)**2*zb(p5,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)
     &    *izb(p1,p8)*s27**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p2,p7)*
     &    za(p1,p8)*za(p3,p5)*za(p3,p8)*zb(p2,p1)*zb(p2,p4)*zb(p3,p6)*
     &    izb(p1,p7)*s18**(-1)*s27**(-1)*s356**(-1) - 4.D0*za(p2,p7)*
     &    za(p1,p8)*za(p3,p5)*za(p5,p8)*zb(p2,p1)*zb(p2,p4)*zb(p5,p6)*
     &    izb(p1,p7)*s18**(-1)*s27**(-1)*s356**(-1) + 4.D0*za(p2,p7)*
     &    za(p1,p8)*za(p3,p8)*za(p4,p5)*zb(p2,p1)*zb(p2,p4)*zb(p4,p6)*
     &    izb(p1,p7)*s18**(-1)*s27**(-1)*s456**(-1) )
      bmp87(jdu,2,1,h56,1,1) = bmp87(jdu,2,1,h56,1,1) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p2,p7)*za(p1,p8)*za(p3,p8)*za(
     &    p5,p6)*zb(p2,p1)*zb(p2,p6)*zb(p4,p6)*izb(p1,p7)*s18**(-1)*
     &    s27**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p3,p5)*za(p7,p8)*zb(
     &    p2,p1)*zb(p3,p6)*zab2(p7,p2,p8,p4)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p3)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p3,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p3)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p4,p6)*zab2(p7,p2,p8,p4)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p3)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p1)*zb(p4,p6)*zab2(p7,p2,p8,p6)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p6)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) )
      bmp87(jdu,2,1,h56,1,1) = bmp87(jdu,2,1,h56,1,1) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * ( 4.D0*za(p1,p5)*za(p3,p5)*za(p7,p8)*zb(p2,
     &    p1)*zb(p5,p6)*zab2(p7,p2,p8,p4)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p5)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p5,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p7)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p3,p6)*zab2(p3,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p1,p7)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p1,p7)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p1,p7)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p6)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p1,p8)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p3,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) )
      bmp87(jdu,2,1,h56,1,1) = bmp87(jdu,2,1,h56,1,1) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p1,p8)*za(p3,p5)*za(p7,p8)*zb(
     &    p2,p4)*zb(p3,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p1,p8)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p1,p8)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p4)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p1,p8)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p6)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p1,p8)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p6)*zb(p4,p6)*zab2(p3,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s456**(-1) )
      bmp87(jdu,2,1,h56,1,2)= + gamz3456(jdu,2,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p8)*za(p1,p3)*za(p3,p5)*zb(p2,p7)**2*zb(p3,p6)*
     &    zab2(p8,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p8)*za(p1,p3)*za(p4,p5)*zb(p2,p7)**2*
     &    zb(p4,p6)*zab2(p8,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p2,p8)*za(p1,p3)*za(p5,p6)*
     &    zb(p2,p7)**2*zb(p4,p6)*zab2(p8,p2,p7,p6)*iza(p7,p8)*izb(p7,p8
     &    )*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p8)*za(p1,p5)*
     &    za(p3,p5)*zb(p2,p7)**2*zb(p5,p6)*zab2(p8,p2,p7,p4)*iza(p7,p8)
     &    *izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p1,p8)
     &    **2*za(p3,p5)*zb(p2,p4)*zb(p1,p7)*zb(p3,p6)*zab2(p3,p1,p8,p7)
     &    *iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) - 4.D0
     &    *za(p1,p8)**2*za(p3,p5)*zb(p2,p4)*zb(p1,p7)*zb(p5,p6)*zab2(p5
     &    ,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s356**(-1) + 4.D0*za(p1,p8)**2*za(p4,p5)*zb(p2,p4)*zb(p1,p7)*
     &    zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s456**(-1) )
      bmp87(jdu,2,1,h56,1,2) = bmp87(jdu,2,1,h56,1,2) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p1,p8)**2*za(p5,p6)*zb(p2,p6)*
     &    zb(p1,p7)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p1,p8)*za(p3,p5)*
     &    zb(p2,p7)*zb(p3,p6)*zab2(p3,p1,p8,p7)*zab2(p8,p2,p7,p4)*iza(
     &    p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s356**(-1) - 4.D0*za(p1
     &    ,p8)*za(p3,p5)*zb(p2,p7)*zb(p5,p6)*zab2(p5,p1,p8,p7)*zab2(p8,
     &    p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s356**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*zb(p2,p7)*zb(p4,p6)*
     &    zab2(p3,p1,p8,p7)*zab2(p8,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s27**(-1)*s456**(-1) - 4.D0*za(p1,p8)*za(p5,p6)*zb(
     &    p2,p7)*zb(p4,p6)*zab2(p3,p1,p8,p7)*zab2(p8,p2,p7,p6)*iza(p7,
     &    p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s456**(-1) )
      bmp87(jdu,2,1,h56,2,1)= + gamz3456(jdu,2,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p7)*za(p1,p3)*za(p1,p7)*za(p3,p5)*zb(p2,p4)*zb(
     &    p2,p8)*zb(p1,p8)*zb(p3,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s27**(-1)*s356**(-1) + 4.D0*za(p2,p7)*za(p1,p3)*za(p1,p7)*za(
     &    p4,p5)*zb(p2,p4)*zb(p2,p8)*zb(p1,p8)*zb(p4,p6)*iza(p7,p8)*
     &    izb(p7,p8)*s18**(-1)*s27**(-1)*s456**(-1) - 4.D0*za(p2,p7)*
     &    za(p1,p3)*za(p1,p7)*za(p5,p6)*zb(p2,p6)*zb(p2,p8)*zb(p1,p8)*
     &    zb(p4,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s456**(-1) - 4.D0*za(p2,p7)*za(p1,p3)*za(p3,p5)*zb(p2,p8)**2*
     &    zb(p3,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p2,p7)*za(p1,p3)*za(p4,p5)*
     &    zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p7,p8
     &    )*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p1,p3)*
     &    za(p5,p6)*zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,p6)*iza(p7,p8)
     &    *izb(p7,p8)*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*
     &    za(p1,p5)*za(p1,p7)*za(p3,p5)*zb(p2,p4)*zb(p2,p8)*zb(p1,p8)*
     &    zb(p5,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s356**(-1) )
      bmp87(jdu,2,1,h56,2,1) = bmp87(jdu,2,1,h56,2,1) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p2,p7)*za(p1,p5)*za(p3,p5)*zb(
     &    p2,p8)**2*zb(p5,p6)*zab2(p7,p2,p8,p4)*iza(p7,p8)*izb(p7,p8)*
     &    s27**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p1,p7)**2*za(p3,p5)
     &    *zb(p2,p4)*zb(p1,p8)*zb(p3,p6)*zab2(p3,p1,p7,p8)*iza(p7,p8)*
     &    izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) - 4.D0*za(p1,p7)**
     &    2*za(p3,p5)*zb(p2,p4)*zb(p1,p8)*zb(p5,p6)*zab2(p5,p1,p7,p8)*
     &    iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*
     &    za(p1,p7)**2*za(p4,p5)*zb(p2,p4)*zb(p1,p8)*zb(p4,p6)*zab2(p3,
     &    p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s456**(-1) - 4.D0*za(p1,p7)**2*za(p5,p6)*zb(p2,p6)*zb(p1,p8)*
     &    zb(p4,p6)*zab2(p3,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s456**(-1) )
      bmp87(jdu,2,1,h56,2,2)= + gamz3456(jdu,2,1)*gamll56(1,h56) * ( 
     &     - 4.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,p4)*zb(p1,p8)*zb(p3,p6)*
     &    zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*s18**(-1)*s178**(-1)*
     &    s356**(-1) - 4.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,p4)*zb(p1,p8)*
     &    zb(p5,p6)*zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*s18**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p1)**2*za(p4,p5)*zb(p2,p4)
     &    *zb(p1,p8)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*
     &    s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p1)**2*za(p5,p6)
     &    *zb(p2,p6)*zb(p1,p8)*zb(p4,p6)*zab2(p3,p1,p8,p7)*iza(p2,p7)*
     &    iza(p2,p8)*s18**(-1)*s178**(-1)*s456**(-1) + 4.D0*za(p2,p1)*
     &    za(p1,p3)*za(p3,p5)*zb(p2,p7)*zb(p1,p8)*zb(p3,p6)*zb(p4,p7)*
     &    iza(p2,p8)*s18**(-1)*s27**(-1)*s356**(-1) - 4.D0*za(p2,p1)*
     &    za(p1,p3)*za(p4,p5)*zb(p2,p7)*zb(p1,p8)*zb(p4,p6)*zb(p4,p7)*
     &    iza(p2,p8)*s18**(-1)*s27**(-1)*s456**(-1) + 4.D0*za(p2,p1)*
     &    za(p1,p3)*za(p5,p6)*zb(p2,p7)*zb(p1,p8)*zb(p4,p6)*zb(p6,p7)*
     &    iza(p2,p8)*s18**(-1)*s27**(-1)*s456**(-1) )
      bmp87(jdu,2,1,h56,2,2) = bmp87(jdu,2,1,h56,2,2) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * ( 4.D0*za(p2,p1)*za(p1,p5)*za(p3,p5)*zb(p2,
     &    p7)*zb(p1,p8)*zb(p4,p7)*zb(p5,p6)*iza(p2,p8)*s18**(-1)*
     &    s27**(-1)*s356**(-1) + 4.D0*za(p2,p1)*za(p3,p5)*zb(p2,p4)*zb(
     &    p3,p6)*zb(p7,p8)*zab2(p3,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p1)*za(p3,p5)*zb(p2,p4)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p3,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p1)*za(p3,p5)*zb(p2,p4)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p5,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p1)*za(p3,p5)*zb(p2,p4)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p5,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p2,p1)*za(p4,p5)*zb(p2,p4)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p3,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p2,p1)*za(p4,p5)*zb(p2,p4)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p3,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) )
      bmp87(jdu,2,1,h56,2,2) = bmp87(jdu,2,1,h56,2,2) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * ( 4.D0*za(p2,p1)*za(p5,p6)*zb(p2,p6)*zb(p4,
     &    p6)*zb(p7,p8)*zab2(p3,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p2,p1)*za(p5,p6)*zb(p2,p6)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p3,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p3,p5)*zb(p2,p7)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p3)*za(p3,p5)*zb(p2,p7)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p3)*za(p3,p5)*zb(p2,p8)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p3)*za(p4,p5)*zb(p2,p7)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p3)*za(p4,p5)*zb(p2,p7)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) )
      bmp87(jdu,2,1,h56,2,2) = bmp87(jdu,2,1,h56,2,2) + gamz3456(jdu,2,
     & 1)*gamll56(1,h56) * (  - 4.D0*za(p1,p3)*za(p4,p5)*zb(p2,p8)*zb(
     &    p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p6)*zb(p2,p7)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p6)*zb(p2,p7)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p1,p3)*za(p5,p6)*zb(p2,p8)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p1,p5)*za(p3,p5)*zb(p2,p7)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p5)*za(p3,p5)*zb(p2,p7)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p5)*za(p3,p5)*zb(p2,p8)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p4)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) )

      bmp87(jdu,1,2,h56,1,1)= + gamz3456(jdu,1,2)*gamll56(2,h56) * ( 4.D
     &    0*za(p2,p7)*za(p3,p5)*za(p7,p8)*zb(p1,p3)*zb(p3,p6)*izb(p1,p8
     &    )*zba3(p1,p2,p7,p8,p4)*s78**(-1)*s278**(-1)*s356**(-1) + 4.D0
     &    *za(p2,p7)*za(p3,p5)*zb(p1,p3)*zb(p3,p6)*zba2(p1,p2,p7,p8)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,p8,p4)*s27**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p2,p7)*za(p4,p5)*za(p7,p8)*
     &    zb(p1,p3)*zb(p4,p6)*izb(p1,p8)*zba3(p1,p2,p7,p8,p4)*s78**(-1)
     &    *s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p4,p5)*za(p7,p8)*
     &    zb(p1,p3)*zb(p5,p6)*izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*s78**(-1)
     &    *s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*za(p4,p5)*zb(p1,p3)*
     &    zb(p4,p6)*zba2(p1,p2,p7,p8)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,
     &    p7,p8,p4)*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*
     &    za(p4,p5)*zb(p1,p3)*zb(p5,p6)*zba2(p1,p2,p7,p8)*izb(p1,p7)*
     &    izb(p1,p8)*zba3(p1,p2,p7,p8,p5)*s27**(-1)*s278**(-1)*
     &    s456**(-1) - 4.D0*za(p2,p7)*za(p5,p6)*za(p7,p8)*zb(p1,p6)*zb(
     &    p3,p6)*izb(p1,p8)*zba3(p1,p2,p7,p8,p4)*s78**(-1)*s278**(-1)*
     &    s356**(-1) )
      bmp87(jdu,1,2,h56,1,1) = bmp87(jdu,1,2,h56,1,1) + gamz3456(jdu,1,
     & 2)*gamll56(2,h56) * (  - 4.D0*za(p2,p7)*za(p5,p6)*zb(p1,p6)*zb(
     &    p3,p6)*zba2(p1,p2,p7,p8)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,
     &    p8,p4)*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0*za(p2,p8)*za(p3
     &    ,p5)*za(p7,p8)*zb(p1,p3)*zb(p3,p6)*izb(p1,p7)*zba3(p1,p2,p7,
     &    p8,p4)*s78**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p2,p8)*za(p4
     &    ,p5)*za(p7,p8)*zb(p1,p3)*zb(p4,p6)*izb(p1,p7)*zba3(p1,p2,p7,
     &    p8,p4)*s78**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p8)*za(p4
     &    ,p5)*za(p7,p8)*zb(p1,p3)*zb(p5,p6)*izb(p1,p7)*zba3(p1,p2,p7,
     &    p8,p5)*s78**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p8)*za(p5
     &    ,p6)*za(p7,p8)*zb(p1,p6)*zb(p3,p6)*izb(p1,p7)*zba3(p1,p2,p7,
     &    p8,p4)*s78**(-1)*s278**(-1)*s356**(-1) )
      bmp87(jdu,1,2,h56,1,2)= + gamz3456(jdu,1,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p4)*za(p2,p8)*za(p1,p8)*za(p3,p5)*zb(p2,p7)*zb(
     &    p1,p3)*zb(p1,p7)*zb(p3,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s27**(-1)*s356**(-1) + 4.D0*za(p2,p4)*za(p2,p8)*za(p1,p8)*za(
     &    p4,p5)*zb(p2,p7)*zb(p1,p3)*zb(p1,p7)*zb(p4,p6)*iza(p7,p8)*
     &    izb(p7,p8)*s18**(-1)*s27**(-1)*s456**(-1) + 4.D0*za(p2,p4)*
     &    za(p2,p8)*za(p1,p8)*za(p5,p6)*zb(p2,p7)*zb(p1,p6)*zb(p1,p7)*
     &    zb(p3,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s356**(-1) - 4.D0*za(p2,p4)*za(p1,p8)*za(p3,p5)*zb(p1,p7)**2*
     &    zb(p3,p6)*zba2(p3,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p4)*za(p1,p8)*za(p4,p5)*
     &    zb(p1,p7)**2*zb(p4,p6)*zba2(p3,p1,p7,p8)*iza(p7,p8)*izb(p7,p8
     &    )*s18**(-1)*s178**(-1)*s456**(-1) + 4.D0*za(p2,p4)*za(p1,p8)*
     &    za(p5,p6)*zb(p1,p7)**2*zb(p3,p6)*zba2(p6,p1,p7,p8)*iza(p7,p8)
     &    *izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p5)*
     &    za(p2,p8)*za(p1,p8)*za(p4,p5)*zb(p2,p7)*zb(p1,p3)*zb(p1,p7)*
     &    zb(p5,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s456**(-1) )
      bmp87(jdu,1,2,h56,1,2) = bmp87(jdu,1,2,h56,1,2) + gamz3456(jdu,1,
     & 2)*gamll56(2,h56) * ( 4.D0*za(p2,p5)*za(p1,p8)*za(p4,p5)*zb(p1,
     &    p7)**2*zb(p5,p6)*zba2(p3,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p8)**2*za(p3,p5)
     &    *zb(p2,p7)*zb(p1,p3)*zb(p3,p6)*zba2(p7,p2,p8,p4)*iza(p7,p8)*
     &    izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0*za(p2,p8)**
     &    2*za(p4,p5)*zb(p2,p7)*zb(p1,p3)*zb(p4,p6)*zba2(p7,p2,p8,p4)*
     &    iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*s456**(-1) + 4.D0*
     &    za(p2,p8)**2*za(p4,p5)*zb(p2,p7)*zb(p1,p3)*zb(p5,p6)*zba2(p7,
     &    p2,p8,p5)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s456**(-1) + 4.D0*za(p2,p8)**2*za(p5,p6)*zb(p2,p7)*zb(p1,p6)*
     &    zb(p3,p6)*zba2(p7,p2,p8,p4)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) )
      bmp87(jdu,1,2,h56,2,1)= + gamz3456(jdu,1,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p4)*za(p1,p7)*za(p3,p5)*zb(p1,p8)**2*zb(p3,p6)*
     &    zba2(p3,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p4)*za(p1,p7)*za(p4,p5)*zb(p1,p8)**2*
     &    zb(p4,p6)*zba2(p3,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p2,p4)*za(p1,p7)*za(p5,p6)*
     &    zb(p1,p8)**2*zb(p3,p6)*zba2(p6,p1,p8,p7)*iza(p7,p8)*izb(p7,p8
     &    )*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p5)*za(p1,p7)*
     &    za(p4,p5)*zb(p1,p8)**2*zb(p5,p6)*zba2(p3,p1,p8,p7)*iza(p7,p8)
     &    *izb(p7,p8)*s18**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p7)
     &    **2*za(p3,p5)*zb(p2,p8)*zb(p1,p3)*zb(p3,p6)*zba2(p8,p2,p7,p4)
     &    *iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0
     &    *za(p2,p7)**2*za(p4,p5)*zb(p2,p8)*zb(p1,p3)*zb(p4,p6)*zba2(p8
     &    ,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s456**(-1) + 4.D0*za(p2,p7)**2*za(p4,p5)*zb(p2,p8)*zb(p1,p3)*
     &    zb(p5,p6)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) )
      bmp87(jdu,1,2,h56,2,1) = bmp87(jdu,1,2,h56,2,1) + gamz3456(jdu,1,
     & 2)*gamll56(2,h56) * ( 4.D0*za(p2,p7)**2*za(p5,p6)*zb(p2,p8)*zb(
     &    p1,p6)*zb(p3,p6)*zba2(p8,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*
     &    s27**(-1)*s278**(-1)*s356**(-1) - 4.D0*za(p2,p7)*za(p3,p5)*
     &    zb(p1,p8)*zb(p3,p6)*zba2(p3,p1,p8,p7)*zba2(p8,p2,p7,p4)*iza(
     &    p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s356**(-1) + 4.D0*za(p2
     &    ,p7)*za(p4,p5)*zb(p1,p8)*zb(p4,p6)*zba2(p3,p1,p8,p7)*zba2(p8,
     &    p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s456**(-1) + 4.D0*za(p2,p7)*za(p4,p5)*zb(p1,p8)*zb(p5,p6)*
     &    zba2(p3,p1,p8,p7)*zba2(p8,p2,p7,p5)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s27**(-1)*s456**(-1) + 4.D0*za(p2,p7)*za(p5,p6)*zb(
     &    p1,p8)*zb(p3,p6)*zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p4)*iza(p7,
     &    p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s356**(-1) )
      bmp87(jdu,1,2,h56,2,2)= + gamz3456(jdu,1,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p4)*za(p3,p5)*zb(p1,p7)*zb(p3,p6)*zb(p7,p8)*
     &    iza(p2,p8)*zba3(p3,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*
     &    s356**(-1) - 4.D0*za(p2,p4)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*zb(
     &    p7,p8)*iza(p2,p7)*zba3(p3,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p4)*za(p3,p5)*zb(p1,p8)*zb(p3,p6)*
     &    zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*zba3(p3,p1,p7,p8,p2)*
     &    s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p4)*za(p4,p5)*
     &    zb(p1,p7)*zb(p4,p6)*zb(p7,p8)*iza(p2,p8)*zba3(p3,p1,p7,p8,p2)
     &    *s78**(-1)*s178**(-1)*s456**(-1) + 4.D0*za(p2,p4)*za(p4,p5)*
     &    zb(p1,p8)*zb(p4,p6)*zb(p7,p8)*iza(p2,p7)*zba3(p3,p1,p7,p8,p2)
     &    *s78**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p4)*za(p4,p5)*
     &    zb(p1,p8)*zb(p4,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*
     &    zba3(p3,p1,p7,p8,p2)*s18**(-1)*s178**(-1)*s456**(-1) + 4.D0*
     &    za(p2,p4)*za(p5,p6)*zb(p1,p7)*zb(p3,p6)*zb(p7,p8)*iza(p2,p8)*
     &    zba3(p6,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*s356**(-1) )
      bmp87(jdu,1,2,h56,2,2) = bmp87(jdu,1,2,h56,2,2) + gamz3456(jdu,1,
     & 2)*gamll56(2,h56) * ( 4.D0*za(p2,p4)*za(p5,p6)*zb(p1,p8)*zb(p3,
     &    p6)*zb(p7,p8)*iza(p2,p7)*zba3(p6,p1,p7,p8,p2)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p2,p4)*za(p5,p6)*zb(p1,p8)*
     &    zb(p3,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*zba3(p6,p1,
     &    p7,p8,p2)*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p5)*
     &    za(p4,p5)*zb(p1,p7)*zb(p5,p6)*zb(p7,p8)*iza(p2,p8)*zba3(p3,p1
     &    ,p7,p8,p2)*s78**(-1)*s178**(-1)*s456**(-1) + 4.D0*za(p2,p5)*
     &    za(p4,p5)*zb(p1,p8)*zb(p5,p6)*zb(p7,p8)*iza(p2,p7)*zba3(p3,p1
     &    ,p7,p8,p2)*s78**(-1)*s178**(-1)*s456**(-1) - 4.D0*za(p2,p5)*
     &    za(p4,p5)*zb(p1,p8)*zb(p5,p6)*zba2(p7,p1,p8,p2)*iza(p2,p7)*
     &    iza(p2,p8)*zba3(p3,p1,p7,p8,p2)*s18**(-1)*s178**(-1)*
     &    s456**(-1) )
      bmp87(jdu,2,2,h56,1,1)= + gamz3456(jdu,2,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p7)*za(p1,p4)*za(p3,p5)*zb(p2,p1)**2*zb(p3,p6)*
     &    zab2(p8,p2,p7,p3)*izb(p1,p7)*izb(p1,p8)*s27**(-1)*s278**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p7)*za(p1,p4)*za(p4,p5)*zb(p2,p1)**2*
     &    zb(p4,p6)*zab2(p8,p2,p7,p3)*izb(p1,p7)*izb(p1,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p2,p7)*za(p1,p4)*za(p5,p6)*
     &    zb(p2,p1)**2*zb(p3,p6)*zab2(p8,p2,p7,p6)*izb(p1,p7)*izb(p1,p8
     &    )*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0*za(p2,p7)*za(p1,p5)*
     &    za(p4,p5)*zb(p2,p1)**2*zb(p5,p6)*zab2(p8,p2,p7,p3)*izb(p1,p7)
     &    *izb(p1,p8)*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p2,p7)*
     &    za(p1,p8)*za(p3,p5)*za(p4,p8)*zb(p2,p1)*zb(p2,p3)*zb(p3,p6)*
     &    izb(p1,p7)*s18**(-1)*s27**(-1)*s356**(-1) + 4.D0*za(p2,p7)*
     &    za(p1,p8)*za(p4,p5)*za(p4,p8)*zb(p2,p1)*zb(p2,p3)*zb(p4,p6)*
     &    izb(p1,p7)*s18**(-1)*s27**(-1)*s456**(-1) + 4.D0*za(p2,p7)*
     &    za(p1,p8)*za(p4,p5)*za(p5,p8)*zb(p2,p1)*zb(p2,p3)*zb(p5,p6)*
     &    izb(p1,p7)*s18**(-1)*s27**(-1)*s456**(-1) )
      bmp87(jdu,2,2,h56,1,1) = bmp87(jdu,2,2,h56,1,1) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * ( 4.D0*za(p2,p7)*za(p1,p8)*za(p4,p8)*za(p5,
     &    p6)*zb(p2,p1)*zb(p2,p6)*zb(p3,p6)*izb(p1,p7)*s18**(-1)*
     &    s27**(-1)*s356**(-1) + 4.D0*za(p1,p4)*za(p3,p5)*za(p7,p8)*zb(
     &    p2,p1)*zb(p3,p6)*zab2(p7,p2,p8,p3)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p4)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p3,p6)*zab2(p8,p2,p7,p3)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p4)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p4,p6)*zab2(p7,p2,p8,p3)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p4)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p4,p6)*zab2(p8,p2,p7,p3)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p4)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p1)*zb(p3,p6)*zab2(p7,p2,p8,p6)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p4)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p1)*zb(p3,p6)*zab2(p8,p2,p7,p6)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) )
      bmp87(jdu,2,2,h56,1,1) = bmp87(jdu,2,2,h56,1,1) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * (  - 4.D0*za(p1,p5)*za(p4,p5)*za(p7,p8)*zb(
     &    p2,p1)*zb(p5,p6)*zab2(p7,p2,p8,p3)*izb(p1,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p5)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p1)*zb(p5,p6)*zab2(p8,p2,p7,p3)*izb(p1,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p7)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p3,p6)*zab2(p4,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p1,p7)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p4,p6)*zab2(p4,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p7)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p7)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p6)*zb(p3,p6)*zab2(p4,p7,p8,p1)*izb(p1,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p1,p8)*za(p3,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p3,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) )
      bmp87(jdu,2,2,h56,1,1) = bmp87(jdu,2,2,h56,1,1) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * (  - 4.D0*za(p1,p8)*za(p3,p5)*za(p7,p8)*zb(
     &    p2,p3)*zb(p3,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p4,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p4,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*za(p7,p8)*
     &    zb(p2,p3)*zb(p5,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p1,p8)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p6)*zb(p3,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p1,p8)*za(p5,p6)*za(p7,p8)*
     &    zb(p2,p6)*zb(p3,p6)*zab2(p4,p7,p8,p1)*izb(p1,p7)*s18**(-1)*
     &    s178**(-1)*s356**(-1) )
      bmp87(jdu,2,2,h56,1,2)= + gamz3456(jdu,2,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p8)*za(p1,p4)*za(p3,p5)*zb(p2,p7)**2*zb(p3,p6)*
     &    zab2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p8)*za(p1,p4)*za(p4,p5)*zb(p2,p7)**2*
     &    zb(p4,p6)*zab2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) + 4.D0*za(p2,p8)*za(p1,p4)*za(p5,p6)*
     &    zb(p2,p7)**2*zb(p3,p6)*zab2(p8,p2,p7,p6)*iza(p7,p8)*izb(p7,p8
     &    )*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0*za(p2,p8)*za(p1,p5)*
     &    za(p4,p5)*zb(p2,p7)**2*zb(p5,p6)*zab2(p8,p2,p7,p3)*iza(p7,p8)
     &    *izb(p7,p8)*s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p1,p8)
     &    **2*za(p3,p5)*zb(p2,p3)*zb(p1,p7)*zb(p3,p6)*zab2(p4,p1,p8,p7)
     &    *iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0
     &    *za(p1,p8)**2*za(p4,p5)*zb(p2,p3)*zb(p1,p7)*zb(p4,p6)*zab2(p4
     &    ,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s456**(-1) + 4.D0*za(p1,p8)**2*za(p4,p5)*zb(p2,p3)*zb(p1,p7)*
     &    zb(p5,p6)*zab2(p5,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s456**(-1) )
      bmp87(jdu,2,2,h56,1,2) = bmp87(jdu,2,2,h56,1,2) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * ( 4.D0*za(p1,p8)**2*za(p5,p6)*zb(p2,p6)*zb(
     &    p1,p7)*zb(p3,p6)*zab2(p4,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s178**(-1)*s356**(-1) - 4.D0*za(p1,p8)*za(p3,p5)*
     &    zb(p2,p7)*zb(p3,p6)*zab2(p4,p1,p8,p7)*zab2(p8,p2,p7,p3)*iza(
     &    p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s356**(-1) + 4.D0*za(p1
     &    ,p8)*za(p4,p5)*zb(p2,p7)*zb(p4,p6)*zab2(p4,p1,p8,p7)*zab2(p8,
     &    p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s456**(-1) + 4.D0*za(p1,p8)*za(p4,p5)*zb(p2,p7)*zb(p5,p6)*
     &    zab2(p5,p1,p8,p7)*zab2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*
     &    s18**(-1)*s27**(-1)*s456**(-1) + 4.D0*za(p1,p8)*za(p5,p6)*zb(
     &    p2,p7)*zb(p3,p6)*zab2(p4,p1,p8,p7)*zab2(p8,p2,p7,p6)*iza(p7,
     &    p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*s356**(-1) )
      bmp87(jdu,2,2,h56,2,1)= + gamz3456(jdu,2,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p7)*za(p1,p4)*za(p1,p7)*za(p3,p5)*zb(p2,p3)*zb(
     &    p2,p8)*zb(p1,p8)*zb(p3,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s27**(-1)*s356**(-1) + 4.D0*za(p2,p7)*za(p1,p4)*za(p1,p7)*za(
     &    p4,p5)*zb(p2,p3)*zb(p2,p8)*zb(p1,p8)*zb(p4,p6)*iza(p7,p8)*
     &    izb(p7,p8)*s18**(-1)*s27**(-1)*s456**(-1) + 4.D0*za(p2,p7)*
     &    za(p1,p4)*za(p1,p7)*za(p5,p6)*zb(p2,p6)*zb(p2,p8)*zb(p1,p8)*
     &    zb(p3,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s356**(-1) - 4.D0*za(p2,p7)*za(p1,p4)*za(p3,p5)*zb(p2,p8)**2*
     &    zb(p3,p6)*zab2(p7,p2,p8,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p2,p7)*za(p1,p4)*za(p4,p5)*
     &    zb(p2,p8)**2*zb(p4,p6)*zab2(p7,p2,p8,p3)*iza(p7,p8)*izb(p7,p8
     &    )*s27**(-1)*s278**(-1)*s456**(-1) + 4.D0*za(p2,p7)*za(p1,p4)*
     &    za(p5,p6)*zb(p2,p8)**2*zb(p3,p6)*zab2(p7,p2,p8,p6)*iza(p7,p8)
     &    *izb(p7,p8)*s27**(-1)*s278**(-1)*s356**(-1) + 4.D0*za(p2,p7)*
     &    za(p1,p5)*za(p1,p7)*za(p4,p5)*zb(p2,p3)*zb(p2,p8)*zb(p1,p8)*
     &    zb(p5,p6)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s27**(-1)*
     &    s456**(-1) )
      bmp87(jdu,2,2,h56,2,1) = bmp87(jdu,2,2,h56,2,1) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * ( 4.D0*za(p2,p7)*za(p1,p5)*za(p4,p5)*zb(p2,
     &    p8)**2*zb(p5,p6)*zab2(p7,p2,p8,p3)*iza(p7,p8)*izb(p7,p8)*
     &    s27**(-1)*s278**(-1)*s456**(-1) - 4.D0*za(p1,p7)**2*za(p3,p5)
     &    *zb(p2,p3)*zb(p1,p8)*zb(p3,p6)*zab2(p4,p1,p7,p8)*iza(p7,p8)*
     &    izb(p7,p8)*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p1,p7)**
     &    2*za(p4,p5)*zb(p2,p3)*zb(p1,p8)*zb(p4,p6)*zab2(p4,p1,p7,p8)*
     &    iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*s456**(-1) + 4.D0*
     &    za(p1,p7)**2*za(p4,p5)*zb(p2,p3)*zb(p1,p8)*zb(p5,p6)*zab2(p5,
     &    p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s456**(-1) + 4.D0*za(p1,p7)**2*za(p5,p6)*zb(p2,p6)*zb(p1,p8)*
     &    zb(p3,p6)*zab2(p4,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s356**(-1) )
      bmp87(jdu,2,2,h56,2,2)= + gamz3456(jdu,2,2)*gamll56(2,h56) * ( 
     &     - 4.D0*za(p2,p1)**2*za(p3,p5)*zb(p2,p3)*zb(p1,p8)*zb(p3,p6)*
     &    zab2(p4,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*s18**(-1)*s178**(-1)*
     &    s356**(-1) + 4.D0*za(p2,p1)**2*za(p4,p5)*zb(p2,p3)*zb(p1,p8)*
     &    zb(p4,p6)*zab2(p4,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*s18**(-1)*
     &    s178**(-1)*s456**(-1) + 4.D0*za(p2,p1)**2*za(p4,p5)*zb(p2,p3)
     &    *zb(p1,p8)*zb(p5,p6)*zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*
     &    s18**(-1)*s178**(-1)*s456**(-1) + 4.D0*za(p2,p1)**2*za(p5,p6)
     &    *zb(p2,p6)*zb(p1,p8)*zb(p3,p6)*zab2(p4,p1,p8,p7)*iza(p2,p7)*
     &    iza(p2,p8)*s18**(-1)*s178**(-1)*s356**(-1) + 4.D0*za(p2,p1)*
     &    za(p1,p4)*za(p3,p5)*zb(p2,p7)*zb(p1,p8)*zb(p3,p6)*zb(p3,p7)*
     &    iza(p2,p8)*s18**(-1)*s27**(-1)*s356**(-1) - 4.D0*za(p2,p1)*
     &    za(p1,p4)*za(p4,p5)*zb(p2,p7)*zb(p1,p8)*zb(p3,p7)*zb(p4,p6)*
     &    iza(p2,p8)*s18**(-1)*s27**(-1)*s456**(-1) - 4.D0*za(p2,p1)*
     &    za(p1,p4)*za(p5,p6)*zb(p2,p7)*zb(p1,p8)*zb(p3,p6)*zb(p6,p7)*
     &    iza(p2,p8)*s18**(-1)*s27**(-1)*s356**(-1) )
      bmp87(jdu,2,2,h56,2,2) = bmp87(jdu,2,2,h56,2,2) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * (  - 4.D0*za(p2,p1)*za(p1,p5)*za(p4,p5)*zb(
     &    p2,p7)*zb(p1,p8)*zb(p3,p7)*zb(p5,p6)*iza(p2,p8)*s18**(-1)*
     &    s27**(-1)*s456**(-1) + 4.D0*za(p2,p1)*za(p3,p5)*zb(p2,p3)*zb(
     &    p3,p6)*zb(p7,p8)*zab2(p4,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p2,p1)*za(p3,p5)*zb(p2,p3)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p4,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p2,p1)*za(p4,p5)*zb(p2,p3)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p4,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p2,p1)*za(p4,p5)*zb(p2,p3)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p4,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p2,p1)*za(p4,p5)*zb(p2,p3)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p5,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s456**(-1) - 4.D0*za(p2,p1)*za(p4,p5)*zb(p2,p3)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p5,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s456**(-1) )
      bmp87(jdu,2,2,h56,2,2) = bmp87(jdu,2,2,h56,2,2) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * (  - 4.D0*za(p2,p1)*za(p5,p6)*zb(p2,p6)*zb(
     &    p3,p6)*zb(p7,p8)*zab2(p4,p1,p7,p8)*iza(p2,p7)*s78**(-1)*
     &    s178**(-1)*s356**(-1) - 4.D0*za(p2,p1)*za(p5,p6)*zb(p2,p6)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p4,p1,p8,p7)*iza(p2,p8)*s78**(-1)*
     &    s178**(-1)*s356**(-1) + 4.D0*za(p1,p4)*za(p3,p5)*zb(p2,p7)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p4)*za(p3,p5)*zb(p2,p7)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) + 4.D0*za(p1,p4)*za(p3,p5)*zb(p2,p8)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p4)*za(p4,p5)*zb(p2,p7)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p4)*za(p4,p5)*zb(p2,p7)*
     &    zb(p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) )
      bmp87(jdu,2,2,h56,2,2) = bmp87(jdu,2,2,h56,2,2) + gamz3456(jdu,2,
     & 2)*gamll56(2,h56) * (  - 4.D0*za(p1,p4)*za(p4,p5)*zb(p2,p8)*zb(
     &    p4,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p4)*za(p5,p6)*zb(p2,p7)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p4)*za(p5,p6)*zb(p2,p7)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p4)*za(p5,p6)*zb(p2,p8)*
     &    zb(p3,p6)*zb(p7,p8)*zab2(p2,p7,p8,p6)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s356**(-1) - 4.D0*za(p1,p5)*za(p4,p5)*zb(p2,p7)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p8)*s78**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p5)*za(p4,p5)*zb(p2,p7)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p8)*s27**(-1)*
     &    s278**(-1)*s456**(-1) - 4.D0*za(p1,p5)*za(p4,p5)*zb(p2,p8)*
     &    zb(p5,p6)*zb(p7,p8)*zab2(p2,p7,p8,p3)*iza(p2,p7)*s78**(-1)*
     &    s278**(-1)*s456**(-1) )
      endif
      do h34=1,2
      if (h34 == 1) then
      p3=i3
      p4=i4
      elseif (h34 == 2) then
      p3=i4
      p4=i3
      endif


      if (j78==1) then
C---  amp78(jdu,h12,h34,h56,h7,h8)
      amp78(jdu,1,h34,h56,1,1)= + gamz34(jdu,1,h34)*gamz56(jdu,1,h56)
     &  * (  - 4.D0*za(p2,p3)*za(p5,p6)*za(p7,p8)*zb(p1,p6)**2*zba2(p4,
     &    p2,p3,p7)*izb(p1,p8)*s78**(-1)*s156**(-1)*s234**(-1) - 4.D0*
     &    za(p2,p3)*za(p5,p6)*za(p7,p8)*zb(p1,p6)**2*zba2(p4,p2,p3,p8)*
     &    izb(p1,p7)*s78**(-1)*s156**(-1)*s234**(-1) - 4.D0*za(p2,p3)*
     &    za(p5,p6)*zb(p1,p6)**2*zba2(p4,p2,p3,p7)*izb(p1,p7)*izb(p1,p8
     &    )*zba4(p1,p2,p3,p4,p7,p8)*s156**(-1)*s234**(-1)*s2347**(-1)
     &     - 4.D0*za(p2,p7)*za(p5,p6)*zb(p1,p6)**2*zba2(p1,p2,p7,p3)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p4,p2,p3,p7,p8)*s27**(-1)*
     &    s156**(-1)*s2347**(-1) + 4.D0*za(p2,p7)*za(p7,p8)*zb(p1,p6)*
     &    zba2(p4,p1,p6,p5)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*s78**(-1)*
     &    s278**(-1)*s156**(-1) + 4.D0*za(p2,p7)*zb(p1,p6)*zba2(p1,p2,
     &    p7,p8)*zba2(p4,p1,p6,p5)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,
     &    p8,p3)*s27**(-1)*s278**(-1)*s156**(-1) + 4.D0*za(p2,p8)*za(p7
     &    ,p8)*zb(p1,p6)*zba2(p4,p1,p6,p5)*izb(p1,p7)*zba3(p1,p2,p7,p8,
     &    p3)*s78**(-1)*s278**(-1)*s156**(-1) )
      amp78(jdu,1,h34,h56,2,1)= + gamz34(jdu,1,h34)*gamz56(jdu,1,h56)
     &  * ( 4.D0*za(p2,p3)*za(p2,p8)*za(p1,p8)*zb(p2,p7)*zb(p1,p6)*zb(
     &    p1,p7)*iza(p7,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,p5)*s18**(-1)*
     &    s27**(-1)*s2347**(-1) - 4.D0*za(p2,p3)*za(p2,p8)*zb(p2,p7)*
     &    zb(p1,p6)*zba2(p7,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*zba3(p4,p2,
     &    p3,p7,p8)*s27**(-1)*s156**(-1)*s2347**(-1) - 4.D0*za(p2,p3)*
     &    za(p1,p8)*zb(p1,p6)*zb(p1,p7)*zba2(p4,p2,p3,p8)*iza(p7,p8)*
     &    izb(p7,p8)*zba3(p7,p2,p3,p4,p5)*s18**(-1)*s234**(-1)*
     &    s2347**(-1) + 4.D0*za(p2,p3)*za(p1,p8)*zb(p1,p7)**2*zba2(p4,
     &    p2,p3,p5)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s234**(-1) + 4.D0*za(p2,p3)*zb(p1,p6)*zba2(p4,p2,
     &    p3,p8)*zba2(p7,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*zba3(p7,p2,p3,
     &    p4,p8)*s156**(-1)*s234**(-1)*s2347**(-1) - 4.D0*za(p2,p8)**2*
     &    zb(p2,p7)*zb(p1,p6)*zba2(p4,p1,p6,p5)*zba2(p7,p2,p8,p3)*iza(
     &    p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*s156**(-1) )
      amp78(jdu,1,h34,h56,1,2)= + gamz34(jdu,1,h34)*gamz56(jdu,1,h56)
     &  * ( 4.D0*za(p2,p3)*za(p1,p7)*zb(p1,p8)**2*zba2(p4,p2,p3,p5)*
     &    zba2(p6,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s234**(-1) + 4.D0*za(p2,p3)*zb(p1,p6)*zba2(p4,p2,p3,p7)*zba2(
     &    p8,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*zba3(p8,p2,p3,p4,p7)*
     &    s156**(-1)*s234**(-1)*s2347**(-1) + 4.D0*za(p2,p3)*zb(p1,p8)*
     &    zba2(p4,p2,p3,p7)*zba2(p6,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*
     &    zba4(p8,p2,p3,p4,p7,p5)*s18**(-1)*s234**(-1)*s2347**(-1) - 4.D
     &    0*za(p2,p7)**2*zb(p2,p8)*zb(p1,p6)*zba2(p4,p1,p6,p5)*zba2(p8,
     &    p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s156**(-1) + 4.D0*za(p2,p7)*zb(p1,p6)*zba2(p4,p2,p3,p7)*zba2(
     &    p8,p2,p7,p3)*zba2(p8,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*
     &    s27**(-1)*s156**(-1)*s2347**(-1) + 4.D0*za(p2,p7)*zb(p1,p8)*
     &    zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*
     &    zba3(p4,p2,p3,p7,p5)*s18**(-1)*s27**(-1)*s2347**(-1) )
      amp78(jdu,1,h34,h56,2,2)= + gamz34(jdu,1,h34)*gamz56(jdu,1,h56)
     &  * ( 4.D0*za(p2,p3)**2*zb(p1,p6)*zb(p3,p4)*zb(p7,p8)*zba2(p7,p1,
     &    p6,p5)*iza(p2,p8)*s78**(-1)*s156**(-1)*s234**(-1) + 4.D0*za(
     &    p2,p3)**2*zb(p1,p6)*zb(p3,p4)*zb(p7,p8)*zba2(p8,p1,p6,p5)*
     &    iza(p2,p7)*s78**(-1)*s156**(-1)*s234**(-1) + 4.D0*za(p2,p3)**
     &    2*zb(p1,p6)*zb(p3,p4)*zba2(p7,p3,p4,p2)*zba2(p8,p1,p6,p5)*
     &    iza(p2,p7)*iza(p2,p8)*s156**(-1)*s234**(-1)*s2347**(-1) + 4.D0
     &    *za(p2,p3)**2*zb(p1,p8)*zb(p3,p4)*zba2(p6,p1,p8,p2)*iza(p2,p7
     &    )*iza(p2,p8)*zba3(p7,p2,p3,p4,p5)*s18**(-1)*s234**(-1)*
     &    s2347**(-1) + 4.D0*za(p2,p3)*zb(p1,p7)*zb(p7,p8)*zba2(p4,p2,
     &    p3,p5)*iza(p2,p8)*zba3(p6,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*
     &    s234**(-1) + 4.D0*za(p2,p3)*zb(p1,p8)*zb(p7,p8)*zba2(p4,p2,p3
     &    ,p5)*iza(p2,p7)*zba3(p6,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*
     &    s234**(-1) - 4.D0*za(p2,p3)*zb(p1,p8)*zba2(p4,p2,p3,p5)*zba2(
     &    p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*zba3(p6,p1,p7,p8,p2)*
     &    s18**(-1)*s178**(-1)*s234**(-1) )
      amp78(jdu,2,h34,h56,1,1)= + gamz34(jdu,2,h34)*gamz56(jdu,2,h56)
     &  * (  - 4.D0*za(p2,p7)*za(p1,p5)*zb(p2,p1)**2*zab2(p3,p1,p5,p6)*
     &    zab2(p8,p2,p7,p4)*izb(p1,p7)*izb(p1,p8)*s27**(-1)*s278**(-1)*
     &    s156**(-1) - 4.D0*za(p2,p7)*za(p1,p5)*zb(p2,p1)*zb(p2,p4)*
     &    zab2(p8,p1,p5,p6)*izb(p1,p7)*izb(p1,p8)*zab3(p3,p2,p4,p7,p1)*
     &    s27**(-1)*s156**(-1)*s2347**(-1) + 4.D0*za(p2,p7)*za(p1,p8)*
     &    za(p5,p8)*zb(p2,p1)*zb(p2,p4)*izb(p1,p7)*zab3(p3,p2,p4,p7,p6)
     &    *s18**(-1)*s27**(-1)*s2347**(-1) + 4.D0*za(p1,p5)*za(p7,p8)*
     &    zb(p2,p1)*zab2(p3,p1,p5,p6)*zab2(p7,p2,p8,p4)*izb(p1,p8)*
     &    s78**(-1)*s278**(-1)*s156**(-1) + 4.D0*za(p1,p5)*za(p7,p8)*
     &    zb(p2,p1)*zab2(p3,p1,p5,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*
     &    s78**(-1)*s278**(-1)*s156**(-1) + 4.D0*za(p1,p5)*za(p7,p8)*
     &    zb(p2,p4)*zab2(p3,p2,p4,p1)*zab2(p7,p1,p5,p6)*izb(p1,p8)*
     &    s78**(-1)*s156**(-1)*s234**(-1) + 4.D0*za(p1,p5)*za(p7,p8)*
     &    zb(p2,p4)*zab2(p3,p2,p4,p1)*zab2(p8,p1,p5,p6)*izb(p1,p7)*
     &    s78**(-1)*s156**(-1)*s234**(-1) )
      amp78(jdu,2,h34,h56,1,1) = amp78(jdu,2,h34,h56,1,1) + gamz34(jdu,
     & 2,h34)*gamz56(jdu,2,h56) * ( 4.D0*za(p1,p5)*zb(p2,p4)*zab2(p3,p2
     &    ,p4,p1)*zab2(p8,p1,p5,p6)*izb(p1,p7)*izb(p1,p8)*zab3(p7,p2,p3
     &    ,p4,p1)*s156**(-1)*s234**(-1)*s2347**(-1) + 4.D0*za(p1,p7)*
     &    za(p7,p8)*zb(p2,p4)*zab2(p3,p2,p4,p6)*zab2(p5,p7,p8,p1)*izb(
     &    p1,p8)*s78**(-1)*s178**(-1)*s234**(-1) - 4.D0*za(p1,p8)*za(p5
     &    ,p8)*zb(p2,p4)*zab2(p3,p2,p4,p1)*izb(p1,p7)*zab3(p7,p2,p3,p4,
     &    p6)*s18**(-1)*s234**(-1)*s2347**(-1) + 4.D0*za(p1,p8)*za(p7,
     &    p8)*zb(p2,p4)*zab2(p3,p2,p4,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*
     &    s78**(-1)*s178**(-1)*s234**(-1) + 4.D0*za(p1,p8)*za(p7,p8)*
     &    zb(p2,p4)*zab2(p3,p2,p4,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*
     &    s18**(-1)*s178**(-1)*s234**(-1) )
      amp78(jdu,2,h34,h56,2,1)= + gamz34(jdu,2,h34)*gamz56(jdu,2,h56)
     &  * (  - 4.D0*za(p2,p8)*za(p1,p5)*zb(p2,p7)**2*zab2(p3,p1,p5,p6)*
     &    zab2(p8,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s156**(-1) + 4.D0*za(p1,p5)*zb(p2,p4)*zab2(p3,p2,p4,p7)*zab2(
     &    p8,p1,p5,p6)*iza(p7,p8)*izb(p7,p8)*zab3(p8,p2,p3,p4,p7)*
     &    s156**(-1)*s234**(-1)*s2347**(-1) + 4.D0*za(p1,p5)*zb(p2,p7)*
     &    zab2(p3,p2,p4,p7)*zab2(p8,p2,p7,p4)*zab2(p8,p1,p5,p6)*iza(p7,
     &    p8)*izb(p7,p8)*s27**(-1)*s156**(-1)*s2347**(-1) + 4.D0*za(p1,
     &    p8)**2*zb(p2,p4)*zb(p1,p7)*zab2(p3,p2,p4,p6)*zab2(p5,p1,p8,p7
     &    )*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*s234**(-1) + 4.D0
     &    *za(p1,p8)*zb(p2,p4)*zab2(p3,p2,p4,p7)*zab2(p5,p1,p8,p7)*iza(
     &    p7,p8)*izb(p7,p8)*zab4(p8,p2,p3,p4,p7,p6)*s18**(-1)*
     &    s234**(-1)*s2347**(-1) + 4.D0*za(p1,p8)*zb(p2,p7)*zab2(p5,p1,
     &    p8,p7)*zab2(p8,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*zab3(p3,p2,p4,
     &    p7,p6)*s18**(-1)*s27**(-1)*s2347**(-1) )
      amp78(jdu,2,h34,h56,1,2)= + gamz34(jdu,2,h34)*gamz56(jdu,2,h56)
     &  * ( 4.D0*za(p2,p7)*za(p1,p5)*za(p1,p7)*zb(p2,p4)*zb(p2,p8)*zb(
     &    p1,p8)*iza(p7,p8)*izb(p7,p8)*zab3(p3,p2,p4,p7,p6)*s18**(-1)*
     &    s27**(-1)*s2347**(-1) - 4.D0*za(p2,p7)*za(p1,p5)*zb(p2,p4)*
     &    zb(p2,p8)*zab2(p7,p1,p5,p6)*iza(p7,p8)*izb(p7,p8)*zab3(p3,p2,
     &    p4,p7,p8)*s27**(-1)*s156**(-1)*s2347**(-1) - 4.D0*za(p2,p7)*
     &    za(p1,p5)*zb(p2,p8)**2*zab2(p3,p1,p5,p6)*zab2(p7,p2,p8,p4)*
     &    iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*s156**(-1) - 4.D0*
     &    za(p1,p5)*za(p1,p7)*zb(p2,p4)*zb(p1,p8)*zab2(p3,p2,p4,p8)*
     &    iza(p7,p8)*izb(p7,p8)*zab3(p7,p2,p3,p4,p6)*s18**(-1)*
     &    s234**(-1)*s2347**(-1) + 4.D0*za(p1,p5)*zb(p2,p4)*zab2(p3,p2,
     &    p4,p8)*zab2(p7,p1,p5,p6)*iza(p7,p8)*izb(p7,p8)*zab3(p7,p2,p3,
     &    p4,p8)*s156**(-1)*s234**(-1)*s2347**(-1) + 4.D0*za(p1,p7)**2*
     &    zb(p2,p4)*zb(p1,p8)*zab2(p3,p2,p4,p6)*zab2(p5,p1,p7,p8)*iza(
     &    p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*s234**(-1) )
      amp78(jdu,2,h34,h56,2,2)= + gamz34(jdu,2,h34)*gamz56(jdu,2,h56)
     &  * ( 4.D0*za(p2,p1)**2*zb(p2,p4)*zb(p1,p8)*zab2(p3,p2,p4,p6)*
     &    zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*s18**(-1)*s178**(-1)*
     &    s234**(-1) + 4.D0*za(p2,p1)*za(p1,p5)*zb(p2,p4)*zb(p1,p8)*
     &    zab2(p3,p2,p4,p7)*iza(p2,p7)*iza(p2,p8)*zab3(p2,p3,p4,p7,p6)*
     &    s18**(-1)*s234**(-1)*s2347**(-1) - 4.D0*za(p2,p1)*za(p1,p5)*
     &    zb(p2,p7)*zb(p1,p8)*zb(p4,p7)*iza(p2,p8)*zab3(p3,p2,p4,p7,p6)
     &    *s18**(-1)*s27**(-1)*s2347**(-1) - 4.D0*za(p2,p1)*zb(p2,p4)*
     &    zb(p7,p8)*zab2(p3,p2,p4,p6)*zab2(p5,p1,p7,p8)*iza(p2,p7)*
     &    s78**(-1)*s178**(-1)*s234**(-1) - 4.D0*za(p2,p1)*zb(p2,p4)*
     &    zb(p7,p8)*zab2(p3,p2,p4,p6)*zab2(p5,p1,p8,p7)*iza(p2,p8)*
     &    s78**(-1)*s178**(-1)*s234**(-1) + 4.D0*za(p1,p5)*zb(p2,p4)*
     &    zb(p7,p8)*zab2(p2,p1,p5,p6)*zab2(p3,p2,p4,p7)*iza(p2,p8)*
     &    s78**(-1)*s156**(-1)*s234**(-1) + 4.D0*za(p1,p5)*zb(p2,p4)*
     &    zb(p7,p8)*zab2(p2,p1,p5,p6)*zab2(p3,p2,p4,p8)*iza(p2,p7)*
     &    s78**(-1)*s156**(-1)*s234**(-1) )
      amp78(jdu,2,h34,h56,2,2) = amp78(jdu,2,h34,h56,2,2) + gamz34(jdu,
     & 2,h34)*gamz56(jdu,2,h56) * ( 4.D0*za(p1,p5)*zb(p2,p4)*zab2(p2,p1
     &    ,p5,p6)*zab2(p3,p2,p4,p7)*iza(p2,p7)*iza(p2,p8)*zab3(p2,p3,p4
     &    ,p7,p8)*s156**(-1)*s234**(-1)*s2347**(-1) - 4.D0*za(p1,p5)*
     &    zb(p2,p7)*zb(p4,p7)*zab2(p2,p1,p5,p6)*iza(p2,p8)*zab3(p3,p2,
     &    p4,p7,p8)*s27**(-1)*s156**(-1)*s2347**(-1) + 4.D0*za(p1,p5)*
     &    zb(p2,p7)*zb(p7,p8)*zab2(p2,p7,p8,p4)*zab2(p3,p1,p5,p6)*iza(
     &    p2,p8)*s78**(-1)*s278**(-1)*s156**(-1) + 4.D0*za(p1,p5)*zb(p2
     &    ,p7)*zb(p7,p8)*zab2(p2,p7,p8,p4)*zab2(p3,p1,p5,p6)*iza(p2,p8)
     &    *s27**(-1)*s278**(-1)*s156**(-1) + 4.D0*za(p1,p5)*zb(p2,p8)*
     &    zb(p7,p8)*zab2(p2,p7,p8,p4)*zab2(p3,p1,p5,p6)*iza(p2,p7)*
     &    s78**(-1)*s278**(-1)*s156**(-1) )
      elseif (j78==2) then
      amp87(jdu,1,h34,h56,1,1)= + gamz34(jdu,1,h34)*gamz56(jdu,1,h56)
     &  * (  - 4.D0*za(p2,p3)*za(p5,p6)*za(p7,p8)*zb(p1,p6)**2*zba2(p4,
     &    p2,p3,p7)*izb(p1,p8)*s78**(-1)*s156**(-1)*s234**(-1) - 4.D0*
     &    za(p2,p3)*za(p5,p6)*za(p7,p8)*zb(p1,p6)**2*zba2(p4,p2,p3,p8)*
     &    izb(p1,p7)*s78**(-1)*s156**(-1)*s234**(-1) - 4.D0*za(p2,p3)*
     &    za(p5,p6)*zb(p1,p6)**2*zba2(p4,p2,p3,p7)*izb(p1,p7)*izb(p1,p8
     &    )*zba4(p1,p2,p3,p4,p7,p8)*s156**(-1)*s234**(-1)*s2347**(-1)
     &     - 4.D0*za(p2,p7)*za(p5,p6)*zb(p1,p6)**2*zba2(p1,p2,p7,p3)*
     &    izb(p1,p7)*izb(p1,p8)*zba3(p4,p2,p3,p7,p8)*s27**(-1)*
     &    s156**(-1)*s2347**(-1) + 4.D0*za(p2,p7)*za(p7,p8)*zb(p1,p6)*
     &    zba2(p4,p1,p6,p5)*izb(p1,p8)*zba3(p1,p2,p7,p8,p3)*s78**(-1)*
     &    s278**(-1)*s156**(-1) + 4.D0*za(p2,p7)*zb(p1,p6)*zba2(p1,p2,
     &    p7,p8)*zba2(p4,p1,p6,p5)*izb(p1,p7)*izb(p1,p8)*zba3(p1,p2,p7,
     &    p8,p3)*s27**(-1)*s278**(-1)*s156**(-1) + 4.D0*za(p2,p8)*za(p7
     &    ,p8)*zb(p1,p6)*zba2(p4,p1,p6,p5)*izb(p1,p7)*zba3(p1,p2,p7,p8,
     &    p3)*s78**(-1)*s278**(-1)*s156**(-1) )
      amp87(jdu,1,h34,h56,1,2)= + gamz34(jdu,1,h34)*gamz56(jdu,1,h56)
     &  * ( 4.D0*za(p2,p3)*za(p2,p8)*za(p1,p8)*zb(p2,p7)*zb(p1,p6)*zb(
     &    p1,p7)*iza(p7,p8)*izb(p7,p8)*zba3(p4,p2,p3,p7,p5)*s18**(-1)*
     &    s27**(-1)*s2347**(-1) - 4.D0*za(p2,p3)*za(p2,p8)*zb(p2,p7)*
     &    zb(p1,p6)*zba2(p7,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*zba3(p4,p2,
     &    p3,p7,p8)*s27**(-1)*s156**(-1)*s2347**(-1) - 4.D0*za(p2,p3)*
     &    za(p1,p8)*zb(p1,p6)*zb(p1,p7)*zba2(p4,p2,p3,p8)*iza(p7,p8)*
     &    izb(p7,p8)*zba3(p7,p2,p3,p4,p5)*s18**(-1)*s234**(-1)*
     &    s2347**(-1) + 4.D0*za(p2,p3)*za(p1,p8)*zb(p1,p7)**2*zba2(p4,
     &    p2,p3,p5)*zba2(p6,p1,p7,p8)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*
     &    s178**(-1)*s234**(-1) + 4.D0*za(p2,p3)*zb(p1,p6)*zba2(p4,p2,
     &    p3,p8)*zba2(p7,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*zba3(p7,p2,p3,
     &    p4,p8)*s156**(-1)*s234**(-1)*s2347**(-1) - 4.D0*za(p2,p8)**2*
     &    zb(p2,p7)*zb(p1,p6)*zba2(p4,p1,p6,p5)*zba2(p7,p2,p8,p3)*iza(
     &    p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*s156**(-1) )
      amp87(jdu,1,h34,h56,2,1)= + gamz34(jdu,1,h34)*gamz56(jdu,1,h56)
     &  * ( 4.D0*za(p2,p3)*za(p1,p7)*zb(p1,p8)**2*zba2(p4,p2,p3,p5)*
     &    zba2(p6,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*
     &    s234**(-1) + 4.D0*za(p2,p3)*zb(p1,p6)*zba2(p4,p2,p3,p7)*zba2(
     &    p8,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*zba3(p8,p2,p3,p4,p7)*
     &    s156**(-1)*s234**(-1)*s2347**(-1) + 4.D0*za(p2,p3)*zb(p1,p8)*
     &    zba2(p4,p2,p3,p7)*zba2(p6,p1,p8,p7)*iza(p7,p8)*izb(p7,p8)*
     &    zba4(p8,p2,p3,p4,p7,p5)*s18**(-1)*s234**(-1)*s2347**(-1) - 4.D
     &    0*za(p2,p7)**2*zb(p2,p8)*zb(p1,p6)*zba2(p4,p1,p6,p5)*zba2(p8,
     &    p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s156**(-1) + 4.D0*za(p2,p7)*zb(p1,p6)*zba2(p4,p2,p3,p7)*zba2(
     &    p8,p2,p7,p3)*zba2(p8,p1,p6,p5)*iza(p7,p8)*izb(p7,p8)*
     &    s27**(-1)*s156**(-1)*s2347**(-1) + 4.D0*za(p2,p7)*zb(p1,p8)*
     &    zba2(p6,p1,p8,p7)*zba2(p8,p2,p7,p3)*iza(p7,p8)*izb(p7,p8)*
     &    zba3(p4,p2,p3,p7,p5)*s18**(-1)*s27**(-1)*s2347**(-1) )
      amp87(jdu,1,h34,h56,2,2)= + gamz34(jdu,1,h34)*gamz56(jdu,1,h56)
     &  * ( 4.D0*za(p2,p3)**2*zb(p1,p6)*zb(p3,p4)*zb(p7,p8)*zba2(p7,p1,
     &    p6,p5)*iza(p2,p8)*s78**(-1)*s156**(-1)*s234**(-1) + 4.D0*za(
     &    p2,p3)**2*zb(p1,p6)*zb(p3,p4)*zb(p7,p8)*zba2(p8,p1,p6,p5)*
     &    iza(p2,p7)*s78**(-1)*s156**(-1)*s234**(-1) + 4.D0*za(p2,p3)**
     &    2*zb(p1,p6)*zb(p3,p4)*zba2(p7,p3,p4,p2)*zba2(p8,p1,p6,p5)*
     &    iza(p2,p7)*iza(p2,p8)*s156**(-1)*s234**(-1)*s2347**(-1) + 4.D0
     &    *za(p2,p3)**2*zb(p1,p8)*zb(p3,p4)*zba2(p6,p1,p8,p2)*iza(p2,p7
     &    )*iza(p2,p8)*zba3(p7,p2,p3,p4,p5)*s18**(-1)*s234**(-1)*
     &    s2347**(-1) + 4.D0*za(p2,p3)*zb(p1,p7)*zb(p7,p8)*zba2(p4,p2,
     &    p3,p5)*iza(p2,p8)*zba3(p6,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*
     &    s234**(-1) + 4.D0*za(p2,p3)*zb(p1,p8)*zb(p7,p8)*zba2(p4,p2,p3
     &    ,p5)*iza(p2,p7)*zba3(p6,p1,p7,p8,p2)*s78**(-1)*s178**(-1)*
     &    s234**(-1) - 4.D0*za(p2,p3)*zb(p1,p8)*zba2(p4,p2,p3,p5)*zba2(
     &    p7,p1,p8,p2)*iza(p2,p7)*iza(p2,p8)*zba3(p6,p1,p7,p8,p2)*
     &    s18**(-1)*s178**(-1)*s234**(-1) )
      amp87(jdu,2,h34,h56,1,1)= + gamz34(jdu,2,h34)*gamz56(jdu,2,h56)
     &  * (  - 4.D0*za(p2,p7)*za(p1,p5)*zb(p2,p1)**2*zab2(p3,p1,p5,p6)*
     &    zab2(p8,p2,p7,p4)*izb(p1,p7)*izb(p1,p8)*s27**(-1)*s278**(-1)*
     &    s156**(-1) - 4.D0*za(p2,p7)*za(p1,p5)*zb(p2,p1)*zb(p2,p4)*
     &    zab2(p8,p1,p5,p6)*izb(p1,p7)*izb(p1,p8)*zab3(p3,p2,p4,p7,p1)*
     &    s27**(-1)*s156**(-1)*s2347**(-1) + 4.D0*za(p2,p7)*za(p1,p8)*
     &    za(p5,p8)*zb(p2,p1)*zb(p2,p4)*izb(p1,p7)*zab3(p3,p2,p4,p7,p6)
     &    *s18**(-1)*s27**(-1)*s2347**(-1) + 4.D0*za(p1,p5)*za(p7,p8)*
     &    zb(p2,p1)*zab2(p3,p1,p5,p6)*zab2(p7,p2,p8,p4)*izb(p1,p8)*
     &    s78**(-1)*s278**(-1)*s156**(-1) + 4.D0*za(p1,p5)*za(p7,p8)*
     &    zb(p2,p1)*zab2(p3,p1,p5,p6)*zab2(p8,p2,p7,p4)*izb(p1,p7)*
     &    s78**(-1)*s278**(-1)*s156**(-1) + 4.D0*za(p1,p5)*za(p7,p8)*
     &    zb(p2,p4)*zab2(p3,p2,p4,p1)*zab2(p7,p1,p5,p6)*izb(p1,p8)*
     &    s78**(-1)*s156**(-1)*s234**(-1) + 4.D0*za(p1,p5)*za(p7,p8)*
     &    zb(p2,p4)*zab2(p3,p2,p4,p1)*zab2(p8,p1,p5,p6)*izb(p1,p7)*
     &    s78**(-1)*s156**(-1)*s234**(-1) )
      amp87(jdu,2,h34,h56,1,1) = amp87(jdu,2,h34,h56,1,1) + gamz34(jdu,
     & 2,h34)*gamz56(jdu,2,h56) * ( 4.D0*za(p1,p5)*zb(p2,p4)*zab2(p3,p2
     &    ,p4,p1)*zab2(p8,p1,p5,p6)*izb(p1,p7)*izb(p1,p8)*zab3(p7,p2,p3
     &    ,p4,p1)*s156**(-1)*s234**(-1)*s2347**(-1) + 4.D0*za(p1,p7)*
     &    za(p7,p8)*zb(p2,p4)*zab2(p3,p2,p4,p6)*zab2(p5,p7,p8,p1)*izb(
     &    p1,p8)*s78**(-1)*s178**(-1)*s234**(-1) - 4.D0*za(p1,p8)*za(p5
     &    ,p8)*zb(p2,p4)*zab2(p3,p2,p4,p1)*izb(p1,p7)*zab3(p7,p2,p3,p4,
     &    p6)*s18**(-1)*s234**(-1)*s2347**(-1) + 4.D0*za(p1,p8)*za(p7,
     &    p8)*zb(p2,p4)*zab2(p3,p2,p4,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*
     &    s78**(-1)*s178**(-1)*s234**(-1) + 4.D0*za(p1,p8)*za(p7,p8)*
     &    zb(p2,p4)*zab2(p3,p2,p4,p6)*zab2(p5,p7,p8,p1)*izb(p1,p7)*
     &    s18**(-1)*s178**(-1)*s234**(-1) )
      amp87(jdu,2,h34,h56,1,2)= + gamz34(jdu,2,h34)*gamz56(jdu,2,h56)
     &  * (  - 4.D0*za(p2,p8)*za(p1,p5)*zb(p2,p7)**2*zab2(p3,p1,p5,p6)*
     &    zab2(p8,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*
     &    s156**(-1) + 4.D0*za(p1,p5)*zb(p2,p4)*zab2(p3,p2,p4,p7)*zab2(
     &    p8,p1,p5,p6)*iza(p7,p8)*izb(p7,p8)*zab3(p8,p2,p3,p4,p7)*
     &    s156**(-1)*s234**(-1)*s2347**(-1) + 4.D0*za(p1,p5)*zb(p2,p7)*
     &    zab2(p3,p2,p4,p7)*zab2(p8,p2,p7,p4)*zab2(p8,p1,p5,p6)*iza(p7,
     &    p8)*izb(p7,p8)*s27**(-1)*s156**(-1)*s2347**(-1) + 4.D0*za(p1,
     &    p8)**2*zb(p2,p4)*zb(p1,p7)*zab2(p3,p2,p4,p6)*zab2(p5,p1,p8,p7
     &    )*iza(p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*s234**(-1) + 4.D0
     &    *za(p1,p8)*zb(p2,p4)*zab2(p3,p2,p4,p7)*zab2(p5,p1,p8,p7)*iza(
     &    p7,p8)*izb(p7,p8)*zab4(p8,p2,p3,p4,p7,p6)*s18**(-1)*
     &    s234**(-1)*s2347**(-1) + 4.D0*za(p1,p8)*zb(p2,p7)*zab2(p5,p1,
     &    p8,p7)*zab2(p8,p2,p7,p4)*iza(p7,p8)*izb(p7,p8)*zab3(p3,p2,p4,
     &    p7,p6)*s18**(-1)*s27**(-1)*s2347**(-1) )
      amp87(jdu,2,h34,h56,2,1)= + gamz34(jdu,2,h34)*gamz56(jdu,2,h56)
     &  * ( 4.D0*za(p2,p7)*za(p1,p5)*za(p1,p7)*zb(p2,p4)*zb(p2,p8)*zb(
     &    p1,p8)*iza(p7,p8)*izb(p7,p8)*zab3(p3,p2,p4,p7,p6)*s18**(-1)*
     &    s27**(-1)*s2347**(-1) - 4.D0*za(p2,p7)*za(p1,p5)*zb(p2,p4)*
     &    zb(p2,p8)*zab2(p7,p1,p5,p6)*iza(p7,p8)*izb(p7,p8)*zab3(p3,p2,
     &    p4,p7,p8)*s27**(-1)*s156**(-1)*s2347**(-1) - 4.D0*za(p2,p7)*
     &    za(p1,p5)*zb(p2,p8)**2*zab2(p3,p1,p5,p6)*zab2(p7,p2,p8,p4)*
     &    iza(p7,p8)*izb(p7,p8)*s27**(-1)*s278**(-1)*s156**(-1) - 4.D0*
     &    za(p1,p5)*za(p1,p7)*zb(p2,p4)*zb(p1,p8)*zab2(p3,p2,p4,p8)*
     &    iza(p7,p8)*izb(p7,p8)*zab3(p7,p2,p3,p4,p6)*s18**(-1)*
     &    s234**(-1)*s2347**(-1) + 4.D0*za(p1,p5)*zb(p2,p4)*zab2(p3,p2,
     &    p4,p8)*zab2(p7,p1,p5,p6)*iza(p7,p8)*izb(p7,p8)*zab3(p7,p2,p3,
     &    p4,p8)*s156**(-1)*s234**(-1)*s2347**(-1) + 4.D0*za(p1,p7)**2*
     &    zb(p2,p4)*zb(p1,p8)*zab2(p3,p2,p4,p6)*zab2(p5,p1,p7,p8)*iza(
     &    p7,p8)*izb(p7,p8)*s18**(-1)*s178**(-1)*s234**(-1) )
      amp87(jdu,2,h34,h56,2,2)= + gamz34(jdu,2,h34)*gamz56(jdu,2,h56)
     &  * ( 4.D0*za(p2,p1)**2*zb(p2,p4)*zb(p1,p8)*zab2(p3,p2,p4,p6)*
     &    zab2(p5,p1,p8,p7)*iza(p2,p7)*iza(p2,p8)*s18**(-1)*s178**(-1)*
     &    s234**(-1) + 4.D0*za(p2,p1)*za(p1,p5)*zb(p2,p4)*zb(p1,p8)*
     &    zab2(p3,p2,p4,p7)*iza(p2,p7)*iza(p2,p8)*zab3(p2,p3,p4,p7,p6)*
     &    s18**(-1)*s234**(-1)*s2347**(-1) - 4.D0*za(p2,p1)*za(p1,p5)*
     &    zb(p2,p7)*zb(p1,p8)*zb(p4,p7)*iza(p2,p8)*zab3(p3,p2,p4,p7,p6)
     &    *s18**(-1)*s27**(-1)*s2347**(-1) - 4.D0*za(p2,p1)*zb(p2,p4)*
     &    zb(p7,p8)*zab2(p3,p2,p4,p6)*zab2(p5,p1,p7,p8)*iza(p2,p7)*
     &    s78**(-1)*s178**(-1)*s234**(-1) - 4.D0*za(p2,p1)*zb(p2,p4)*
     &    zb(p7,p8)*zab2(p3,p2,p4,p6)*zab2(p5,p1,p8,p7)*iza(p2,p8)*
     &    s78**(-1)*s178**(-1)*s234**(-1) + 4.D0*za(p1,p5)*zb(p2,p4)*
     &    zb(p7,p8)*zab2(p2,p1,p5,p6)*zab2(p3,p2,p4,p7)*iza(p2,p8)*
     &    s78**(-1)*s156**(-1)*s234**(-1) + 4.D0*za(p1,p5)*zb(p2,p4)*
     &    zb(p7,p8)*zab2(p2,p1,p5,p6)*zab2(p3,p2,p4,p8)*iza(p2,p7)*
     &    s78**(-1)*s156**(-1)*s234**(-1) )
      amp87(jdu,2,h34,h56,2,2) = amp87(jdu,2,h34,h56,2,2) + gamz34(jdu,
     & 2,h34)*gamz56(jdu,2,h56) * ( 4.D0*za(p1,p5)*zb(p2,p4)*zab2(p2,p1
     &    ,p5,p6)*zab2(p3,p2,p4,p7)*iza(p2,p7)*iza(p2,p8)*zab3(p2,p3,p4
     &    ,p7,p8)*s156**(-1)*s234**(-1)*s2347**(-1) - 4.D0*za(p1,p5)*
     &    zb(p2,p7)*zb(p4,p7)*zab2(p2,p1,p5,p6)*iza(p2,p8)*zab3(p3,p2,
     &    p4,p7,p8)*s27**(-1)*s156**(-1)*s2347**(-1) + 4.D0*za(p1,p5)*
     &    zb(p2,p7)*zb(p7,p8)*zab2(p2,p7,p8,p4)*zab2(p3,p1,p5,p6)*iza(
     &    p2,p8)*s78**(-1)*s278**(-1)*s156**(-1) + 4.D0*za(p1,p5)*zb(p2
     &    ,p7)*zb(p7,p8)*zab2(p2,p7,p8,p4)*zab2(p3,p1,p5,p6)*iza(p2,p8)
     &    *s27**(-1)*s278**(-1)*s156**(-1) + 4.D0*za(p1,p5)*zb(p2,p8)*
     &    zb(p7,p8)*zab2(p2,p7,p8,p4)*zab2(p3,p1,p5,p6)*iza(p2,p7)*
     &    s78**(-1)*s278**(-1)*s156**(-1) )
      endif
      enddo  !jdu
      enddo  !h34
      enddo  !h56
      enddo  !j78
      amp78=amp78+bmp78
      amp87=amp87+bmp87
      return
      end
