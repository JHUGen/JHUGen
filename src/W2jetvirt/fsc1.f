      double complex function Fsc1(j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4,j5,j6
      double complex L0,L1,Lsm1_2mh,I3m,Lnrat,Ls1
      double precision t  

      Fsc1= 
     .-(za(j1,j5)*(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .((s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j1,j5)*zb(j1,j2)-
     .(-s(j1,j2)-s(j3,j4)+s(j5,j6))*za(j5,j6)*zb(j2,j6)))/
     .(2d0*(s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2d0*s(j1,j2)*
     .s(j5,j6)-
     .2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j3,j4)*za(j5,j6)*zb(j1,j2)*
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))-
     .(za(j1,j5)*zb(j1,j3)*(-(za(j1,j5)*
     .(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3)))-
     .za(j2,j4)*za(j5,j6)*zb(j2,j6)))/
     .(2d0*za(j3,j4)*za(j5,j6)*zb(j1,j2)*
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))+
     .((-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*zb(j4,j6)*
     .(-((-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j3,j4)*zb(j4,j6))+
     .(-s(j1,j2)-s(j3,j4)+s(j5,j6))*za(j3,j5)*zb(j5,j6)))/
     .(2d0*(s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2d0*s(j1,j2)*
     .s(j5,j6)-
     .2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j3,j4)*zb(j1,j2)*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*zb(j5,j6))-
     .(za(j2,j4)*zb(j4,j6)*(-((za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*
     .zb(j4,j6))-za(j3,j5)*zb(j1,j3)*zb(j5,j6)))/
     .(2d0*za(j3,j4)*zb(j1,j2)*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j
     .4))*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*zb(j5,j6))+
     .(L1(-t(j1,j2,j3),-s(j1,j2))*za(j1,j2)**2*zb(j1,j3)**2*
     .(-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6))**2)/
     .(2d0*s(j1,j2)**2*za(j1,j3)*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2
     .,j4))*
     .zb(j5,j6)*t(j1,j2,j3))-(zb(j1,j3)**2*
     .(za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6))**2)/
     .(2d0*za(j1,j3)*zb(j1,j2)*zb(j2,j3)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*zb(j5,j6)*t(j1,j2,j3))+
     .(za(j1,j2)*za(j2,j3)*zb(j1,j3)**2*
     .(za(j1,j2)*zb(j2,j6)+za(j1,j3)*zb(j3,j6))**2*
     .(-L1(-t(j1,j2,j3),-s(j2,j3))/(2d0*s(j2,j3)**2)+
     .Ls1(-s(j1,j2),-t(j1,j2,j3),-s(j2,j3),-t(j1,j2,j3))/t(j1,j2,j3)**2)
     .)/
     .(za(j1,j3)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*zb(j5,j6)*
     .t(j1,j2,j3))+(L0(-t(j1,j2,j3),-s(j1,j2))*za(j1,j2)*za(j2,j3)*zb(j1
     .,j3)*
     .zb(j4,j6)**2*t(j1,j2,j3))/
     .(s(j1,j2)*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**2*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*zb(j5,j6))+
     .(Lsm1_2mh(s(j3,j4),t(j1,j2,j3),s(j1,j2),s(j5,j6))*za(j2,j3)*
     .(-(za(j1,j2)*zb(j1,j4))+za(j2,j3)*zb(j3,j4))*zb(j4,j6)**2*t(j1,j2,
     .j3))/
     .((-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**3*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*zb(j5,j6))+
     .(za(j2,j4)*(-(za(j1,j2)*zb(j1,j4))+za(j2,j3)*zb(j3,j4))*zb(j4,j6)*
     .*2*
     .(-(L1(-s(j5,j6),-t(j1,j2,j3))*za(j2,j4))/(2d0*za(j2,j3)*t(j1,j2,j3
     .)**2)+
     .L0(-s(j5,j6),-t(j1,j2,j3))/
     .((-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*t(j1,j2,j3)))*
     .t(j1,j2,j3))/
     .((-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*zb(j5,j6))+
     .(Lnrat(-s(j3,j4),-s(j5,j6))*zb(j4,j6)*
     .((-3d0*za(j1,j2)*(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*zb(j3,j4)*
     .((s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j3,j5)+
     .2d0*za(j3,j4)*za(j5,j6)*zb(j4,j6)))/
     .(s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2d0*s(j1,j2)*s(j5,j
     .6)-
     .2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2-
     .(za(j2,j4)*zb(j1,j3)*zb(j4,j6))/(2d0*zb(j1,j2)*zb(j5,j6))-
     .((-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-((-(za(j1,j3)*zb(j1,j6))-za(j2,j3)*zb(j2,j6))*
     .zb(j3,j4))+(-3d0*s(j1,j2)-s(j3,j4)+s(j5,j6))*zb(j4,j6))
     .+(-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j2,j4)*zb(j3,j4)*
     .(-((za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*zb(j4,j6))-
     .za(j3,j5)*zb(j1,j3)*zb(j5,j6))-
     .2d0*zb(j4,j6)*t(j1,j2,j3)*
     .(2d0*(-(za(j2,j3)*zb(j1,j3))+za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))+
     .(za(j2,j4)*zb(j1,j3)+
     .(za(j2,j3)*zb(j1,j4)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3)))/
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4)))*
     .(t(j1,j2,j3)-t(j1,j2,j4))))/
     .(2d0*(s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2d0*s(j1,j2)*s(j5,j6)-2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*zb(j1,j2)
     .*
     .zb(j5,j6))))/
     .((-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))


      Fsc1=Fsc1
     .+(Lnrat(-s(j1,j2),-s(j5,j6))*za(j1,j2)*zb(j4,j6)*
     .((3d0*(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-((-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j3,j5)*zb(j3,j4))+
     .(-s(j1,j2)-s(j3,j4)+s(j5,j6))*za(j5,j6)*zb(j4,j6)))/
     .((s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2d0*s(j1,j2)*s(j5,j6)-2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4)))+
     .((-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*zb(j4,j6))/
     .(2d0*(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**2*zb(j5,j6))+
     .((-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(za(j3,j4)*zb(j4,j6)-za(j3,j5)*zb(j5,j6))-
     .(s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j2,j4)*
     .(-((za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3))*zb(j4,j6))-
     .za(j3,j5)*zb(j1,j3)*zb(j5,j6))-
     .((-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j3,j4)*
     .(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*zb(j4,j6)*
     .(t(j1,j2,j3)-t(j1,j2,j4)))/
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4)))/
     .(2d0*(s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2d0*s(j1,j2)*s(j5,j6)-2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j3,j4)
     .*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))*zb(j5,j6))))/
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))


      Fsc1=Fsc1
     .+(I3m(s(j1,j2),s(j3,j4),s(j5,j6))*za(j1,j2)*zb(j4,j6)*
     .(-(za(j2,j4)*za(j5,j6)*zb(j1,j3)*zb(j3,j6))+
     .((-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(za(j3,j5)*(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .zb(j3,j4)+3d0*za(j2,j3)*za(j5,j6)*zb(j1,j3)*zb(j4,j6)))/
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))+
     .(3d0*za(j5,j6)*(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*zb(j3,j4)*
     .(-((-s(j1,j2)+s(j3,j4)-s(j5,j6))*za(j3,j4)*zb(j4,j6))+
     .(-s(j1,j2)-s(j3,j4)+s(j5,j6))*za(j3,j5)*zb(j5,j6)))/
     .((s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2d0*s(j1,j2)*s(j5,j6)-2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4)))-
     .(za(j2,j3)*za(j5,j6)*zb(j1,j4)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*zb(j4,j6)*
     .(t(j1,j2,j3)-t(j1,j2,j4)))/
     .(-(za(j1,j3)*zb(j1,j4))-za(j2,j3)*zb(j2,j4))**2))/
     .((s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2d0*s(j1,j2)*s(j5,
     .j6)-
     .2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))


      Fsc1=Fsc1
     .+(L1(-t(j2,j3,j4),-s(j3,j4))*za(j2,j4)**2*
     .(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4))**2*zb(j3,j4)**2)/
     .(2d0*s(j3,j4)**2*za(j5,j6)*zb(j2,j4)*
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*t(j2,j3,j4))


      Fsc1=Fsc1
     .-(za(j2,j4)**2*(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4))**2)/
     .(2d0*za(j2,j3)*za(j3,j4)*za(j5,j6)*zb(j2,j4)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*t(j2,j3,j4))


      Fsc1=Fsc1
     .+(za(j2,j4)**2*zb(j2,j3)*zb(j3,j4)*
     .(-(za(j2,j5)*zb(j2,j4))-za(j3,j5)*zb(j3,j4))**2*
     .(-L1(-t(j2,j3,j4),-s(j2,j3))/(2d0*s(j2,j3)**2)+
     .Ls1(-s(j3,j4),-t(j2,j3,j4),-s(j2,j3),-t(j2,j3,j4))/t(j2,j3,j4)**2)
     .)/
     .(za(j5,j6)*zb(j2,j4)*(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))*
     .t(j2,j3,j4))


      Fsc1=Fsc1
     .+(L0(-t(j2,j3,j4),-s(j3,j4))*za(j1,j5)**2*za(j2,j4)*zb
     .(j2,j3)*
     .zb(j3,j4)*t(j2,j3,j4))/
     .(s(j3,j4)*za(j5,j6)*
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))**2*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))


      Fsc1=Fsc1
     .+(Lsm1_2mh(s(j1,j2),t(j2,j3,j4),s(j3,j4),s(j5,j6))*za(j1,j5)**2*
     .zb(j2,j3)*
     .(za(j1,j2)*zb(j2,j3)-za(j1,j4)*zb(j3,j4))*t(j2,j3,j4))/
     .(za(j5,j6)*(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))**3*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))


      Fsc1=Fsc1
     .+(za(j1,j5)**2*zb(j1,j3)*
     .(za(j1,j2)*zb(j2,j3)-za(j1,j4)*zb(j3,j4))*
     .(-(L1(-s(j5,j6),-t(j2,j3,j4))*zb(j1,j3))/
     .(2d0*zb(j2,j3)*t(j2,j3,j4)**2)+
     .L0(-s(j5,j6),-t(j2,j3,j4))/
     .((-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*t(j2,j3,j4)))*
     .t(j2,j3,j4))/
     .(za(j5,j6)*(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))


      Fsc1=Fsc1
     .-(I3m(s(j3,j4),s(j1,j2),s(j5,j6))*za(j1,j5)*zb(j3,j4)*
     .(za(j2,j4)*za(j2,j5)*zb(j1,j3)*zb(j5,j6)+
     .(3d0*za(j1,j2)*(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .((s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j1,j5)*zb(j1,j2)-
     .(-s(j1,j2)-s(j3,j4)+s(j5,j6))*za(j5,j6)*zb(j2,j6))*zb(j5,j6)
     .)/((s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2d0*s(j1,j2)*
     .s(j5,j6)-
     .2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4)))+
     .((-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j2)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .zb(j2,j6))-3d0*za(j1,j5)*za(j2,j4)*zb(j2,j3)*zb(j5,j6)))/
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))+
     .(za(j1,j4)*za(j1,j5)*(-(za(j2,j3)*zb(j1,j3))-
     .za(j2,j4)*zb(j1,j4))*zb(j2,j3)*zb(j5,j6)*
     .(-t(j1,j3,j4)+t(j2,j3,j4)))/
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))**2))/
     .((s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2d0*s(j1,j2)*s(j5,
     .j6)-
     .2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))


      Fsc1=Fsc1
     .-(Lnrat(-s(j3,j4),-s(j5,j6))*za(j1,j5)*zb(j3,j4)*
     .(-(za(j1,j5)*(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3)))/
     .(2d0*za(j5,j6)*(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))**2)+
     .(3d0*(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .((s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j1,j2)*zb(j2,j6)-
     .(-s(j1,j2)-s(j3,j4)+s(j5,j6))*za(j1,j5)*zb(j5,j6)))/
     .((s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2d0*s(j1,j2)*s(j5,j6)-2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2*
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4)))+
     .((-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-(za(j1,j5)*zb(j1,j2))+za(j5,j6)*zb(j2,j6))+
     .(-s(j1,j2)+s(j3,j4)-s(j5,j6))*zb(j1,j3)*
     .(-(za(j1,j5)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3)))-
     .za(j2,j4)*za(j5,j6)*zb(j2,j6))+
     .((s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j1,j5)*zb(j1,j2)*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .(-t(j1,j3,j4)+t(j2,j3,j4)))/
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4)))/
     .(2d0*(s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2d0*s(j1,j2)*s(j5,j6)-2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j5,j6)
     .*
     .zb(j1,j2)*(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4)))))/
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4))


      Fsc1=Fsc1
     .+(Lnrat(-s(j1,j2),-s(j5,j6))*za(j1,j5)*
     .(-(za(j1,j5)*za(j2,j4)*zb(j1,j3))/(2d0*za(j3,j4)*za(j5,j6))-
     .(3d0*za(j1,j2)*(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*zb(j3,j4)*
     .((-s(j1,j2)+s(j3,j4)-s(j5,j6))*zb(j2,j6)+
     .2d0*za(j1,j5)*zb(j1,j2)*zb(j5,j6)))/
     .(s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-2d0*s(j1,j2)*s(j5,j
     .6)-
     .2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)**2-
     .((-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(-(za(j1,j4)*zb(j1,j3))-za(j2,j4)*zb(j2,j3))*
     .((-s(j1,j2)-3d0*s(j3,j4)+s(j5,j6))*za(j1,j5)+
     .za(j1,j2)*(za(j3,j5)*zb(j2,j3)+za(j4,j5)*zb(j2,j4)))+
     .(s(j1,j2)-s(j3,j4)-s(j5,j6))*za(j1,j2)*zb(j1,j3)*
     .(-(za(j1,j5)*(za(j2,j4)*zb(j1,j2)+za(j3,j4)*zb(j1,j3)))-
     .za(j2,j4)*za(j5,j6)*zb(j2,j6))-
     .2d0*za(j1,j5)*t(j2,j3,j4)*
     .(2d0*(-(za(j2,j3)*zb(j1,j3))-za(j2,j4)*zb(j1,j4))*
     .(za(j1,j4)*zb(j1,j3)-za(j2,j4)*zb(j2,j3))+
     .(za(j2,j4)*zb(j1,j3)+
     .(za(j1,j4)*(-(za(j2,j3)*zb(j1,j3))-
     .za(j2,j4)*zb(j1,j4))*zb(j2,j3))/
     .(-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4)))*
     .(-t(j1,j3,j4)+t(j2,j3,j4))))/
     .(2d0*(s(j1,j2)**2-2d0*s(j1,j2)*s(j3,j4)+s(j3,j4)**2-
     .2d0*s(j1,j2)*s(j5,j6)-2d0*s(j3,j4)*s(j5,j6)+s(j5,j6)**2)*za(j3,j4)
     .*
     .za(j5,j6))))/
     .((-(za(j1,j3)*zb(j2,j3))-za(j1,j4)*zb(j2,j4))*
     .(za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)))

      return
      end
