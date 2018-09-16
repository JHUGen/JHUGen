      double complex function a6vQLslc(st,j1,j2,j3,j4,j5,j6,za,zb)
****************************************************
* virtual amplitude for
* 0->q(p1)+qb(p2)+glu(p3)+gam(p4)+lb*p5)+l(p6)
* where one photon is coming from quark line
* and another photon coming from lepton line
****************************************************
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      character*14 st
      integer j1,j2,j3,j4,j5,j6
      double complex a6treeQLslc,vQLslc,fQLslc
c-----
      a6vQLslc = a6treeQLslc(st,j1,j2,j3,j4,j5,j6,za,zb)
     .           *vQLslc(st,j1,j2,j3,j4,j5,j6)
c-----
      a6vQLslc = a6vQLslc+fQLslc(st,j1,j2,j3,j4,j5,j6,za,zb)
c-----
      return
      end


      double complex function vQLslc(st,j1,j2,j3,j4,j5,j6)
c-----divergent part
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'epinv2.f'
      integer j1,j2,j3,j4,j5,j6
      double complex Lnrat,xl12,xl456
      character*14 st
      double precision t
c-----
      xl12=Lnrat(musq,-s(j1,j2))
      xl456=Lnrat(musq,-t(j4,j5,j6))
c-----
      vQLslc = -(7D0/2D0)
     .         -(epinv*epinv2+epinv*xl12+half*xl12**2)
     .         -(3D0/2D0)*(epinv+xl456)
c-----
      return
      end
      

      double complex function fQLslc(st,j1,j2,j3,j4,j5,j6,za,zb)
c-----finite part
      implicit none
      character*14 st
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      double complex Lsm1,L0,L1
      double precision t
c-----
      if(st.eq.'q+qb-g+g+lb-l+') then
c-----(+-++-+)
      fQLslc = 
     .      + za(j2,j5)**2/(za(j1,j3)*za(j2,j3)*za(j4,j5)*za(j4,j6))
     .        *Lsm1(-s(j1,j2),-t(j4,j5,j6),-s(j1,j3),-t(j4,j5,j6))
     .      + za(j1,j2)**2*za(j3,j5)**2
     .        /(za(j1,j3)**3*za(j2,j3)*za(j4,j5)*za(j4,j6))
     .        *Lsm1(-s(j1,j2),-t(j4,j5,j6),-s(j2,j3),-t(j4,j5,j6))
     .      - ( 2D0*s(j1,j3)*za(j1,j5)*za(j2,j5)
     .         -za(j2,j3)*zb(j3,j1)*za(j1,j5)**2 )
     .        /(za(j1,j3)**2*za(j4,j5)*za(j4,j6))
     .        *L0(-t(j4,j5,j6),-s(j2,j3))/s(j2,j3)
     .      - zb(j1,j3)**2*za(j2,j3)*za(j1,j5)**2
     .        /(2D0*za(j1,j3)*za(j4,j5)*za(j4,j6))
     .        *L1(-t(j4,j5,j6),-s(j2,j3))/s(j2,j3)**2
     .      - za(j2,j1)*zb(j1,j3)*za(j1,j5)*za(j3,j5)
     .        /(za(j1,j3)**2*za(j4,j5)*za(j4,j6))
     .        *L0(-t(j4,j5,j6),-s(j1,j2))/s(j1,j2)
     .      + za(j2,j1)*zb(j1,j3)*za(j5,j3)
     .        *( za(j5,j4)*zb(j4,j3) + za(j5,j6)*zb(j6,j3) )
     .        /(za(j1,j3)*za(j4,j5)*za(j4,j6))
     .        *L1(-t(j4,j5,j6),-s(j1,j2))/s(j1,j2)**2
     .      -  ( za(j5,j1)*zb(j1,j3)+za(j5,j2)*zb(j2,j3) )
     .        *( zb(j1,j3)*( za(j5,j4)*zb(j4,j2) + za(j5,j6)*zb(j6,j2) )
     .          +zb(j2,j3)*( za(j5,j4)*zb(j4,j1) + za(j5,j6)*zb(j6,j1) )
     .         )/(2D0*t(j4,j5,j6)
     .            *zb(j1,j2)*zb(j2,j3)*za(j1,j3)*za(j4,j5)*za(j4,j6))
c-----
      elseif(st.eq.'q+qb-g+g-lb-l+') then
c-----(+-+--+)
      fQLslc =
     .      -  (za(j2,j1)*zb(j1,j6)+ za(j2,j3)*zb(j3,j6))
     .        *(za(j2,j4)*zb(j4,j6)+ za(j2,j5)*zb(j5,j6))
     .        /(t(j4,j5,j6)*za(j1,j3)*za(j2,j3)*zb(j4,j5)*zb(j4,j6))
     .        *Lsm1(-s(j1,j2),-t(j4,j5,j6),-s(j1,j3),-t(j4,j5,j6))
     .      - za(j1,j2)**2*(za(j3,j1)*zb(j1,j6)+ za(j3,j2)*zb(j2,j6))
     .        *(za(j3,j4)*zb(j4,j6)+ za(j3,j5)*zb(j5,j6))
     .        /(t(j4,j5,j6)*za(j1,j3)**3*za(j2,j3)*zb(j4,j5)*zb(j4,j6))
     .        *Lsm1(-s(j1,j2),-t(j4,j5,j6),-s(j2,j3),-t(j4,j5,j6))
     .      + (za(j1,j2)*zb(j2,j6)+ za(j1,j3)*zb(j3,j6))
     .        *( 2D0*s(j1,j3)*(za(j2,j4)*zb(j4,j6)+ za(j2,j5)*zb(j5,j6))
     .          -za(j2,j3)*zb(j3,j1)
     .           *(za(j1,j4)*zb(j4,j6)+ za(j1,j5)*zb(j5,j6)) )
     .        /(t(j4,j5,j6)*za(j1,j3)**2*zb(j4,j5)*zb(j4,j6))
     .        *L0(-t(j4,j5,j6),-s(j2,j3))/s(j2,j3)
     .      + zb(j1,j3)**2*za(j2,j3)
     .        *(za(j1,j2)*zb(j2,j6)+ za(j1,j3)*zb(j3,j6))
     .        *(za(j1,j4)*zb(j4,j6)+ za(j1,j5)*zb(j5,j6))
     .        /(2D0*t(j4,j5,j6)*za(j1,j3)*zb(j4,j5)*zb(j4,j6))
     .        *L1(-t(j4,j5,j6),-s(j2,j3))/s(j2,j3)**2
     .      + za(j2,j1)*zb(j1,j3)
     .        *(za(j3,j1)*zb(j1,j6)+ za(j3,j2)*zb(j2,j6))
     .        *(za(j1,j4)*zb(j4,j6)+ za(j1,j5)*zb(j5,j6))
     .        /(t(j4,j5,j6)*za(j1,j3)**2*zb(j4,j5)*zb(j4,j6))
     .        *L0(-t(j4,j5,j6),-s(j1,j2))/s(j1,j2)
     .      + za(j2,j1)*zb(j1,j3)*zb(j6,j3)
     .        *( za(j3,j4)*zb(j4,j6) + za(j3,j5)*zb(j5,j6) )
     .        /(za(j1,j3)*zb(j4,j5)*zb(j4,j6))
     .        *L1(-t(j4,j5,j6),-s(j1,j2))/s(j1,j2)**2
     .      - zb(j3,j6)*( zb(j1,j3)*zb(j6,j2) + zb(j2,j3)*zb(j6,j1) )
     .        /(2D0*zb(j1,j2)*zb(j2,j3)*za(j1,j3)*zb(j4,j5)*zb(j4,j6))
c-----
      else
      write(6,*) 'unimplemented st',st
      stop
      endif
c-----
      return
      end

