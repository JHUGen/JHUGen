      function a6vQLlc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6vQLlc
****************************************************
* virtual amplitude for leading color
* 0->q(p1)+qb(p2)+glu(p3)+gam(p4)+lb*p5)+l(p6)
* where one photon is coming from quark line
****************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      character*14 st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a6treeQLlc,vQLlc,fQLlc
c-----
      a6vQLlc = a6treeQLlc(st,j1,j2,j3,j4,j5,j6,za,zb)
     &         *vQLlc(st,j1,j2,j3,j4,j5,j6)
c-----
      a6vQLlc = a6vQLlc+fQLlc(st,j1,j2,j3,j4,j5,j6,za,zb)
c-----
      return
      end


      function vQLlc(st,j1,j2,j3,j4,j5,j6)
      implicit none
      include 'types.f'
      complex(dp):: vQLlc
c-----divergent part

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scale.f'
      include 'epinv2.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: Lnrat,xl13,xl23,xl456
      character*14 st
      real(dp):: t
c-----
      xl13=Lnrat(musq,-s(j1,j3))
      xl23=Lnrat(musq,-s(j2,j3))
      xl456=Lnrat(musq,-t(j4,j5,j6))
c-----
      vQLlc = -3._dp
     &        -(epinv*epinv2+epinv*xl13+half*xl13**2)
     &        -(epinv*epinv2+epinv*xl23+half*xl23**2)
     &        -(3._dp/2._dp)*(epinv+xl456)
c-----
      return
      end


      function fQLlc(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: fQLlc
c-----finite part

      character*14 st
      integer:: j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      complex(dp):: Lsm1,L0,L1
      real(dp):: t
c-----
      if(st=='q+qb-g+g+lb-l+') then
c-----(+-++-+)
      fQLlc =
     &       - za(j2,j5)**2/(za(j1,j3)*za(j2,j3)*za(j4,j5)*za(j4,j6))
     &         *Lsm1(-s(j1,j3),-t(j4,j5,j6),-s(j2,j3),-t(j4,j5,j6))
     &       - za(j1,j2)*za(j2,j5)
     &         *( za(j5,j4)*zb(j4,j1) + za(j5,j6)*zb(j6,j1) )
     &         /(za(j1,j3)*za(j2,j3)*za(j4,j5)*za(j4,j6))
     &         *L0(-s(j2,j3),-t(j4,j5,j6))/t(j4,j5,j6)
     &       + za(j1,j2)**2
     &         *( za(j5,j2)*zb(j2,j1) + za(j5,j3)*zb(j3,j1) )
     &         *( za(j5,j4)*zb(j4,j1) + za(j5,j6)*zb(j6,j1) )
     &         /(2._dp*za(j1,j3)*za(j2,j3)*za(j4,j5)*za(j4,j6))
     &         *L1(-s(j2,j3),-t(j4,j5,j6))/t(j4,j5,j6)**2
c-----
      elseif(st=='q+qb-g+g-lb-l+') then
c-----(+-+--+)
      fQLlc =
     &       +  (za(j2,j1)*zb(j1,j6)+ za(j2,j3)*zb(j3,j6))
     &         *(za(j2,j4)*zb(j4,j6)+ za(j2,j5)*zb(j5,j6))
     &         /(t(j4,j5,j6)*za(j1,j3)*za(j2,j3)*zb(j4,j5)*zb(j4,j6))
     &         *Lsm1(-s(j1,j3),-t(j4,j5,j6),-s(j2,j3),-t(j4,j5,j6))
     &       - za(j2,j1)*zb(j1,j6)
     &         *(za(j2,j4)*zb(j4,j6) + za(j2,j5)*zb(j5,j6))
     &         /(za(j1,j3)*za(j2,j3)*zb(j4,j5)*zb(j4,j6))
     &         *L0(-s(j2,j3),-t(j4,j5,j6))/t(j4,j5,j6)
     &       - zb(j1,j2)**2*zb(j1,j6)**2
     &         /(2._dp*za(j1,j3)*za(j2,j3)*zb(j4,j5)*zb(j4,j6))
     &         *L1(-s(j2,j3),-t(j4,j5,j6))/t(j4,j5,j6)**2
c-----
      else
      write(6,*) 'unimplemented st',st
      stop
      endif
c-----
      return
      end

