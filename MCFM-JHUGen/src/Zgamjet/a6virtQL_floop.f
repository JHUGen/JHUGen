
      function a6vQLfloop(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: a6vQLfloop
****************************************************
* virtual amplitude for 
* 0->q(p1)+qb(p2)+glu(p3)+gam(p4)+lb*p5)+l(p6)
* where one photon is coming from quark line
* fermion loop contributions
****************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      character*14 st
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: fQLfloop
c-----
      a6vQLfloop = fQLfloop(st,j1,j2,j3,j4,j5,j6,za,zb)
c-----
      return
      end


      function fQLfloop(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: fQLfloop
c-----finite part
      
      character*14 st
      integer:: j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      complex(dp):: L1
      real(dp):: t
c-----
      if(st=='q+qb-g+g+lb-l+') then
c-----(+-++-+)
      fQLfloop =
     &      + zb(j1,j3)*za(j2,j5)
     &        *( za(j5,j4)*zb(j4,j3) + za(j5,j6)*zb(j6,j3) )
     &        /(za(j4,j5)*za(j4,j6))
     &        *( L1(-s(j1,j2),-t(j4,j5,j6))/t(j4,j5,j6)**2
     &          -1._dp/(12._dp*t(j4,j5,j6)*mt**2) )
c-----
      elseif(st=='q+qb-g+g-lb-l+') then
c-----(+-+--+)
      fQLfloop =
     &      + zb(j1,j3)*zb(j3,j6)
     &        *( za(j2,j4)*zb(j4,j6) + za(j2,j5)*zb(j5,j6) )
     &        /(zb(j4,j5)*zb(j4,j6))
     &        *( L1(-s(j1,j2),-t(j4,j5,j6))/t(j4,j5,j6)**2
     &          -1._dp/(12._dp*t(j4,j5,j6)*mt**2) )
c-----
      else
      write(6,*) 'unimplemented st',st
      stop
      endif
c-----
      return
      end


