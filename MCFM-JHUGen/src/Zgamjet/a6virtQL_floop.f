
      double complex function a6vQLfloop(st,j1,j2,j3,j4,j5,j6,za,zb)
****************************************************
* virtual amplitude for 
* 0->q(p1)+qb(p2)+glu(p3)+gam(p4)+lb*p5)+l(p6)
* where one photon is coming from quark line
* fermion loop contributions
****************************************************
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      character*14 st
      integer j1,j2,j3,j4,j5,j6
      double complex fQLfloop
c-----
      a6vQLfloop = fQLfloop(st,j1,j2,j3,j4,j5,j6,za,zb)
c-----
      return
      end


      double complex function fQLfloop(st,j1,j2,j3,j4,j5,j6,za,zb)
c-----finite part
      implicit none
      character*14 st
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      double complex L1
      double precision t
c-----
      if(st.eq.'q+qb-g+g+lb-l+') then
c-----(+-++-+)
      fQLfloop =
     .      + zb(j1,j3)*za(j2,j5)
     .        *( za(j5,j4)*zb(j4,j3) + za(j5,j6)*zb(j6,j3) )
     .        /(za(j4,j5)*za(j4,j6))
     .        *( L1(-s(j1,j2),-t(j4,j5,j6))/t(j4,j5,j6)**2
     .          -1D0/(12D0*t(j4,j5,j6)*mt**2) )
c-----
      elseif(st.eq.'q+qb-g+g-lb-l+') then
c-----(+-+--+)
      fQLfloop =
     .      + zb(j1,j3)*zb(j3,j6)
     .        *( za(j2,j4)*zb(j4,j6) + za(j2,j5)*zb(j5,j6) )
     .        /(zb(j4,j5)*zb(j4,j6))
     .        *( L1(-s(j1,j2),-t(j4,j5,j6))/t(j4,j5,j6)**2
     .          -1D0/(12D0*t(j4,j5,j6)*mt**2) )
c-----
      else
      write(6,*) 'unimplemented st',st
      stop
      endif
c-----
      return
      end


