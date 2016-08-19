
      double complex function a6treeQLlc(st,j1,j2,j3,j4,j5,j6,za,zb)
c-----Tree amplitude for leading color contributions
c----- 0 -> q(p1) + qb(p2) + glu(p3) + gam(p4) + lb(p5) + l(p6)
c-----where the photon is coming from  the quark line
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      character*14 st
      integer j1,j2,j3,j4,j5,j6
      double precision t
c-----Multiplied by (-i), factor out sqrt(2) 
      if(st.eq.'q+qb-g+g+lb-l+') then
c-----(+-++-+)
      a6treeQLlc=
     .za(j2,j5)**2/(za(j1,j3)*za(j2,j3)*za(j4,j5)*za(j4,j6))
      elseif(st.eq.'q+qb-g+g-lb-l+') then
c-----(+-+--+)
      a6treeQLlc=
     .(za(j2,j4)*zb(j4,j6)+za(j2,j5)*zb(j5,j6))**2/
     .(t(j4,j5,j6)*za(j1,j3)*za(j2,j3)*zb(j4,j5)*zb(j4,j6))
      else
      write(6,*) 'unimplemented st',st
      stop
      endif
c-----done
      return
      end


      double complex function a6treeQLslc(st,j1,j2,j3,j4,j5,j6,za,zb)
c-----Tree amplitude for subleading color contributions
c----- 0 -> q(p1) + qb(p2) + glu(p3) + gam(p4) + lb(p5) + l(p6)
c-----where the photon is coming from  the lepton line
      implicit none
      integer j1,j2,j3,j4,j5,j6
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      character*14 st
      double precision t
c-----Multiplied by (-i), factor out sqrt(2) 
      if(st.eq.'q+qb-g+g+lb-l+') then
c-----(+-++-+)
      a6treeQLslc=
     .-za(j2,j5)**2/(za(j1,j3)*za(j2,j3)*za(j4,j5)*za(j4,j6))
      elseif(st.eq.'q+qb-g+g-lb-l+') then
c-----(+-+--+)
      a6treeQLslc= 
     .-(za(j2,j4)*zb(j4,j6)+za(j2,j5)*zb(j5,j6))**2/
     . (t(j4,j5,j6)*za(j1,j3)*za(j2,j3)*zb(j4,j5)*zb(j4,j6))
      else 
      write(6,*) 'unimplemented st',st
      stop
      endif
      return
      end



