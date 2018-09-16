      double complex function atree(st,j1,j2,j3,j4,j5,j6,za,zb) 
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
************************************************************************
c---- the basic process is atreepp which is the amplitude for
c---- q-(p4)+Q(-p2)+l-(-p5) ---> q+(p1)+Q-(p3)+l+(p6)
c---- All outgoing particles are right-handed, 
c---- except Q(p3) which is left-handed
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4,j5,j6
      character*2 st
      double precision t
      double complex tree
c-----Note all tree amplitudes divided by im wrt to BDKW paper
c-----Statement function for pp (others derived by crossing)
c-----This is an implementation of BDKW Eq.(3.10)
      tree(j1,j2,j3,j4,j5,j6)=
     . +zb(j1,j2)*za(j5,j4)*(za(j3,j1)*zb(j1,j6)+za(j3,j2)*zb(j2,j6))
     .      /(s(j2,j3)*s(j5,j6)*t(j1,j2,j3))
     . +za(j3,j4)*zb(j6,j1)*(za(j5,j3)*zb(j3,j2)+za(j5,j4)*zb(j4,j2))
     .       /(s(j2,j3)*s(j5,j6)*t(j2,j3,j4))

      if (st .eq. 'pp') then
         atree=+tree(j1,j2,j3,j4,j5,j6)
      elseif (st .eq. 'pm') then
         atree=-tree(j1,j3,j2,j4,j5,j6)
      elseif (st .eq. 'sl') then
c---this is not a real tree amplitude but it is an auxiliary quantity.
         atree=-tree(j3,j1,j2,j4,j5,j6)
      endif

      return 
      end
