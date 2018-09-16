****************************************************************
*   Color ordered amplitudes for:
*   0 -> q(-p1) + qb(-p4) + g(p2) + a(p3) + lb(p5) + l(p6)
****************************************************************
      subroutine zaj_a60h(j1,j2,j3,j4,j5,j6,za,zb,a60h)
      implicit none
      include 'types.f'
c-----the order of momentum in the argument
c-----(q,g,ph,qb,lb,l,za,zb,a60h)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a6treeg,a6treeQLlc
      complex(dp):: a60h(4,16)
      character*9 st1,st2,st3
      character*14 st4,st5 
      integer:: i,j,k
c-----helicity stamp
      st1='q+qb-g+g+'
      st2='q+qb-g+g-'
      st3='q+qb-g-g+'
      st4='q+qb-g+g+lb-l+'
      st5='q+qb-g+g-lb-l+'
c-----QQ contributions
      a60h(1,1) = a6treeg(st1,j1,j4,j2,j3,j5,j6,za,zb)
     &           +a6treeg(st1,j1,j4,j3,j2,j5,j6,za,zb)
      a60h(1,2) = a6treeg(st2,j1,j4,j2,j3,j5,j6,za,zb)
     &           +a6treeg(st3,j1,j4,j3,j2,j5,j6,za,zb)
      a60h(1,3) = a6treeg(st3,j1,j4,j2,j3,j5,j6,za,zb)
     &           +a6treeg(st2,j1,j4,j3,j2,j5,j6,za,zb)
      a60h(1,4) = a6treeg(st1,j4,j1,j2,j3,j6,j5,zb,za)
     &           +a6treeg(st1,j4,j1,j3,j2,j6,j5,zb,za)
      a60h(1,5) = a6treeg(st1,j1,j4,j2,j3,j6,j5,za,zb)
     &           +a6treeg(st1,j1,j4,j3,j2,j6,j5,za,zb)
      a60h(1,6) = a6treeg(st2,j1,j4,j2,j3,j6,j5,za,zb)
     &           +a6treeg(st3,j1,j4,j3,j2,j6,j5,za,zb)
      a60h(1,7) = a6treeg(st3,j1,j4,j2,j3,j6,j5,za,zb)
     &           +a6treeg(st2,j1,j4,j3,j2,j6,j5,za,zb)
      a60h(1,8) = a6treeg(st1,j4,j1,j2,j3,j5,j6,zb,za)
     &           +a6treeg(st1,j4,j1,j3,j2,j5,j6,zb,za)
      a60h(1,9) = a6treeg(st1,j1,j4,j2,j3,j5,j6,zb,za)
     &           +a6treeg(st1,j1,j4,j3,j2,j5,j6,zb,za)
      a60h(1,10) = a6treeg(st3,j4,j1,j2,j3,j6,j5,za,zb)
     &            +a6treeg(st2,j4,j1,j3,j2,j6,j5,za,zb)
      a60h(1,11) = a6treeg(st2,j4,j1,j2,j3,j6,j5,za,zb)
     &            +a6treeg(st3,j4,j1,j3,j2,j6,j5,za,zb)
      a60h(1,12) = a6treeg(st1,j4,j1,j2,j3,j6,j5,za,zb)
     &            +a6treeg(st1,j4,j1,j3,j2,j6,j5,za,zb)
      a60h(1,13) = a6treeg(st1,j1,j4,j2,j3,j6,j5,zb,za)
     &            +a6treeg(st1,j1,j4,j3,j2,j6,j5,zb,za)
      a60h(1,14) = a6treeg(st3,j4,j1,j2,j3,j5,j6,za,zb)
     &            +a6treeg(st2,j4,j1,j3,j2,j5,j6,za,zb)
      a60h(1,15) = a6treeg(st2,j4,j1,j2,j3,j5,j6,za,zb)
     &            +a6treeg(st3,j4,j1,j3,j2,j5,j6,za,zb)
      a60h(1,16) = a6treeg(st1,j4,j1,j2,j3,j5,j6,za,zb)
     &            +a6treeg(st1,j4,j1,j3,j2,j5,j6,za,zb)
c-----QL contributions
      a60h(3,1) = a6treeQLlc(st4,j1,j4,j2,j3,j5,j6,za,zb)
      a60h(3,2) = a6treeQLlc(st5,j1,j4,j2,j3,j5,j6,za,zb)
      a60h(3,3) = a6treeQLlc(st5,j4,j1,j2,j3,j6,j5,zb,za)
      a60h(3,4) = a6treeQLlc(st4,j4,j1,j2,j3,j6,j5,zb,za)
      a60h(3,5) = a6treeQLlc(st4,j1,j4,j2,j3,j6,j5,za,zb)
      a60h(3,6) = a6treeQLlc(st5,j1,j4,j2,j3,j6,j5,za,zb)
      a60h(3,7) = a6treeQLlc(st5,j4,j1,j2,j3,j5,j6,zb,za)
      a60h(3,8) = a6treeQLlc(st4,j4,j1,j2,j3,j5,j6,zb,za)
      a60h(3,9) = a6treeQLlc(st4,j1,j4,j2,j3,j5,j6,zb,za)
      a60h(3,10) = a6treeQLlc(st5,j1,j4,j2,j3,j5,j6,zb,za)
      a60h(3,11) = a6treeQLlc(st5,j4,j1,j2,j3,j6,j5,za,zb)
      a60h(3,12) = a6treeQLlc(st4,j4,j1,j2,j3,j6,j5,za,zb)
      a60h(3,13) = a6treeQLlc(st4,j1,j4,j2,j3,j6,j5,zb,za)
      a60h(3,14) = a6treeQLlc(st5,j1,j4,j2,j3,j6,j5,zb,za)
      a60h(3,15) = a6treeQLlc(st5,j4,j1,j2,j3,j5,j6,za,zb)
      a60h(3,16) = a6treeQLlc(st4,j4,j1,j2,j3,j5,j6,za,zb)
c-----donehere  
      return
      end
