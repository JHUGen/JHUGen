      subroutine a6treemass(k1,k2,k3,k4,k5,k6,mq,
     & a6treemm,a6treemp,a6treepm,a6treepp)
      implicit none
      include 'types.f'
c--- Routine to compute the tree level amplitude A6tree for the Wbb
c--- process, where the mass of the b-quark is kept non-zero.
c--- The labels on this routine refer to the momenta in "mom" (passed
c--- via common block), in which the massive momenta have been made massless
c--- a la Rodrigo and appear in positions 2 and 3
c---
c---     0 -> q(k1) + qb(k4) + W(->e(k6)+nubar(k5)) + Q(k3) + Qbar(k2)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
c      include 'momwbbm.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: k1,k2,k3,k4,k5,k6
      complex(dp):: a6treemm,a6treemp,a6treepm,a6treepp,
     & zba5243,zba5342,zab6123,zab6132
      real(dp):: s123,s234,mq
c      zab2(k1,k2,k3,k4)=za(k1,k2)*zb(k2,k4)+za(k1,k3)*zb(k3,k4)
c      zba2(k1,k2,k3,k4)=zb(k1,k2)*za(k2,k4)+zb(k1,k3)*za(k3,k4)
      zba5243=zb(k5,k2)*za(k2,k3)+zb(k5,k4)*za(k4,k3)
      zba5342=zb(k5,k3)*za(k3,k2)+zb(k5,k4)*za(k4,k2)
      zab6132=za(k6,k1)*zb(k1,k2)+za(k6,k3)*zb(k3,k2)
      zab6123=za(k6,k1)*zb(k1,k3)+za(k6,k2)*zb(k2,k3)
      s123=s(k1,k2)+s(k2,k3)+s(k3,k1)
      s234=s(k2,k3)+s(k3,k4)+s(k4,k2)
      
c      propagators s(k2,k3) and s(k5,k6) removed
      a6treemp=za(k1,k3)*zb(k4,k5)*zab6132/s123
     &        -zb(k4,k2)*za(k1,k6)*zba5243/s234
      a6treepm=za(k1,k2)*zb(k4,k5)*zab6123/s123
     &        -zb(k4,k3)*za(k1,k6)*zba5342/s234

c -original expression
c      a6tree(1,1)=mq/zb(k2,k3)*(
c     &  zb(k4,k5)*(zab6123*za(k3,k1)-zab6132*za(k2,k1))/s123
c     & +za(k1,k6)*(zba5342*zb(k2,k4)-zba5243*zb(k3,k4))/s234)

      a6treemm=two*mq/zb(k2,k3)*(
     &  zb(k4,k5)*zab6123*za(k3,k1)/s123
     & -za(k1,k6)*zba5243*zb(k3,k4)/s234)

      a6treepp=a6treemm*zb(k3,k2)/za(k3,k2)

c      add propagators s(k2,k3) and s(k5,k6)
      a6treemm=a6treemm/(s(k2,k3)*s(k5,k6))
      a6treemp=a6treemp/(s(k2,k3)*s(k5,k6))
      a6treepm=a6treepm/(s(k2,k3)*s(k5,k6))
      a6treepp=a6treepp/(s(k2,k3)*s(k5,k6))

      return
      end
      
      
