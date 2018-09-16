      subroutine A5NLO_VH(j1,j2,j3,j4,j5,za,zb,A5LOm,A5NLOm)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5
      complex(dp):: A51_VH,A52_VH,A5NLOm,A5LOm,zab2
      real(dp):: s125

* As originally written, the functions A51, A52 correspond to
* 0 --> q_R(1)+qb_L(3)+g_R(2)+ebar_L(4)+e_R(5)
* with all RH couplings
* However we want it in our
* standard form
*       0--> qb_R(1)+q_L(2)++e_L(3)+ebar_R(4)+g_L(5)
* with all LH couplings

* so we have made the changes
*
*                    'q+g+qb-'   (A51)
*                   (1 ---> 2)
*                   (2 ---> 5)
*                   (3 ---> 1)
*                   (4 ---> 4)
*                   (5 ---> 3)

*                    'q+qb-g+'   (A52)
*                   (1 ---> 2)
*                   (2 ---> 1)
*                   (3 ---> 5)
*                   (4 ---> 4)
*                   (5 ---> 3)

*  and also exchanged za and zb.

C--- corresponds to (1V.1) times minus i, with the (A51) change
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

!      A5LOm=-zb(j1,j4)**2/(zb(j2,j5)*zb(j5,j1)*zb(j4,j3))
      s125=s(j1,j2)+s(j1,j5)+s(j2,j5)
      A5LOm=-zb(j1,j4)*zab2(j3,j2,j5,j1)/(zb(j2,j5)*zb(j5,j1)*s125)
      A5NLOm=A51_VH(j2,j5,j1,j4,j3,zb,za)
     &     +A52_VH(j2,j1,j5,j4,j3,zb,za)/xnsq

      return
      end
