      subroutine A5NLO(j1,j2,j3,j4,j5,za,zb,A5LOm,A5NLOm)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5
      double complex A51,A52,A5NLOm,A5LOm

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
      A5LOm=-zb(j1,j4)**2/(zb(j2,j5)*zb(j5,j1)*zb(j4,j3))
      A5NLOm=A51(j2,j5,j1,j4,j3,zb,za)+A52(j2,j1,j5,j4,j3,zb,za)/xnsq

      return
      end
