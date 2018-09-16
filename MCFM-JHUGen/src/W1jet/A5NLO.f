      subroutine A5NLO(j1,j2,j3,j4,j5,za,zb,A5LOm,A5NLOm,A5ax)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'kprocess.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5
      complex(dp):: A51,A52,A5NLOm,A5LOm,A5ax,zab2,F1anom
!      complex(dp):: L1
      real(dp):: s125,s12,mt2,musq

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


      A5ax=czip

!      A5LOm=-zb(j1,j4)**2/(zb(j2,j5)*zb(j5,j1)*zb(j4,j3))
      s125=s(j1,j2)+s(j1,j5)+s(j2,j5)
      A5LOm=-zb(j1,j4)*zab2(j3,j2,j5,j1)/(zb(j2,j5)*zb(j5,j1)*s125)
      A5NLOm=A51(j2,j5,j1,j4,j3,zb,za)+A52(j2,j1,j5,j4,j3,zb,za)/xnsq

      if (kcase .eq. kZ_1jet) then
      s12=s(j1,j2)
      mt2=mt**2
!     Musq is irrelevant,set it to some value
      musq=abs(s125)
      A5ax=-za(j2,j5)*za(j3,j5)*zb(j1,j4)/s125
     & *(F1anom(s12,s125,mt2,musq)-F1anom(s12,s125,zero,musq))

!      write(6,*) 'A5ax exact ',A5ax
!      A5ax=za(j2,j5)*za(j3,j5)*zb(j1,j4)/s125
!     & *(-two*s125)*(L1(-s12,-s125)/s125**2-one/(twelve*s125*mt2))
!      write(6,*) 'A5ax approx',A5ax
!      pause

      endif

      return
      end
