      double precision function virt5(ip,za,zb)
      implicit none
************************************************************************ 
*     Author: R.K. Ellis                                               *
*     July, 1999.                                                      *
*   Given za and zb calculate the                                      *
*   the interference of the amplitude for the process                  *
*   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_L/R(5)                         *
*   at one loop with the corresponding lowest order amplitude          *
*   summed over the polarizations of the emitted gluon                 *
*   Virtual terms are in units of 
*   (as/4/pi) (4 pi)^ep Gamma(1+ep)*Gamma(1-ep)^2/Gamma(1-2*ep)
************************************************************************
      include 'constants.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      integer ip(5)
      double complex A5LOm,A5NLOm,A5LOp,A5NLOp

c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_L(5)
      call A5NLO(ip(1),ip(2),ip(3),ip(4),ip(5),za,zb,A5LOm,A5NLOm)
c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_R(5)
      call A5NLO(ip(2),ip(1),ip(4),ip(3),ip(5),zb,za,A5LOp,A5NLOp)

      virt5=
     . +ason2pi*(Dble(Dconjg(A5LOp)*A5NLOp)+Dble(Dconjg(A5LOm)*A5NLOm))

      return
      end

