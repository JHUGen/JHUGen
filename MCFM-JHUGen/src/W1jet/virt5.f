      function virt5(ip,za,zb)
      implicit none
      include 'types.f'
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
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'kprocess.f'
      include 'zprods_decl.f'
      real(dp):: virt5,virt5ax
      common/virt5ax/virt5ax
      integer:: ip(5)
      complex(dp):: A5LOm,A5NLOm,A5LOp,A5NLOp,A5axp,A5axm
!$omp threadprivate(/virt5ax/)

c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_L(5)
      call A5NLO(ip(1),ip(2),ip(3),ip(4),ip(5),za,zb,A5LOm,A5NLOm,A5axm)
c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_R(5)
      call A5NLO(ip(2),ip(1),ip(4),ip(3),ip(5),zb,za,A5LOp,A5NLOp,A5axp)

      virt5ax=zero
      if (kcase == kZ_1jet) then
      virt5ax=
     & +ason2pi*(real(conjg(A5LOp)*A5axp,dp)
     &          +real(conjg(A5LOm)*A5axm,dp))
      endif

      virt5=
     & +ason2pi*(real(conjg(A5LOp)*A5NLOp,dp)
     &          +real(conjg(A5LOm)*A5NLOm,dp))

      return
      end

