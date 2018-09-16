      function virt5_VHtop(ip,za,zb)
      implicit none
      include 'types.f'
      real(dp):: virt5_VHtop
************************************************************************ 
*     Author: C. Williams                                              *
*     July, 2015.                                                      *
*     I've defined My amplitudes as 1_L + 2_R +3_L +4_R                *
*     virt5_VH, does 1_R + 2_L + 3_L +4_R, so to use the same arrays I *
*     swap 2 and 1 w.r.t to the calling in that routine,               *
************************************************************************
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: ip(5)
      complex(dp):: A5LOm,A5NLOm,A5LOp,A5NLOp

c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_L(5)
      call A5NLO_VHtop(ip(2),ip(1),ip(3),ip(4),ip(5),za,zb,A5LOp,A5NLOp)
c   0--> qb_R(1)+q_L(2)+l_L(3)+a_R(4)+g_R(5)
      call A5NLO_VHtop(ip(1),ip(2),ip(4),ip(3),ip(5),zb,za,A5LOm,A5NLOm)

      virt5_VHtop= +real(conjg(A5LOp)*A5NLOp,dp)
     &             +real(conjg(A5LOm)*A5NLOm,dp)

      return
      end

