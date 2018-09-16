      subroutine qqb_w2jet_gs(p,msqc)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J. M. Campbell                                           *
*     August, 2002.                                                    *
*                                                                      *
*     This is merely a wrapper routine to qqb_w(m/p)2jet_gs            *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'nwz.f'
      include 'ptilde.f'
      real(dp):: p(mxpart,4),msqc(maxd,-nf:nf,-nf:nf)

      if     (nwz == +1) then
        call qqb_wp2jet_gs(p,msqc)
      elseif (nwz == -1) then
        call qqb_wm2jet_gs(p,msqc)
      else
        write(6,*) 'nwz not equal to +1 or -1 in'
        write(6,*) 'qqb_w2jet_gs.f'
      endif

      return
      end

