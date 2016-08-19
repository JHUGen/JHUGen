      subroutine qqb_w2jet_gs(p,msqc)
************************************************************************
*     Author: J. M. Campbell                                           *
*     August, 2002.                                                    *
*                                                                      *
*     This is merely a wrapper routine to qqb_w(m/p)2jet_gs            *
************************************************************************
      implicit none 
      include 'constants.f'
      include 'nwz.f'
      include 'ptilde.f'
      double precision p(mxpart,4),msqc(maxd,-nf:nf,-nf:nf)

      if     (nwz .eq. +1) then
        call qqb_wp2jet_gs(p,msqc)
      elseif (nwz .eq. -1) then
        call qqb_wm2jet_gs(p,msqc)
      else
        write(6,*) 'nwz not equal to +1 or -1 in'
        write(6,*) 'qqb_w2jet_gs.f'
      endif
      
      return
      end
      
