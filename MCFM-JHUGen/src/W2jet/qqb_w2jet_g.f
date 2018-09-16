      subroutine qqb_w2jet_g(p,msq)
************************************************************************
*     Author: J. M. Campbell                                           *
*     August, 2002.                                                    *
*                                                                      *
*     This is merely a wrapper routine to qqb_w(m/p)2jet_g             *
************************************************************************
      implicit none
      include 'constants.f'
      include 'nwz.f'
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf)

      if     (nwz .eq. +1) then
        call qqb_wp2jet_g(p,msq)
!      write(*,*) msq(-1,2),p(7,:)
      elseif (nwz .eq. -1) then
        call qqb_wm2jet_g(p,msq)
      else
        write(6,*) 'nwz not equal to +1 or -1 in'
        write(6,*) 'qqb_w2jet_g.f'
      endif
      return
      end

