      subroutine qqb_wbjet_gs(p,msqc)
************************************************************************
*     Author: J. M. Campbell                                           *
*     January, 2004.                                                   *
*                                                                      *
*     This is merely a wrapper routine to qqb_w(m/p)bjet_gs            *
************************************************************************
      implicit none 
      include 'constants.f'
      include 'nwz.f'
      include 'ptilde.f'
      double precision p(mxpart,4),msqc(maxd,-nf:nf,-nf:nf)

      if     (nwz .eq. +1) then
        call qqb_wpbjet_gs(p,msqc)
      elseif (nwz .eq. -1) then
        call qqb_wmbjet_gs(p,msqc)
      else
        write(6,*) 'nwz not equal to +1 or -1 in'
        write(6,*) 'qqb_wbjet_gs.f'
      endif
      
      return
      end
      
