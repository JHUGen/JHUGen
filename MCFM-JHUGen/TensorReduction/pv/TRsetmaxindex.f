      subroutine TRsetmaxindex(cindex,dindex,eindex)
c--- sets the maximum number of tensor indices calculated in the
c--- routines pvCfill (cindex), pvDfill (dindex) and pvEfill (eindex);
c--- should be set on a per-calculation basis
      implicit none
      include 'TRmaxindex.f'
      integer cindex,dindex,eindex
      maxcindex=cindex
      maxdindex=dindex
      maxeindex=eindex
      return
      end
