      subroutine wtgen(npart,q,wt)
      include 'constants.f'
      double precision q(mxpart,4),wt
      integer npart
      if (npart .eq. 5) then
      call wt4gen(q,wt)
      elseif (npart .eq. 3) then
      call wt2gen(q,wt)
      endif
      return
      end
