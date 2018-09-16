      subroutine wtgen(npart,q,wt)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      real(dp):: q(mxpart,4),wt
      integer:: npart
      if (npart == 5) then
      call wt4gen(q,wt)
      elseif (npart == 3) then
      call wt2gen(q,wt)
      endif
      return
      end
