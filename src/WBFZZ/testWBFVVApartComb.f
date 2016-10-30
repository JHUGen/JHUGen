      subroutine testWBFVVApartComb(j1,j2,j7,j8,partSwapOk)
      implicit none
      include 'constants.f'
      include 'plabel.f'
      integer j1,j2,j7,j8
      logical partSwapOk

      partSwapOk=.true.
      if (
     & (
     & (j1.eq.7 .or. j2.eq.7) .and. (
     &      (plabel(7) .eq. 'uq')
     & .or. (plabel(7) .eq. 'dq')
     & .or. (plabel(7) .eq. 'cq')
     & .or. (plabel(7) .eq. 'sq')
     & .or. (plabel(7) .eq. 'bq')
     & .or. (plabel(7) .eq. 'nl')
     & .or. (plabel(7) .eq. 'nm')
     & .or. (plabel(7) .eq. 'nt')
     & .or. (plabel(7) .eq. 'el')
     & .or. (plabel(7) .eq. 'ml')
     & .or. (plabel(7) .eq. 'tl')
     & )
     & ) .or. (
     & (j1.eq.8 .or. j2.eq.8) .and. (
     &      (plabel(8) .eq. 'uq')
     & .or. (plabel(8) .eq. 'dq')
     & .or. (plabel(8) .eq. 'cq')
     & .or. (plabel(8) .eq. 'sq')
     & .or. (plabel(8) .eq. 'bq')
     & .or. (plabel(8) .eq. 'nl')
     & .or. (plabel(8) .eq. 'nm')
     & .or. (plabel(8) .eq. 'nt')
     & .or. (plabel(8) .eq. 'el')
     & .or. (plabel(8) .eq. 'ml')
     & .or. (plabel(8) .eq. 'tl')
     & )
     & ) .or. (
     & (j7.eq.7 .or. j8.eq.7) .and. (
     &      (plabel(7) .eq. 'ua')
     & .or. (plabel(7) .eq. 'da')
     & .or. (plabel(7) .eq. 'ca')
     & .or. (plabel(7) .eq. 'sa')
     & .or. (plabel(7) .eq. 'ba')
     & .or. (plabel(7) .eq. 'na')
     & .or. (plabel(7) .eq. 'ea')
     & .or. (plabel(7) .eq. 'ma')
     & .or. (plabel(7) .eq. 'ta')
     & )
     & ) .or. (
     & (j7.eq.8 .or. j8.eq.8) .and. (
     &      (plabel(8) .eq. 'ua')
     & .or. (plabel(8) .eq. 'da')
     & .or. (plabel(8) .eq. 'ca')
     & .or. (plabel(8) .eq. 'sa')
     & .or. (plabel(8) .eq. 'ba')
     & .or. (plabel(8) .eq. 'na')
     & .or. (plabel(8) .eq. 'ea')
     & .or. (plabel(8) .eq. 'ma')
     & .or. (plabel(8) .eq. 'ta')
     & )
     & ) ) then
      partSwapOk = .false.
      endif

      return
      end
