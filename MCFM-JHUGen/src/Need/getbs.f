      subroutine getbs(pjet,nbq,nba)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'jetlabel.f'
      integer:: i,nbq,nba
      real(dp):: pjet(mxpart,4)

c--- note: this function ASSUMES that there is at most one b-quark
c--- and one anti-b-quark, returning zero if there are less than this

      nbq=0
      nba=0

      do i=1,jets
        if (jetlabel(i) == 'bq') nbq=i+4
        if (jetlabel(i) == 'ba') nba=i+4
      enddo

      return
      end

