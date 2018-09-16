      function threebee()
       implicit none
      include 'types.f'
      logical:: threebee

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'jetlabel.f'
      integer:: i,nbq

c--- note: this function returns true if there are three or more
c--- b-jets in the event

      nbq=0
      do i=1,jets
        if ((jetlabel(i)=='bq').or.(jetlabel(i)=='ba')) nbq=nbq+1
      enddo
      threebee=(nbq >= 3)
      return
      end

