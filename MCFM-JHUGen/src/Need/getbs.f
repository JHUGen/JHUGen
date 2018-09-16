      subroutine getbs(pjet,nbq,nba)
      implicit none
      include 'constants.f'
      include 'jetlabel.f'
      integer i,nbq,nba
      double precision pjet(mxpart,4)
      
c--- note: this function ASSUMES that there is at most one b-quark
c--- and one anti-b-quark, returning zero if there are less than this

      nbq=0
      nba=0
     
      do i=1,jets
        if (jetlabel(i) .eq. 'bq') nbq=i+4
        if (jetlabel(i) .eq. 'ba') nba=i+4
      enddo

      return
      end
       
