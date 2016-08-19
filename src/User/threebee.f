      logical function threebee()
      implicit none
      include 'constants.f'
      include 'jetlabel.f'
      integer i,nbq
      
c--- note: this function returns true if there are three or more 
c--- b-jets in the event

      nbq=0
      do i=1,jets
        if ((jetlabel(i).eq.'bq').or.(jetlabel(i).eq.'ba')) nbq=nbq+1
      enddo
      threebee=(nbq .ge. 3)
      return
      end
       
