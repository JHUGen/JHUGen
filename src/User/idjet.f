      subroutine idjet(pjet,jetindex,numljets,numbjets,
     &     pljet,pbjet)
      implicit none
      include 'constants.f'
      include 'jetlabel.f'
      double precision pjet(mxpart,4),pljet(mxpart,4),pbjet(mxpart,4)
      integer numbjets,numljets,j,jetindex(mxpart)

      numljets=0
      numbjets=0

      do j=1,jets
         if (jetlabel(j) .eq. 'pp' .or. jetlabel(j) .eq. 'qj') then 
            numljets=numljets+1
            pljet(numljets,:)=pjet(jetindex(j),:)
         elseif (jetlabel(j) .eq. 'bq' .or. jetlabel(j) .eq. 'ba'
     &      .or. jetlabel(j) .eq. 'bb') then
            numbjets=numbjets+1
            pbjet(numbjets,:)=pjet(jetindex(j),:)
         else
            write(*,*) 'In idjet, something wrong in jetlabel', 
     &           jetlabel(j)
            stop
         endif
      enddo
      
      return
      end
