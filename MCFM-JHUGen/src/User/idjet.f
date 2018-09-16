      subroutine idjet(pjet,jetindex,numljets,numbjets,
     &     pljet,pbjet)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'jetlabel.f'
      real(dp):: pjet(mxpart,4),pljet(mxpart,4),pbjet(mxpart,4)
      integer:: numbjets,numljets,j,jetindex(mxpart)

      numljets=0
      numbjets=0

      do j=1,jets
         if (jetlabel(j) == 'pp' .or. jetlabel(j) == 'qj') then
            numljets=numljets+1
            pljet(numljets,:)=pjet(jetindex(j),:)
         elseif (jetlabel(j) == 'bq' .or. jetlabel(j) == 'ba'
     &      .or. jetlabel(j) == 'bb') then
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
