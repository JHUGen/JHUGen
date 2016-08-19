      subroutine swapjet(pjet,jetindex,i,j)
c--- swaps jets i..j in pjet
      implicit none
      include 'constants.f'
      include 'jetlabel.f'
      integer i,j,k,jetindex(mxpart)
      double precision pjet(mxpart,4),tmp
      character*2 chartmp
 
c--- escape if we're trying to swap the same jets
      if (i .eq. j) return

      do k=1,4
        tmp=pjet(i,k)
        pjet(i,k)=pjet(j,k)
        pjet(j,k)=tmp
      enddo
 
      chartmp=jetlabel(i)
      jetlabel(i)=jetlabel(j)
      jetlabel(j)=chartmp
      
      return
      end
