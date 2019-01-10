      subroutine pvclearcache
      implicit none
      include 'TRclear.f'
      include 'pvRespectmaxcindex.f'
      include 'TRbadpoint.f'
      integer j
      do j=1,5
      clear(j)=.true.
      enddo
      pvRespectmaxcindex=.true.
      pvbadpoint=.false.
      return
      end
