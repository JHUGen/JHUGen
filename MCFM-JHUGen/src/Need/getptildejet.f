      subroutine getptildejet(nd,pjet)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'npart.f'
      include 'ptilde.f'
      integer:: nd,i,j
      real(dp):: pjet(mxpart,4)

      do j=1,4
        do i=1,npart+2
        pjet(i,j)=ptildejet(nd,i,j)
        enddo
      enddo
      pjet(npart+3,:)=0._dp   ! ensure next entry is zero

      return
      end

