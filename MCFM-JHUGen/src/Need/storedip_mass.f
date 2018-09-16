      subroutine storedip_mass(msq_dip,msq_dipv)
      implicit none
      include 'types.f'
c--- this routine transfers the information on the colour
c--- structure from a common block into separate arrays for
c--- each parton configuration
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'msq_cs.f'
      include 'msqv_cs.f'
      integer:: i,j,k
      real(dp):: 
     & msq_dip(0:2,-nf:nf,-nf:nf),msq_dipv(0:2,-nf:nf,-nf:nf)
      
      do i=0,2
        do j=-nf,nf
        do k=-nf,nf
          msq_dip(i,j,k)=msq_cs(i,j,k)
          msq_dipv(i,j,k)=msqv_cs(i,j,k)
        enddo
        enddo
      enddo
      
      return
      end
