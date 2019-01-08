      subroutine storedip_mass(msq_dip,msq_dipv)
c--- this routine transfers the information on the colour
c--- structure from a common block into separate arrays for
c--- each parton configuration
      implicit none
      include 'constants.f'
      include 'msq_cs.f'
      include 'msqv_cs.f'
      integer i,j,k
      double precision 
     . msq_dip(0:2,-nf:nf,-nf:nf),msq_dipv(0:2,-nf:nf,-nf:nf)
      
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
