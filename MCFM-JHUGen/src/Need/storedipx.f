      subroutine storedipx(msqx_st,msqvx_st,msqx_dip,msqx_dipv,
     & sub_st,subv_st,sub_dip,sub_dipv,n)
      implicit none
      include 'types.f'
c--- this routine transfers the information on the colour
c--- structure from a common block into separate arrays for
c--- each parton configuration
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: i,j,k,l,m,n
      real(dp):: 
     & msqx_st(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     & msqvx_st(0:2,-1:1,-1:1,-1:1,-1:1),
     & msqx_dip(36,0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     & msqx_dipv(36,0:2,-1:1,-1:1,-1:1,-1:1),
     & sub_st(4),sub_dip(36,4),subv_st,sub_dipv(36)
      

      do i=0,2
c        do j=-nf,nf
c        do k=-nf,nf
c        do l=-nf,nf
c        do m=-nf,nf
        do j=-2,2
        do k=-2,2
        do l=-2,2
        do m=-2,2
          msqx_dip(n,i,j,k,l,m)=msqx_st(i,j,k,l,m)
          if ((abs(j) < 2) .and. (abs(k) < 2)
     &  .and. (abs(l) < 2) .and. (abs(m) < 2)) then
            msqx_dipv(n,i,j,k,l,m)=msqvx_st(i,j,k,l,m)
          endif
        enddo
        enddo
        enddo
        enddo
      enddo
      
      do i=1,4
        sub_dip(n,i)=sub_st(i)
      enddo
      sub_dipv(n)=subv_st
      
      return
      end
