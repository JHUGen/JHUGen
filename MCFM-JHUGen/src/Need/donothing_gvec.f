      subroutine donothing_gvec(p,n,in,msq)
      implicit none
      include 'constants.f'
      integer j,k,in
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),n(4)
     
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      return
      end


