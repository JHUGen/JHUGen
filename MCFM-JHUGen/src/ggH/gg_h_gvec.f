      subroutine gg_h_gvec(p,n,in,msq)
      implicit none
      include 'constants.f'
C  in is the label of the momentum contracted with n
      integer j,k,in
      double precision msq(-nf:nf,-nf:nf),msqt(-nf:nf,-nf:nf)
      double precision n(4),nDn,p(mxpart,4)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2
      call gg_h(p,msqt)

      msq(0,0)=-0.5d0*nDn*msqt(0,0)      

      return
      end


