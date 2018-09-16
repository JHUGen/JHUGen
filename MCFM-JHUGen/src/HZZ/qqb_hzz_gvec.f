      subroutine qqb_hzz_gvec(p,n,in,msq)
      implicit none
      include 'constants.f'
C  ip is the label of the emitting parton
C  kp is the label of the spectator parton
      integer j,k,in
      double precision msq(-nf:nf,-nf:nf),msqt(-nf:nf,-nf:nf)
      double precision n(4),nDn,p(mxpart,4)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2
      call qqb_hzz(p,msqt)

      msq(0,0)=-0.5d0*nDn*msqt(0,0)      

      return
      end

