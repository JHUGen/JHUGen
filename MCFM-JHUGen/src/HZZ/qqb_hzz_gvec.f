      subroutine qqb_hzz_gvec(p,n,in,msq)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
C  ip is the label of the emitting parton
C  kp is the label of the spectator parton
      integer:: j,k,in
      real(dp):: msq(-nf:nf,-nf:nf),msqt(-nf:nf,-nf:nf)
      real(dp):: n(4),nDn,p(mxpart,4)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2
      call qqb_hzz(p,msqt)

      msq(0,0)=-0.5_dp*nDn*msqt(0,0)

      return
      end

