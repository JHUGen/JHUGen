      subroutine qqb_hzz_gs(p,msq)
      implicit none
      include 'types.f'
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c----for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  Z- (mu^-(p5)+mu^+(p6)) + Z (e^-(p3)+e^+(p4))
c     +g(p7)
c---

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      integer:: j,k,nd

      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp):: msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     & msq17_2v(-nf:nf,-nf:nf),msq27_1v(-nf:nf,-nf:nf),
     & sub17_2(4),sub27_1(4),sub17_2v,sub27_1v
      external qqb_hzz,qqb_hzz_gvec

      ndmax=2

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,7,2,sub17_2,sub17_2v,msq17_2,msq17_2v,
     & qqb_hzz,qqb_hzz_gvec)
      call dips(2,p,2,7,1,sub27_1,sub27_1v,msq27_1,msq27_1v,
     & qqb_hzz,qqb_hzz_gvec)

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
      msq(nd,j,k)=0._dp
      enddo


      if     ((j .ne. 0) .and. (k == 0)) then
         msq(1,j,k)=two*cf
     &   *(msq17_2(0,0)*sub17_2(gq)+msq17_2v(0,0)*sub17_2v)
      elseif ((j == 0) .and. (k .ne. 0)) then
         msq(2,j,k)=two*cf
     &   *(msq27_1(0,0)*sub27_1(gq)+msq27_1v(0,0)*sub27_1v)
      elseif ((j == 0) .and. (k == 0)) then
         msq(1,j,k)=two*xn
     &   *(msq17_2(j,k)*sub17_2(gg)+msq17_2v(j,k)*sub17_2v)
         msq(2,j,k)=two*xn
     &   *(msq27_1(j,k)*sub27_1(gg)+msq27_1v(j,k)*sub27_1v)
      endif


      enddo
      enddo

      return
      end


