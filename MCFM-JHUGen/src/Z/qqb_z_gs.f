      subroutine qqb_z_gs(p,msq)
      implicit none
      include 'types.f'
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  Z^0((p3)+a(p4))+g(p5)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      integer:: j,k,nd

      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp):: msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     & sub15_2(4),sub25_1(4),dummyv(-nf:nf,-nf:nf),dsubv
      external qqb_z,donothing_gvec

      ndmax=2

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,
     & qqb_z,donothing_gvec)
      call dips(2,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,
     & qqb_z,donothing_gvec)


      do j=-nf,nf
      do k=-nf,nf

      do nd=1,ndmax
      msq(nd,j,k)=0._dp
      enddo

      if(((j .ne. 0) .and. (k .ne. 0) .and. (j.ne.-k))
     &  .or. ((j == 0) .and. (k == 0))) goto 19

      if     (((j > 0) .and. (k < 0))
     &   .or. ((j < 0) .and. (k > 0))) then

      msq(1,j,k)=2._dp*cf*sub15_2(qq)*msq15_2(j,k)
      msq(2,j,k)=2._dp*cf*sub25_1(qq)*msq25_1(j,k)

      elseif ((j .ne. 0) .and. (k == 0)) then
      msq(1,j,k)=0._dp
      msq(2,j,k)=2._dp*tr*sub25_1(qg)*msq25_1(j,-j)

      elseif ((j == 0) .and. (k .ne. 0)) then
      msq(1,j,k)=2._dp*tr*sub15_2(qg)*msq15_2(-k,k)
      msq(2,j,k)=0._dp

      endif

   19 continue
      enddo
      enddo

      return
      end


