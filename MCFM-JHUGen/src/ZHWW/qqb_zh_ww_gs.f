      subroutine qqb_zh_ww_gs(p,msq)
      implicit none
      include 'types.f'
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  H  + Z +g(p7)
c                           |    |
c                           |    --> fermion(p3)+antifermion(p4)
c                           |
c                           ---> W^+(nu(p5),e^+(p6))W^-(e^-(p7),nub(p8))

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      integer:: j,k,nd

      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp):: msq19_2(-nf:nf,-nf:nf),msq29_1(-nf:nf,-nf:nf),
     & sub19_2(4),sub29_1(4),dummyv(-nf:nf,-nf:nf),dsubv
      external qqb_zh_ww,donothing_gvec

      ndmax=2

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,9,2,sub19_2,dsubv,msq19_2,dummyv,
     & qqb_zh_ww,donothing_gvec)
      call dips(2,p,2,9,1,sub29_1,dsubv,msq29_1,dummyv,
     & qqb_zh_ww,donothing_gvec)


      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
      msq(nd,j,k)=0._dp
      enddo

      if((j .ne. 0 .and. k .ne. 0 .and. j.ne.-k)
     &  .or. (j == 0 .and. k == 0)) goto 19

      if     (((j > 0) .and. (k < 0))
     &   .or. ((j < 0) .and. (k > 0))) then

      msq(1,j,k)=2._dp*cf*sub19_2(qq)*msq19_2(j,k)
      msq(2,j,k)=2._dp*cf*sub29_1(qq)*msq29_1(j,k)

      elseif ((j .ne. 0) .and. (k == 0)) then
      msq(1,j,k)=0._dp
      msq(2,j,k)=2._dp*tr*sub29_1(qg)*msq29_1(j,-j)

      elseif ((j == 0) .and. (k .ne. 0)) then
      msq(1,j,k)=2._dp*tr*sub19_2(qg)*msq19_2(-k,k)
      msq(2,j,k)=0._dp

      endif

   19 continue
      enddo
      enddo

      return
      end


