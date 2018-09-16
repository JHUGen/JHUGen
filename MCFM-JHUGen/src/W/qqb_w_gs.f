      subroutine qqb_w_gs(p,msq)
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  l(p3)+a(p4)+g(p5)
c   positively charged W only
      implicit none 
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,nd

      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     . sub15_2(4),sub25_1(4),dummyv(-nf:nf,-nf:nf),dsubv
      external qqb_w,donothing_gvec

      ndmax=2

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,
     . qqb_w,donothing_gvec)
      call dips(2,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,
     . qqb_w,donothing_gvec)

      do j=-nf,nf
      do k=-nf,nf

      do nd=1,ndmax
      msq(nd,j,k)=0d0
      enddo

      if  ((j .eq. 0) .and. (k .eq. 0)) then
         goto 20
      elseif  ((j .gt. 0) .and. (k .lt. 0)
     .     .or.(j .lt. 0) .and. (k .gt. 0)) then
         msq(1,j,k)=2d0*cf*sub15_2(qq)*msq15_2(j,k)
         msq(2,j,k)=2d0*cf*sub25_1(qq)*msq25_1(j,k)
      elseif ((j .ne. 0) .and. (k .eq. 0)) then
         msq(2,j,k)=2d0*tr*sub25_1(qg)*(msq25_1(j,+1)+msq25_1(j,+2)
     &   +msq25_1(j,+3)+msq25_1(j,+4)+msq25_1(j,+5)
     &                                 +msq25_1(j,-1)+msq25_1(j,-2)
     &   +msq25_1(j,-3)+msq25_1(j,-4)+msq25_1(j,-5))
      elseif ((j .eq. 0) .and. (k .ne. 0)) then
         msq(1,j,k)=2d0*tr*sub15_2(qg)*(msq15_2(+1,k)+msq15_2(+2,k)
     &   +msq15_2(+3,k)+msq15_2(+4,k)+msq15_2(+5,k)
     &                                 +msq15_2(-1,k)+msq15_2(-2,k)
     &   +msq15_2(-3,k)+msq15_2(-4,k)+msq15_2(-5,k))
      endif
 20   continue

      enddo
      enddo

      return      
      end


