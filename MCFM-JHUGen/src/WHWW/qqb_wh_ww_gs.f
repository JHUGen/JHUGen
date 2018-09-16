      subroutine qqb_wh_ww_gs(p,msq)
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c---for nwz=1
c     q(-p1)+qbar(-p2) -->  H  + W +g(p9)
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> W^+(nu(p5),e^+(p6))W^-(e^-(p7),nub(p8))
c---for nwz=-1
c     q(-p1)+qbar(-p2) -->  H  + W +g(p9)
c                           |    |
c                           |    --> e^-(p3)+nubar(p4)
c                           |
c                           ---> W^+(nu(p5),e^+(p6))W^-(e^-(p7),nub(p8))
      implicit none 
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      integer j,k,nd

      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision msq19_2(-nf:nf,-nf:nf),msq29_1(-nf:nf,-nf:nf),
     . sub19_2(4),sub29_1(4),dummyv(-nf:nf,-nf:nf),dsubv
      external qqb_wh_ww,donothing_gvec

      ndmax=2

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,9,2,sub19_2,dsubv,msq19_2,dummyv,
     . qqb_wh_ww,donothing_gvec)
      call dips(2,p,2,9,1,sub29_1,dsubv,msq29_1,dummyv,
     . qqb_wh_ww,donothing_gvec)

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
      msq(nd,j,k)=0d0
      enddo

      if  ((j .eq. 0) .and. (k .eq. 0)) then
         goto 20
      elseif  ((j .gt. 0) .and. (k .lt. 0)
     .     .or.(j .lt. 0) .and. (k .gt. 0)) then
         msq(1,j,k)=2d0*cf*sub19_2(qq)*msq19_2(j,k)
         msq(2,j,k)=2d0*cf*sub29_1(qq)*msq29_1(j,k)

      elseif ((j .ne. 0) .and. (k .eq. 0)) then
         msq(1,j,k)=0d0
         msq(2,j,k)=2d0*tr*sub29_1(qg)*(msq29_1(j,+1)+msq29_1(j,+2)
     &   +msq29_1(j,+3)+msq29_1(j,+4)+msq29_1(j,+5)
     &                                 +msq29_1(j,-1)+msq29_1(j,-2)
     &   +msq29_1(j,-3)+msq29_1(j,-4)+msq29_1(j,-5))

      elseif ((j .eq. 0) .and. (k .ne. 0)) then
         msq(1,j,k)=2d0*tr*sub19_2(qg)*(msq19_2(+1,k)+msq19_2(+2,k)
     &   +msq19_2(+3,k)+msq19_2(+4,k)+msq19_2(+5,k)
     &                                 +msq19_2(-1,k)+msq19_2(-2,k)
     &   +msq19_2(-3,k)+msq19_2(-4,k)+msq19_2(-5,k))
      msq(2,j,k)=0d0

      endif
 20   continue
      enddo
      enddo

      return      
      end


