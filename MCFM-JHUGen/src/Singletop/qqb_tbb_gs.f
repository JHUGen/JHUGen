      subroutine qqb_tbb_gs(p,msq)
      implicit none
      include 'types.f'
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  t(l(p3)+a(p4)+c(p5))+b(p6)+g(p7)
c   positively charged W only

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'nwz.f'
      include 'qqgg.f'
      include 'stopbmass.f'
      include 'breit.f'
      integer:: j,k,ib,nd,isub

      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp)::
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     & msq35_4(-nf:nf,-nf:nf),msq45_3(-nf:nf,-nf:nf),
     & msq14_2(-nf:nf,-nf:nf),msq24_1(-nf:nf,-nf:nf),
     & msq15_4(-nf:nf,-nf:nf),msq25_3(-nf:nf,-nf:nf),
     & msq45_1(-nf:nf,-nf:nf),msq35_2(-nf:nf,-nf:nf),
     & msq25_4(-nf:nf,-nf:nf),msq15_3(-nf:nf,-nf:nf),
     & msq45_2(-nf:nf,-nf:nf),msq35_1(-nf:nf,-nf:nf),
     & sub15_2(4),sub25_1(4),
     & sub35_4(4),sub45_3(4),
     & sub14_2(4),sub24_1(4),
     & sub15_4(4),sub25_3(4),sub45_1(4),sub35_2(4),
     & sub15_3(4),sub25_4(4),sub45_2(4),sub35_1(4),
     & dummyv(-nf:nf,-nf:nf),dsubv
      real(dp):: oldmass2
      external qqb_tbb,bq_tpq,donothing_gvec

      common/isub/isub

c--- section for two b-quarks in the final state (s-channel, plus
c---  some t-channel g-q diagrams)
      if (isub == 2) then

      ndmax=4

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies

      call dips_mass(1,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,
     & qqb_tbb,donothing_gvec)
      call dips_mass(2,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,
     & qqb_tbb,donothing_gvec)
      qqproc=.true.
      call dips_mass(3,p,3,5,4,sub35_4,dsubv,msq35_4,dummyv,
     & qqb_tbb,donothing_gvec)
      oldmass2=mass2
      mass2=0._dp
      call dips_mass(4,p,4,5,3,sub45_3,dsubv,msq45_3,dummyv,
     & qqb_tbb,donothing_gvec)
      mass2=oldmass2
c      call dips_mass(5,p,1,4,2,sub14_2,dsubv,msq14_2,dummyv,
c     & bq_tpq,donothing_gvec)
c      call dips_mass(6,p,2,4,1,sub24_1,dsubv,msq24_1,dummyv,
c     & bq_tpq,donothing_gvec)

      do j=-nf,nf
      do k=-nf,nf

      do nd=1,ndmax
      msq(nd,j,k)=0._dp
      enddo

      if  ((j == 0) .and. (k == 0)) then
         goto 20
      elseif  ((j > 0) .and. (k < 0)
     &     .or.(j < 0) .and. (k > 0)) then
         msq(1,j,k)=2._dp*cf*sub15_2(qq)*msq15_2(j,k)
         msq(2,j,k)=2._dp*cf*sub25_1(qq)*msq25_1(j,k)
         msq(3,j,k)=2._dp*cf*sub35_4(qq)*msq35_4(j,k)
         msq(4,j,k)=2._dp*cf*sub45_3(qq)*msq45_3(j,k)
      elseif ((j .ne. 0) .and. (k == 0)) then
         msq(2,j,k)=2._dp*tr*sub25_1(qg)*(msq25_1(j,+1)+msq25_1(j,+2)
     &   +msq25_1(j,+3)+msq25_1(j,+4)+msq25_1(j,+5)
     &                                 +msq25_1(j,-1)+msq25_1(j,-2)
     &   +msq25_1(j,-3)+msq25_1(j,-4)+msq25_1(j,-5))
c         msq(6,j,k)=2._dp*tr*sub24_1(qg)*msq24_1(j,5)
      elseif ((j == 0) .and. (k .ne. 0)) then
         msq(1,j,k)=2._dp*tr*sub15_2(qg)*(msq15_2(+1,k)+msq15_2(+2,k)
     &   +msq15_2(+3,k)+msq15_2(+4,k)+msq15_2(+5,k)
     &                                 +msq15_2(-1,k)+msq15_2(-2,k)
     &   +msq15_2(-3,k)+msq15_2(-4,k)+msq15_2(-5,k))
c         msq(5,j,k)=2._dp*tr*sub14_2(qg)*msq14_2(5,k)
      endif
 20   continue

      enddo
      enddo

      endif

c--- section for one b-quark in the final state (only t-channel)
      if (isub == 1) then
      ndmax=12

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies

      qqproc=.true.
      oldmass2=mass2
      mass2=0._dp
      call dips_mass(1,p,1,5,4,sub15_4,dsubv,msq15_4,dummyv,
     & bq_tpq,donothing_gvec)
      call dips_mass(2,p,4,5,1,sub45_1,dsubv,msq45_1,dummyv,
     & bq_tpq,donothing_gvec)
      mass2=oldmass2
      call dips_mass(3,p,2,5,3,sub25_3,dsubv,msq25_3,dummyv,
     & bq_tpq,donothing_gvec)
      call dips_mass(4,p,3,5,2,sub35_2,dsubv,msq35_2,dummyv,
     & bq_tpq,donothing_gvec)
      oldmass2=mass2
      mass2=0._dp
      call dips_mass(5,p,2,5,4,sub25_4,dsubv,msq25_4,dummyv,
     & bq_tpq,donothing_gvec)
      call dips_mass(6,p,4,5,2,sub45_2,dsubv,msq45_2,dummyv,
     & bq_tpq,donothing_gvec)
      mass2=oldmass2
      call dips_mass(7,p,1,5,3,sub15_3,dsubv,msq15_3,dummyv,
     & bq_tpq,donothing_gvec)
      call dips_mass(8,p,3,5,1,sub35_1,dsubv,msq35_1,dummyv,
     & bq_tpq,donothing_gvec)

      call dips_mass(9,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,
     & bq_tpq,donothing_gvec)
      call dips_mass(10,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,
     & bq_tpq,donothing_gvec)

      call dips_mass(11,p,1,4,2,sub14_2,dsubv,msq14_2,dummyv,
     & bq_tpq,donothing_gvec)
      call dips_mass(12,p,2,4,1,sub24_1,dsubv,msq24_1,dummyv,
     & bq_tpq,donothing_gvec)

c--- for nwz=+1, initial state is b, for nwz=-1 it is b~
      ib=5*nwz

      do j=-nf,nf
      do k=-nf,nf

      do nd=1,ndmax
      msq(nd,j,k)=0._dp
      enddo

      if  ((j == 0) .and. (k == 0)) then
         goto 21
      elseif ((j == ib) .and. (k .ne. 0)) then
c--- checked 15,25,45 singularities, soft singularity looks ok
        msq(5,j,k)=2._dp*cf*sub25_4(qq)*msq25_4(j,k)
        msq(6,j,k)=2._dp*cf*sub45_2(qq)*msq45_2(j,k)
        msq(7,j,k)=2._dp*cf*sub15_3(qq)*msq15_3(j,k)
        msq(8,j,k)=2._dp*cf*sub35_1(qq)*msq35_1(j,k)
      elseif ((j .ne. 0) .and. (k == ib)) then
c--- checked 15,25,45 singularities, soft singularity looks ok
        msq(1,j,k)=2._dp*cf*sub15_4(qq)*msq15_4(j,k)
        msq(2,j,k)=2._dp*cf*sub45_1(qq)*msq45_1(j,k)
        msq(3,j,k)=2._dp*cf*sub25_3(qq)*msq25_3(j,k)
        msq(4,j,k)=2._dp*cf*sub35_2(qq)*msq35_2(j,k)
c--- checked
      elseif ((j == ib) .and. (k == 0)) then
         msq(10,j,k)=2._dp*tr*sub25_1(qg)*(msq25_1(j,+1)+msq25_1(j,+2)
     &   +msq25_1(j,+3)+msq25_1(j,+4)+msq25_1(j,+5))
         msq(12,j,k)=2._dp*tr*sub24_1(qg)*(msq24_1(j,-1)+msq24_1(j,-2)
     &   +msq24_1(j,-3)+msq24_1(j,-4)+msq24_1(j,-5))
      elseif ((j .ne. 0) .and. (k == 0)) then
         if (masslessb) then
         msq(10,j,k)=2._dp*tr*sub25_1(qg)*msq25_1(j,ib)
         endif
      elseif ((j == 0) .and. (k == ib)) then
         msq(9,j,k)=2._dp*tr*sub15_2(qg)*(msq15_2(+1,k)+msq15_2(+2,k)
     &   +msq15_2(+3,k)+msq15_2(+4,k)+msq15_2(+5,k))
         msq(11,j,k)=2._dp*tr*sub14_2(qg)*(msq14_2(-1,k)+msq14_2(-2,k)
     &   +msq14_2(-3,k)+msq14_2(-4,k)+msq14_2(-5,k))
      elseif ((j == 0) .and. (k .ne. 0)) then
         if (masslessb) then
         msq(9,j,k)=2._dp*tr*sub15_2(qg)*msq15_2(ib,k)
         endif
      endif

 21   continue

      enddo
      enddo

      endif

      return
      end


