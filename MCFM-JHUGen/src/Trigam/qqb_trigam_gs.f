      subroutine qqb_trigam_gs(p,msq)
      implicit none
      include 'types.f'
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c---    q(p1) + q~(p2) --> gam(p3) + gam(p4) + gam(p5) + g(p6)
c--- (and all crossings)
c---
c--- J. M. Campbell, March 2013
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'ewcharge.f'
      include 'frag.f'
      include 'phot_dip.f'
      integer:: j,k,nd
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp):: msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     & sub16_2(4),sub26_1(4),dummyv(-nf:nf,-nf:nf),dsubv,
     &  msq36_1(-nf:nf,-nf:nf), msq46_1(-nf:nf,-nf:nf),
     &  msq56_1(-nf:nf,-nf:nf),sub36_1,sub46_1,sub56_1
      external qqb_trigam,qqb_gmgmjt,donothing_gvec

      if (frag) then
        ndmax=5
      else
        ndmax=2
      endif
      
c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,6,2,sub16_2,dsubv,msq16_2,dummyv,
     & qqb_trigam,donothing_gvec)
      call dips(2,p,2,6,1,sub26_1,dsubv,msq26_1,dummyv,
     & qqb_trigam,donothing_gvec)
c--- extra photon dipoles in the case of fragmentation
      if (frag) then
         call dipsfrag(3,p,3,6,1,sub36_1,msq36_1,qqb_gmgmjt)
         call dipsfrag(4,p,4,6,1,sub46_1,msq46_1,qqb_gmgmjt)
         call dipsfrag(5,p,5,6,1,sub56_1,msq56_1,qqb_gmgmjt)
         phot_dip(3:5)=.true.
      endif
      
      do j=-nf,nf
      do k=-nf,nf

      do nd=1,ndmax
      msq(nd,j,k)=0._dp
      enddo

      if     ((j > 0) .and. (k < 0)
     &    .or.(j < 0) .and. (k > 0)) then
         msq(1,j,k)=2._dp*cf*sub16_2(qq)*msq16_2(j,k)
         msq(2,j,k)=2._dp*cf*sub26_1(qq)*msq26_1(j,k)
      elseif ((j .ne. 0) .and. (k == 0)) then
         msq(2,j,k)=2._dp*tr*sub26_1(qg)*msq26_1(j,-j)
         if (frag) then
         msq(3,j,k)=Q(j)**2*msq36_1(j,k)*sub36_1/3._dp
         msq(4,j,k)=Q(j)**2*msq46_1(j,k)*sub46_1/3._dp
         msq(5,j,k)=Q(j)**2*msq56_1(j,k)*sub56_1/3._dp
         endif
      elseif ((j == 0) .and. (k .ne. 0)) then
         msq(1,j,k)=2._dp*tr*sub16_2(qg)*msq16_2(-k,k)
         if (frag) then
         msq(3,j,k)=Q(k)**2*msq36_1(j,k)*sub36_1/3._dp
         msq(4,j,k)=Q(k)**2*msq46_1(j,k)*sub46_1/3._dp
         msq(5,j,k)=Q(k)**2*msq56_1(j,k)*sub56_1/3._dp
         endif
      endif

      enddo
      enddo

      return      
      end
