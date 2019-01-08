      subroutine qqb_trigam_gs(p,msq)
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c---    q(p1) + q~(p2) --> gam(p3) + gam(p4) + gam(p5) + g(p6)
c--- (and all crossings)
c---
c--- J. M. Campbell, March 2013
      implicit none 
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'ewcharge.f'
      include 'frag.f'
      include 'phot_dip.f'
      integer j,k,nd
      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
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
      msq(nd,j,k)=0d0
      enddo

      if     ((j .gt. 0) .and. (k .lt. 0)
     &    .or.(j .lt. 0) .and. (k .gt. 0)) then
         msq(1,j,k)=2d0*cf*sub16_2(qq)*msq16_2(j,k)
         msq(2,j,k)=2d0*cf*sub26_1(qq)*msq26_1(j,k)
      elseif ((j .ne. 0) .and. (k .eq. 0)) then
         msq(2,j,k)=2d0*tr*sub26_1(qg)*msq26_1(j,-j)
         if (frag) then
         msq(3,j,k)=Q(j)**2*msq36_1(j,k)*sub36_1/3d0
         msq(4,j,k)=Q(j)**2*msq46_1(j,k)*sub46_1/3d0
         msq(5,j,k)=Q(j)**2*msq56_1(j,k)*sub56_1/3d0
         endif
      elseif ((j .eq. 0) .and. (k .ne. 0)) then
         msq(1,j,k)=2d0*tr*sub16_2(qg)*msq16_2(-k,k)
         if (frag) then
         msq(3,j,k)=Q(k)**2*msq36_1(j,k)*sub36_1/3d0
         msq(4,j,k)=Q(k)**2*msq46_1(j,k)*sub46_1/3d0
         msq(5,j,k)=Q(k)**2*msq56_1(j,k)*sub56_1/3d0
         endif
      endif

      enddo
      enddo

      return      
      end
