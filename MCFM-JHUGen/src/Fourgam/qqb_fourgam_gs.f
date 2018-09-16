      subroutine qqb_fourgam_gs(p,msq)
      implicit none
      include 'types.f'
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c---    q(p1) + q~(p2) --> gam(p3) + gam(p4) + gam(p5) +gam(p6) + g(p7)
c--- (and all crossings)
c---
c--- C. W, Sept 2014
       
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
      real(dp):: msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     & sub17_2(4),sub27_1(4),dummyv(-nf:nf,-nf:nf),dsubv,
     &  msq37_1(-nf:nf,-nf:nf), msq47_1(-nf:nf,-nf:nf),
     &  msq57_1(-nf:nf,-nf:nf), msq67_1(-nf:nf,-nf:nf),
     &  sub37_1,sub47_1,sub57_1,sub67_1
      external qqb_fourgam,qqb_trigam_g,donothing_gvec

      if (frag) then
         ndmax=6
      else
        ndmax=2
      endif
      
c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,7,2,sub17_2,dsubv,msq17_2,dummyv,
     & qqb_fourgam,donothing_gvec)
      call dips(2,p,2,7,1,sub27_1,dsubv,msq27_1,dummyv,
     & qqb_fourgam,donothing_gvec)
c--- extra photon dipoles in the case of fragmentation
      if (frag) then
         call dipsfrag(3,p,3,7,1,sub37_1,msq37_1,qqb_trigam_g)
         call dipsfrag(4,p,4,7,1,sub47_1,msq47_1,qqb_trigam_g)
         call dipsfrag(5,p,5,7,1,sub57_1,msq57_1,qqb_trigam_g)
         call dipsfrag(6,p,6,7,1,sub67_1,msq67_1,qqb_trigam_g)
         phot_dip(3:6)=.true.
      endif
      
      do j=-nf,nf
      do k=-nf,nf

      do nd=1,ndmax
      msq(nd,j,k)=0d0
      enddo

      if     ((j > 0) .and. (k < 0)
     &    .or.(j < 0) .and. (k > 0)) then
         msq(1,j,k)=2d0*cf*sub17_2(qq)*msq17_2(j,k)
         msq(2,j,k)=2d0*cf*sub27_1(qq)*msq27_1(j,k)
      elseif ((j .ne. 0) .and. (k == 0)) then
         msq(2,j,k)=2d0*tr*sub27_1(qg)*msq27_1(j,-j)
         if (frag) then
            msq(3,j,k)=Q(j)**2*msq37_1(j,k)*sub37_1/4d0
            msq(4,j,k)=Q(j)**2*msq47_1(j,k)*sub47_1/4d0
            msq(5,j,k)=Q(j)**2*msq57_1(j,k)*sub57_1/4d0
            msq(6,j,k)=Q(j)**2*msq67_1(j,k)*sub67_1/4d0
         endif
      elseif ((j == 0) .and. (k .ne. 0)) then
         msq(1,j,k)=2d0*tr*sub17_2(qg)*msq17_2(-k,k)
         if (frag) then
          msq(3,j,k)=Q(k)**2*msq37_1(j,k)*sub37_1/4d0
          msq(4,j,k)=Q(k)**2*msq47_1(j,k)*sub47_1/4d0
          msq(5,j,k)=Q(k)**2*msq57_1(j,k)*sub57_1/4d0
          msq(6,j,k)=Q(k)**2*msq67_1(j,k)*sub67_1/4d0
         endif
      endif

      enddo
      enddo

      return      
      end
