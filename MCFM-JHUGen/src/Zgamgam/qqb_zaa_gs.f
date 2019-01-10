      subroutine qqb_zaa_gs(p,msq)
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  e-(p3)+e+(p4))+a(p5)+g(p6)
      implicit none 
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'frag.f'
      include 'ewcharge.f'
      include 'phot_dip.f'
      include 'ipsgen.f'
      integer j,k,nd
      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     . sub17_2(4),sub27_1(4),dummyv(-nf:nf,-nf:nf),dsubv
      double precision sub56_1,sub56_2,msq56_1(-nf:nf,-nf:nf),
     . msq56_2(-nf:nf,-nf:nf)
      double precision msq57_1(-nf:nf,-nf:nf),msq57_2(-nf:nf,-nf:nf),
     & msq67_1(-nf:nf,-nf:nf),msq67_2(-nf:nf,-nf:nf),
     & sub57_1,sub57_2,sub67_1,sub67_2
      external qqb_zaa,donothing_gvec
      external qqb_zaj

c---ipsgen=1: photons 5 and 6 radiated from initial quark line 
c---ipsgen=2: photons 5 and 6 radiated in Z decay 
c---ipsgen=3: photon 5 in Z decay, photon 6 radiated from initial quark line
c---ipsgen=4: photon 6 in Z decay, photon 5 radiated from initial quark line

      if ((frag) .and. (ipsgen .ne. 2)) then
         ndmax=4
      else
         ndmax=2
      endif

      do j=1,mxpart
         phot_dip(j)=.false.
      enddo

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,7,2,sub17_2,dsubv,msq17_2,dummyv,
     . qqb_zaa,donothing_gvec)
      call dips(2,p,2,7,1,sub27_1,dsubv,msq27_1,dummyv,
     . qqb_zaa,donothing_gvec)

c-----extra photon dipoles
      if (frag) then 
c--------gamma(5) dipoles: only if 5 is radiated from initial state
         if ((ipsgen .eq. 1) .or. (ipsgen .eq. 4)) then
         call dipsfrag(3,p,5,7,2,sub57_2,msq57_2,qqb_zaj)
         endif
c--------gamma(6) dipoles: only if 6 is radiated from initial state
         if ((ipsgen .eq. 1) .or. (ipsgen .eq. 3)) then
         call dipsfrag(4,p,6,7,2,sub67_2,msq67_2,qqb_zaj)
         endif
         do j=3,4
            phot_dip(j)=.true.
         enddo
      endif  


      do j=-nf,nf
      do k=-nf,nf

      do nd=1,ndmax
      msq(nd,j,k)=0d0
      enddo

      if  ((j .eq. 0) .and. (k .eq. 0)) then
         goto 20
      elseif  ((j .gt. 0) .and. (k .eq. -j)
     .     .or.(j .lt. 0) .and. (k .eq. -j)) then
         msq(1,j,k)=2d0*cf*sub17_2(qq)*msq17_2(j,k)
         msq(2,j,k)=2d0*cf*sub27_1(qq)*msq27_1(j,k)

      elseif ((j .ne. 0) .and. (k .eq. 0)) then
         msq(2,j,k)=2d0*tr*sub27_1(qg)*msq27_1(j,-j)
         if (frag) then 
            if ((ipsgen .eq. 1) .or. (ipsgen .eq. 4)) then
              msq(3,j,k)=Q(j)**2*msq57_2(j,k)*sub57_2*half
            endif
            if ((ipsgen .eq. 1) .or. (ipsgen .eq. 3)) then
              msq(4,j,k)=Q(j)**2*msq67_2(j,k)*sub67_2*half
            endif
         endif

      elseif ((j .eq. 0) .and. (k .ne. 0)) then
         msq(1,j,k)=2d0*tr*sub17_2(qg)*msq17_2(-k,k)         
         if (frag) then 
            if ((ipsgen .eq. 1) .or. (ipsgen .eq. 4)) then
              msq(3,j,k)=Q(k)**2*msq57_2(j,k)*sub57_2*half
            endif
            if ((ipsgen .eq. 1) .or. (ipsgen .eq. 3)) then
              msq(4,j,k)=Q(k)**2*msq67_2(j,k)*sub67_2*half
            endif
         endif
         
      endif
 20   continue

      enddo
      enddo

   
      return      
      end


