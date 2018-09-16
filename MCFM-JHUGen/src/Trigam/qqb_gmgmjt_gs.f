      subroutine qqb_gmgmjt_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*    Author: J.M. Campbell                                             *
*    March, 2013.                                                      *
*    Matrix element SUBTRACTION squared averaged over initial colors   *
*    and spins                                                         *
c     q(-p1)+qbar(-p2) --> gam(p3) + gam(p4) + parton(p5) + parton(p6) *
************************************************************************


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'nflav.f'
      include 'qqgg.f'
      include 'phot_dip.f'
      include 'frag.f'
      include 'ewcharge.f'
      include 'msqbits.f'
      integer:: j,k,nd
c --- remember: nd will count the dipoles

      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf),facgg,
     & pdip(mxpart,4),pswap(mxpart,4)
      real(dp)::
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     & msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     & msq15_6(-nf:nf,-nf:nf),msq26_5(-nf:nf,-nf:nf),
     & msq16_5(-nf:nf,-nf:nf),msq25_6(-nf:nf,-nf:nf),
     & msq56_1v(-nf:nf,-nf:nf),msq56_2v(-nf:nf,-nf:nf),
     & msq26_5v(-nf:nf,-nf:nf),msq26_1v(-nf:nf,-nf:nf),
     & msq15_6v(-nf:nf,-nf:nf),msq16_2v(-nf:nf,-nf:nf),
     & msq16_5v(-nf:nf,-nf:nf),msq25_6v(-nf:nf,-nf:nf),
     & msq25_1v(-nf:nf,-nf:nf),
     & msq15_2v(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),
     & sub15_2(4),sub25_1(4),sub16_2(4),sub26_1(4),
     & sub15_6(4),sub16_5(4),sub25_6(4),sub26_5(4),
     & sub56_1(4),sub56_2(4),sub56_1v,sub56_2v,
     & sub26_5v,sub25_1v,sub26_1v,sub16_5v,sub16_2v,sub15_2v,sub15_6v,
     & sub25_6v,
     & msq35_1(-nf:nf,-nf:nf),msq36_1(-nf:nf,-nf:nf),
     & msq45_1(-nf:nf,-nf:nf),msq46_1(-nf:nf,-nf:nf),
     & msq35_1_swap(-nf:nf,-nf:nf),msq36_1_swap(-nf:nf,-nf:nf),
     & msq45_1_swap(-nf:nf,-nf:nf),msq46_1_swap(-nf:nf,-nf:nf),
     & sub35_1,sub36_1,sub45_1,sub46_1,
     & msqbits35_1(12),msqbits36_1(12),msqbits45_1(12),msqbits46_1(12),
     & msqbits35_1_swap(12),msqbits36_1_swap(12),
     & msqbits45_1_swap(12),msqbits46_1_swap(12)
      external qqb_gmgmjt,qqb_gmgmjt_gvec,qqb_dirgam_g,qqb_dirgam_g_swap

      if (frag) then
        ndmax=10
      else
        ndmax=6
      endif

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & qqb_gmgmjt,qqb_gmgmjt_gvec)
      call dips(2,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & qqb_gmgmjt,qqb_gmgmjt_gvec)
      call dips(3,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     & qqb_gmgmjt,qqb_gmgmjt_gvec)
      call dips(4,p,2,6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
     & qqb_gmgmjt,qqb_gmgmjt_gvec)

c--- now the basic initial final ones
      call dips(5,p,1,5,6,sub15_6,sub15_6v,msq15_6,msq15_6v,
     & qqb_gmgmjt,qqb_gmgmjt_gvec)
c--- called for final initial the routine only supplies new values for
c--- sub... and sub...v and msqv
      call dips(5,p,5,6,1,sub56_1,sub56_1v,dummy,msq56_1v,
     & qqb_gmgmjt,qqb_gmgmjt_gvec)
      call dips(5,p,1,6,5,sub16_5,sub16_5v,msq16_5,msq16_5v,
     & qqb_gmgmjt,qqb_gmgmjt_gvec)

      call dips(6,p,2,6,5,sub26_5,sub26_5v,msq26_5,msq26_5v,
     & qqb_gmgmjt,qqb_gmgmjt_gvec)
      call dips(6,p,5,6,2,sub56_2,sub56_2v,dummy,msq56_2v,
     & qqb_gmgmjt,qqb_gmgmjt_gvec)
      call dips(6,p,2,5,6,sub25_6,sub25_6v,msq25_6,msq25_6v,
     & qqb_gmgmjt,qqb_gmgmjt_gvec)
c--- extra photon dipoles in the case of fragmentation
c--- note: qqb_dirgam_g_swap is just a wrapper routine to
c---       qqb_dirgam_g, with p4 and p5 interchanged;
c---       it is used, for instance, to obtain contributions
c---       from q q~ -> q~ q instead of the usual ordering q q~ -> q q~
      if (frag) then
         msqbits(:)=zip
         call dipsfrag(7,p,3,5,1,sub35_1,msq35_1,qqb_dirgam_g)
         msqbits35_1(:)=msqbits(:)
         call fill_gmgmjt_swap(7,msq35_1_swap)
         msqbits35_1_swap(:)=msqbits(:)

         msqbits(:)=zip
         call dipsfrag(8,p,3,6,1,sub36_1,msq36_1,qqb_dirgam_g)
         msqbits36_1(:)=msqbits(:)
         call fill_gmgmjt_swap(8,msq36_1_swap)
         msqbits36_1_swap(:)=msqbits(:)

         msqbits(:)=zip
         call dipsfrag(9,p,4,5,1,sub45_1,msq45_1,qqb_dirgam_g)
         msqbits45_1(:)=msqbits(:)
         call fill_gmgmjt_swap(9,msq45_1_swap)
         msqbits45_1_swap(:)=msqbits(:)

         msqbits(:)=zip
         call dipsfrag(10,p,4,6,1,sub46_1,msq46_1,qqb_dirgam_g)
         msqbits46_1(:)=msqbits(:)
         call fill_gmgmjt_swap(10,msq46_1_swap)
         msqbits46_1_swap(:)=msqbits(:)

         phot_dip(7:10)=.true.
      endif

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msq(nd,j,k)=0._dp
      enddo
      enddo
      enddo

c--- this is the factor to apply the right couplings
c--- in the gluon-gluon contribution
      facgg=
     &  (2._dp*Q(2)**4+real(nf-2,dp)*Q(1)**4)
     & /(2._dp*Q(2)**2+real(nf-2,dp)*Q(1)**2)

c--- 2-quark, 2-gluon subtraction pieces
      do j=-nf,nf
      do k=-nf,nf

      if ((j .ne. 0) .and. (k .ne. 0) .and. (j.ne.-k)) goto 19

c--- do only q-qb and qb-q cases
      if (  ((j > 0).and.(k < 0))
     & .or. ((j < 0).and.(k > 0))) then
C-----half=statistical factor
      msq(1,j,k)=-half*msq15_2(j,k)*sub15_2(qq)/xn
      msq(2,j,k)=-half*msq25_1(j,k)*sub25_1(qq)/xn
      msq(3,j,k)=-half*msq16_2(j,k)*sub16_2(qq)/xn
      msq(4,j,k)=-half*msq26_1(j,k)*sub26_1(qq)/xn
      msq(5,j,k)=half*xn*(
     &  msq15_6(j,k)*(sub15_6(qq)+0.5_dp*sub56_1(gg))
     & +0.5_dp*msq56_1v(j,k)*sub56_1v
     & +msq16_5(j,k)*(sub16_5(qq)+0.5_dp*sub56_1(gg))
     & +0.5_dp*msq56_1v(j,k)*sub56_1v)
      msq(6,j,k)=half*xn*(
     &  msq26_5(j,k)*(sub26_5(qq)+0.5_dp*sub56_2(gg))
     & +0.5_dp*msq56_2v(j,k)*sub56_2v
     & +msq25_6(j,k)*(sub25_6(qq)+0.5_dp*sub56_2(gg))
     & +0.5_dp*msq56_2v(j,k)*sub56_2v)
      elseif ((k == 0).and.(j.ne.0)) then
c--- q-g and qb-g cases
      msq(2,j,k)=2._dp*tr*msq25_1(j,-j)*sub25_1(qg)
      msq(3,j,k)=xn*msq16_2(j,k)*sub16_2(qq)
      msq(4,j,k)=xn*(msq26_1(j,k)*sub26_1(gg)+msq26_1v(j,k)*sub26_1v)
      msq(5,j,k)=-(msq16_5(j,k)*sub16_5(qq)+msq16_5(j,k)*sub56_1(qq))/xn
      msq(6,j,k)=xn*(msq26_5(j,k)*sub26_5(gg)+msq26_5v(j,k)*sub26_5v
     &              +msq26_5(j,k)*sub56_2(qq))
      if (frag) then
        msq(7,j,k)=Q(j)**2*msq35_1(j,k)*sub35_1/2._dp
        msq(9,j,k)=Q(j)**2*msq45_1(j,k)*sub45_1/2._dp
      endif

      elseif ((j == 0).and.(k.ne.0)) then
c--- g-q and g-qb cases
      msq(1,j,k)=2._dp*tr*msq15_2(-k,k)*sub15_2(qg)
      msq(3,j,k)=xn*(msq16_2(j,k)*sub16_2(gg)+msq16_2v(j,k)*sub16_2v)
      msq(4,j,k)=xn*msq26_1(j,k)*sub26_1(qq)
      msq(5,j,k)=xn*(msq16_5(j,k)*sub16_5(gg)+msq16_5v(j,k)*sub16_5v
     &              +msq15_6(j,k)*sub56_1(qq))
      msq(6,j,k)=-(msq26_5(j,k)*sub26_5(qq)+msq26_5(j,k)*sub56_2(qq))/xn
      if (frag) then
        msq(7,j,k)=Q(k)**2*msq35_1(j,k)*sub35_1/2._dp
        msq(9,j,k)=Q(k)**2*msq45_1(j,k)*sub45_1/2._dp
      endif

      elseif ((j == 0).and.(k == 0)) then
c--- g-g case (real process is g(p1)+g(p2) --> qb(p5)+q(p6)
c---Hence 15 split multiplies q(15)+g(p2)-->Z+q(p6)
c---Hence 25 split multiplies g(p1)+q(p25)-->Z+q(p6)
      msq(1,j,k)=(msq15_2(+1,k)+msq15_2(+2,k)+msq15_2(+3,k)
     &           +msq15_2(+4,k)+msq15_2(+5,k))*sub15_2(qg)*2._dp*tr
      msq(2,j,k)=(msq25_1(k,+1)+msq25_1(k,+2)+msq25_1(k,+3)
     &           +msq25_1(k,+4)+msq25_1(k,+5))*sub25_1(qg)*2._dp*tr
      msq(3,j,k)=(msq16_2(-5,k)+msq16_2(-4,k)+msq16_2(-3,k)
     &           +msq16_2(-2,k)+msq16_2(-1,k))*sub16_2(qg)*2._dp*tr
      msq(4,j,k)=(msq26_1(k,-5)+msq26_1(k,-4)+msq26_1(k,-3)
     &           +msq26_1(k,-2)+msq26_1(k,-1))*sub26_1(qg)*2._dp*tr
      if (frag) then
        msq(7,j,k)=msq35_1(j,k)*sub35_1/2._dp*facgg
        msq(8,j,k)=msq36_1(j,k)*sub36_1/2._dp*facgg
        msq(9,j,k)=msq45_1(j,k)*sub45_1/2._dp*facgg
        msq(10,j,k)=msq46_1(j,k)*sub46_1/2._dp*facgg
      endif

      endif

 19   continue
      enddo
      enddo


c--- 4-quark subtraction pieces
      do j=-nf,nf
      do k=-nf,nf

      if (((j > 0).and.(k > 0)) .or.
     &    ((j < 0).and.(k < 0))) then
c--q-q or qb-qb
      if (j==k) then
      msq(1,j,k)=msq(1,j,k)+0.5_dp*(xn-1._dp/xn)
     &  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
      msq(2,j,k)=msq(2,j,k)+0.5_dp*(xn-1._dp/xn)
     &  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
      msq(3,j,k)=msq(3,j,k)+0.5_dp*(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      msq(4,j,k)=msq(4,j,k)+0.5_dp*(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      else
      msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
     &  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
      msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      endif
      if (frag) then
        msq(7,j,k)=Q(j)**2*msq35_1(j,k)*sub35_1/2._dp
        msq(8,j,k)=Q(k)**2*msq36_1_swap(j,k)*sub36_1/2._dp
        msq(9,j,k)=Q(j)**2*msq45_1(j,k)*sub45_1/2._dp
        msq(10,j,k)=Q(k)**2*msq46_1_swap(j,k)*sub46_1/2._dp
      endif
      elseif (((j > 0).and.(k < 0))
     &    .or.((j < 0).and.(k > 0)))  then
c q-qbar or qbar-q
      if (j==-k) then
      msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
     &  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
      msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      msq(5,j,k)=msq(5,j,k)+tr*real(nflav,dp)
     &  *(msq16_5(j,k)*sub56_1(gq)-msq56_1v(j,k)*sub56_1v)
      msq(6,j,k)=msq(6,j,k)+tr*real(nflav,dp)
     &  *(msq26_5(j,k)*sub56_2(gq)-msq56_2v(j,k)*sub56_2v)
      if (frag) then
c--- note: subtraction terms use symmetry of qqb and qbq amplitudes
        if     ((abs(j) == 2) .or. (abs(j) == 4)) then
        msq(7,j,k)=sub35_1/2._dp*(
     &   Q(2)**2*(msqbits35_1(uub_uub)+msqbits35_1(uub_ccb))
     &  +Q(1)**2*(3._dp*msqbits35_1(uub_ddb)))
        msq(8,j,k)=sub36_1/2._dp*(
     &   Q(2)**2*(msqbits36_1_swap(uub_uub)+msqbits36_1_swap(uub_ccb))
     &  +Q(1)**2*(3._dp*msqbits36_1_swap(uub_ddb)))
        msq(9,j,k)=sub45_1/2._dp*(
     &   Q(2)**2*(msqbits45_1(uub_uub)+msqbits45_1(uub_ccb))
     &  +Q(1)**2*(3._dp*msqbits45_1(uub_ddb)))
        msq(10,j,k)=sub46_1/2._dp*(
     &   Q(2)**2*(msqbits46_1_swap(uub_uub)+msqbits46_1_swap(uub_ccb))
     &  +Q(1)**2*(3._dp*msqbits46_1_swap(uub_ddb)))
        elseif ((abs(j) == 1) .or. (abs(j) == 3)
     &     .or. (abs(j) == 5)) then
        msq(7,j,k)=sub35_1/2._dp*(
     &   Q(1)**2*(msqbits35_1(ddb_ddb)+2._dp*msqbits35_1(ddb_ssb))
     &  +Q(2)**2*(2._dp*msqbits35_1(ddb_uub)))
        msq(8,j,k)=sub36_1/2._dp*(
     &   Q(1)**2*(msqbits36_1_swap(ddb_ddb)
     &           +2._dp*msqbits36_1_swap(ddb_ssb))
     &  +Q(2)**2*(2._dp*msqbits36_1_swap(ddb_uub)))
        msq(9,j,k)=sub45_1/2._dp*(
     &   Q(1)**2*(msqbits45_1(ddb_ddb)+2._dp*msqbits45_1(ddb_ssb))
     &  +Q(2)**2*(2._dp*msqbits45_1(ddb_uub)))
        msq(10,j,k)=sub46_1/2._dp*(
     &   Q(1)**2*(msqbits46_1_swap(ddb_ddb)
     &           +2._dp*msqbits46_1_swap(ddb_ssb))
     &  +Q(2)**2*(2._dp*msqbits46_1_swap(ddb_uub)))
        endif
      endif
      else
      msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
     &  *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
      msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      if (frag) then
        if (j > 0) then
        msq(7,j,k)=Q(j)**2*msq35_1(j,k)*sub35_1/2._dp
        msq(8,j,k)=Q(k)**2*msq36_1_swap(j,k)*sub36_1/2._dp
        msq(9,j,k)=Q(j)**2*msq45_1(j,k)*sub45_1/2._dp
        msq(10,j,k)=Q(k)**2*msq46_1_swap(j,k)*sub46_1/2._dp
        else
        msq(7,j,k)=Q(j)**2*msq35_1_swap(j,k)*sub35_1/2._dp
        msq(8,j,k)=Q(k)**2*msq36_1(j,k)*sub36_1/2._dp
        msq(9,j,k)=Q(j)**2*msq45_1_swap(j,k)*sub45_1/2._dp
        msq(10,j,k)=Q(k)**2*msq46_1(j,k)*sub46_1/2._dp
        endif
      endif
      endif
c--qbar-q
c      elseif ((j < 0).and.(k > 0)) then
c      if (j==-k) then
c      msq(2,j,k)=msq(2,j,k)+(xn-1._dp/xn)
c     &  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
c      msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
c     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
c      msq(6,j,k)=msq(6,j,k)+2._dp*tr*real(nflav,dp)
c     & *(msq26_5(j,k)*sub56_2(gq)-msq56_2v(j,k)*sub56_2v)
c      else
c      msq(2,j,k)=msq(2,j,k)+(xn-1._dp/xn)
c     &  *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
c      msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
c     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
c      endif
      endif


      enddo
      enddo

      return
      end




      subroutine fill_gmgmjt_swap(nd,msq_swap)
      implicit none
      include 'types.f'
c--- routine that calls qqb_dirgam_g_swap and fills matrix elements
c---  input: nd (dipole number)
c--- output:msq_swap (array of matrix elements after 4<->5 swap)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'z_dip.f'
      include 'incldip.f'
      integer:: nd
      real(dp):: pdip(mxpart,4),msq_swap(-nf:nf,-nf:nf)

      if (incldip(nd)) then
        call getptilde(nd,pdip)
        pdip(4,:)=pdip(4,:)/z_dip(nd)
        call qqb_dirgam_g_swap(pdip,msq_swap)
      else
        msq_swap(:,:)=0._dp
      endif

      return
      end

