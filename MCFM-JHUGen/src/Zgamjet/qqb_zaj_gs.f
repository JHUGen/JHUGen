      subroutine qqb_zaj_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Authors: J. Campbell and H. Hartanto                             *
*     March, 2012.                                                     *
*                                                                      *
*    Matrix element SUBTRACTION squared averag'd over init'l colors    *
*    and spins                                                         *
*                                                                      *
*     q(-p1)+qbar(-p2) -->  Z + gamma(p5) + parton(p6)                 *
*                           |                                          *
*                            -->l(p3)+a(p4)                            *
************************************************************************

       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'nflav.f'
      include 'qqgg.f'
      include 'frag.f'
      include 'ewcharge.f'
      include 'phot_dip.f'
      include 'ipsgen.f'        
      integer:: j,k,nd,f
c --- remember: nd will count the dipoles      
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf),nup,ndo
      real(dp):: 
     & msq16_2(-nf:nf,-nf:nf),msq26_1(-nf:nf,-nf:nf),
     & msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     & msq16_7(-nf:nf,-nf:nf),msq27_6(-nf:nf,-nf:nf),
     & msq17_6(-nf:nf,-nf:nf),msq26_7(-nf:nf,-nf:nf),
     & msq67_1v(-nf:nf,-nf:nf),msq67_2v(-nf:nf,-nf:nf),
     & msq27_6v(-nf:nf,-nf:nf),msq27_1v(-nf:nf,-nf:nf),
     & msq16_7v(-nf:nf,-nf:nf),msq17_2v(-nf:nf,-nf:nf),
     & msq17_6v(-nf:nf,-nf:nf),msq26_7v(-nf:nf,-nf:nf),
     & msq26_1v(-nf:nf,-nf:nf),
     & msq16_2v(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),
     & sub16_2(4),sub26_1(4),sub17_2(4),sub27_1(4),
     & sub16_7(4),sub17_6(4),sub26_7(4),sub27_6(4),
     & sub67_1(4),sub67_2(4),sub67_1v,sub67_2v,
     & sub27_6v,sub26_1v,sub27_1v,sub17_6v,sub17_2v,sub16_2v,sub16_7v,
     & sub26_7v
      real(dp):: sub56_1,sub56_2,sub56_7,sub57_6,sub57_1,sub57_2
      real(dp):: msq56_1(-nf:nf,-nf:nf),msq56_2(-nf:nf,-nf:nf),
     &     msq56_7(-nf:nf,-nf:nf),msq57_1(-nf:nf,-nf:nf), 
     &     msq57_2(-nf:nf,-nf:nf),msq57_6(-nf:nf,-nf:nf)  
      real(dp)::
     &  msq56_1x(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     &  msq57_1x(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     &  msq56_7x(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     &  msq57_6x(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp):: 
     &  msq57_1x_swap(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     &  msq57_1_swap(-nf:nf,-nf:nf)
      real(dp):: 
     &  msq56_1x_swap(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf),
     &  msq56_1_swap(-nf:nf,-nf:nf)
      external qqb_zaj,qqb_zaj_gvec
      external qqb_z2jetx

      if((frag) .and. (ipsgen == 1)) then 
         ndmax=8
      else
         ndmax=6 
      endif

!---- Intialise Photon dipoles 
      do j=1,mxpart
        phot_dip(j)=.false.
      enddo


c--- calculate all the initial-initial dipoles
      call dips(1,p,1,6,2,sub16_2,sub16_2v,msq16_2,msq16_2v,
     & qqb_zaj,qqb_zaj_gvec)
      call dips(2,p,2,6,1,sub26_1,sub26_1v,msq26_1,msq26_1v,
     & qqb_zaj,qqb_zaj_gvec)
      call dips(3,p,1,7,2,sub17_2,sub17_2v,msq17_2,msq17_2v,
     & qqb_zaj,qqb_zaj_gvec)
      call dips(4,p,2,7,1,sub27_1,sub27_1v,msq27_1,msq27_1v,
     & qqb_zaj,qqb_zaj_gvec)

c--- now the basic initial final ones
      call dips(5,p,1,6,7,sub16_7,sub16_7v,msq16_7,msq16_7v,
     & qqb_zaj,qqb_zaj_gvec)
c--- called for final initial the routine only supplies new values for
c--- sub... and sub...v and msqv
      call dips(5,p,6,7,1,sub67_1,sub67_1v,dummy,msq67_1v,
     & qqb_zaj,qqb_zaj_gvec)
      call dips(5,p,1,7,6,sub17_6,sub17_6v,msq17_6,msq17_6v,
     & qqb_zaj,qqb_zaj_gvec)

      call dips(6,p,2,7,6,sub27_6,sub27_6v,msq27_6,msq27_6v,
     & qqb_zaj,qqb_zaj_gvec)
      call dips(6,p,6,7,2,sub67_2,sub67_2v,dummy,msq67_2v,
     & qqb_zaj,qqb_zaj_gvec)
      call dips(6,p,2,6,7,sub26_7,sub26_7v,msq26_7,msq26_7v,
     & qqb_zaj,qqb_zaj_gvec)

!--------------------------------------------------------------
!     
!     PHOTON DIPOLES 
!
!--------------------------------------------------------------

      if((frag) .and. (ipsgen == 1)) then 
         call dipsfragx(7,p,5,6,1,sub56_1,msq56_1,msq56_1x,qqb_z2jetx) 
         call fill_z2jet_swap(7,msq56_1_swap,msq56_1x_swap) 
    
         call dipsfragx(8,p,5,7,1,sub57_1,msq57_1,msq57_1x,qqb_z2jetx) 
         call fill_z2jet_swap(8,msq57_1_swap,msq57_1x_swap) 
         do j=7,8
            phot_dip(j)=.true. 
         enddo
      endif

      do j=-nf,nf
      do k=-nf,nf      
      do nd=1,ndmax
        msq(nd,j,k)=0._dp
      enddo
      enddo
      enddo
      
c--- subtraction pieces for real diagrams with 2 gluons and 2 quarks
      do j=-nf,nf
      do k=-nf,nf
      
      if ((j .ne. 0) .and. (k .ne. 0) .and. (j.ne.-k)) goto 19

c--- do only q-qb and qb-q cases      
      if (  ((j > 0).and.(k < 0))
     & .or. ((j < 0).and.(k > 0))) then
C-----half=statistical factor
      msq(1,j,k)=-half*msq16_2(j,k)*sub16_2(qq)/xn
      msq(2,j,k)=-half*msq26_1(j,k)*sub26_1(qq)/xn
      msq(3,j,k)=-half*msq17_2(j,k)*sub17_2(qq)/xn
      msq(4,j,k)=-half*msq27_1(j,k)*sub27_1(qq)/xn
      msq(5,j,k)=half*xn*(
     &  msq16_7(j,k)*(sub16_7(qq)+0.5_dp*sub67_1(gg))
     & +0.5_dp*msq67_1v(j,k)*sub67_1v
     & +msq17_6(j,k)*(sub17_6(qq)+0.5_dp*sub67_1(gg))
     & +0.5_dp*msq67_1v(j,k)*sub67_1v)
      msq(6,j,k)=half*xn*(
     &  msq27_6(j,k)*(sub27_6(qq)+0.5_dp*sub67_2(gg))
     & +0.5_dp*msq67_2v(j,k)*sub67_2v
     & +msq26_7(j,k)*(sub26_7(qq)+0.5_dp*sub67_2(gg))
     & +0.5_dp*msq67_2v(j,k)*sub67_2v)

      elseif ((k == 0).and.(j.ne.0)) then
c--- q-g and qb-g cases
      msq(2,j,k)=2._dp*tr*msq26_1(j,-j)*sub26_1(qg)
      msq(3,j,k)=xn*msq17_2(j,k)*sub17_2(qq)
      msq(4,j,k)=xn*(msq27_1(j,k)*sub27_1(gg)+msq27_1v(j,k)*sub27_1v)
      msq(5,j,k)=-(msq17_6(j,k)*sub17_6(qq)+msq17_6(j,k)*sub67_1(qq))/xn
      msq(6,j,k)=xn*(msq27_6(j,k)*sub27_6(gg)+msq27_6v(j,k)*sub27_6v
     &              +msq27_6(j,k)*sub67_2(qq))
 
      if((frag) .and. (ipsgen == 1)) then 
         msq(7,j,k)=Q(abs(j))**2*msq56_1(j,k)*sub56_1
      endif

      elseif ((j == 0).and.(k.ne.0)) then
c--- g-q and g-qb cases
      msq(1,j,k)=2._dp*tr*msq16_2(-k,k)*sub16_2(qg)
      msq(3,j,k)=xn*(msq17_2(j,k)*sub17_2(gg)+msq17_2v(j,k)*sub17_2v)
      msq(4,j,k)=xn*msq27_1(j,k)*sub27_1(qq)
      msq(5,j,k)=xn*(msq17_6(j,k)*sub17_6(gg)+msq17_6v(j,k)*sub17_6v
     &              +msq16_7(j,k)*sub67_1(qq))
      msq(6,j,k)=-(msq27_6(j,k)*sub27_6(qq)+msq27_6(j,k)*sub67_2(qq))/xn

      if((frag) .and. (ipsgen == 1)) then 
         msq(7,j,k)=Q(abs(k))**2*msq56_1(j,k)*sub56_1
      endif

      elseif ((j == 0).and.(k == 0)) then
c--- g-g case (real process is g(p1)+g(p2) --> qb(p6)+q(p7)
c---Hence 16 split multiplies q(16)+g(p2)-->Z+q(p7)
c---Hence 26 split multiplies g(p1)+q(p26)-->Z+q(p7)
      msq(1,j,k)=(msq16_2(-1,k)+msq16_2(-2,k)+msq16_2(-3,k)
     &           +msq16_2(-4,k)+msq16_2(-5,k))*sub16_2(qg)*2._dp*tr
      msq(2,j,k)=(msq26_1(k,-1)+msq26_1(k,-2)+msq26_1(k,-3)
     &           +msq26_1(k,-4)+msq26_1(k,-5))*sub26_1(qg)*2._dp*tr
      msq(3,j,k)=(msq17_2(+5,k)+msq17_2(+4,k)+msq17_2(+3,k)
     &           +msq17_2(+2,k)+msq17_2(+1,k))*sub17_2(qg)*2._dp*tr
      msq(4,j,k)=(msq27_1(k,+5)+msq27_1(k,+4)+msq27_1(k,+3)
     &           +msq27_1(k,+2)+msq27_1(k,+1))*sub27_1(qg)*2._dp*tr

      if((frag) .and. (ipsgen == 1)) then 
         msq(7,j,k)=2._dp*Q(2)**2*sub56_1
     & *(msq56_1x(0,j,k,2,-2)+msq56_1x(1,j,k,2,-2)+msq56_1x(2,j,k,2,-2))
     &             +3._dp*Q(1)**2*sub56_1
     & *(msq56_1x(0,j,k,1,-1)+msq56_1x(1,j,k,1,-1)+msq56_1x(2,j,k,1,-1))
         msq(8,j,k)=2._dp*Q(2)**2*sub57_1
     & *(msq57_1x(0,j,k,-2,2)+msq57_1x(1,j,k,-2,2)+msq57_1x(2,j,k,-2,2))
     &             +3._dp*Q(1)**2*sub57_1
     & *(msq57_1x(0,j,k,-1,1)+msq57_1x(1,j,k,-1,1)+msq57_1x(2,j,k,-1,1))
      endif

      endif

 19   continue
      enddo
      enddo

c--- subtraction pieces for real diagrams with 4 quarks
      do j=-nf,nf
      do k=-nf,nf      

      if (((j > 0).and.(k > 0)) .or. 
     &    ((j < 0).and.(k < 0))) then
c--q-q or qb-qb
      if (j==k) then
      msq(1,j,k)=msq(1,j,k)+0.5_dp*(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      msq(2,j,k)=msq(2,j,k)+0.5_dp*(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      msq(3,j,k)=msq(3,j,k)+0.5_dp*(xn-1._dp/xn)
     &  *(msq17_2(0,k)*sub17_2(gq)+msq17_2v(0,k)*sub17_2v)
      msq(4,j,k)=msq(4,j,k)+0.5_dp*(xn-1._dp/xn)
     &  *(msq27_1(j,0)*sub27_1(gq)+msq27_1v(j,0)*sub27_1v)
      else
      msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &  *(msq27_1(j,0)*sub27_1(gq)+msq27_1v(j,0)*sub27_1v)
      endif
      if((frag) .and. (ipsgen == 1)) then 
         msq(7,j,k)=Q(abs(j))**2*sub56_1
     & *(msq56_1x(0,j,k,j,k)+msq56_1x(1,j,k,j,k)+msq56_1x(2,j,k,j,k))
         msq(8,j,k)=Q(abs(k))**2*sub57_1
     & *(msq57_1x(0,j,k,k,j)+msq57_1x(1,j,k,k,j)+msq57_1x(2,j,k,k,j))
      endif

      elseif ((j > 0).and.(k < 0)) then
c q-qbar
      if (j==-k) then
      msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &  *(msq27_1(j,0)*sub27_1(gq)+msq27_1v(j,0)*sub27_1v)
      msq(5,j,k)=msq(5,j,k)+tr*real(nflav,dp)
     & *(msq17_6(j,k)*sub67_1(gq)-msq67_1v(j,k)*sub67_1v)
      msq(6,j,k)=msq(6,j,k)+tr*real(nflav,dp)
     & *(msq27_6(j,k)*sub67_2(gq)-msq67_2v(j,k)*sub67_2v)
      else 
      msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
     &  *(msq16_2(0,k)*sub16_2(gq)+msq16_2v(0,k)*sub16_2v)
      msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &  *(msq27_1(j,0)*sub27_1(gq)+msq27_1v(j,0)*sub27_1v)
      endif

      if((frag) .and. (ipsgen == 1)) then 
         if (j == -k) then
           do f=1,5
           msq(7,j,k)=msq(7,j,k)+Q(f)**2*sub56_1*(msq56_1x(0,j,k,f,-f)
     &                      +msq56_1x(1,j,k,f,-f)+msq56_1x(2,j,k,f,-f))
           msq(8,j,k)=msq(8,j,k)+Q(f)**2*sub57_1*(msq57_1x(0,j,k,-f,f)
     &                      +msq57_1x(1,j,k,-f,f)+msq57_1x(2,j,k,-f,f))
           enddo
         else
           msq(7,j,k)=Q(abs(j))**2*sub56_1
     &   *(msq56_1x(0,j,k,j,k)+msq56_1x(1,j,k,j,k)+msq56_1x(2,j,k,j,k))
           msq(8,j,k)=Q(abs(k))**2*sub57_1
     &   *(msq57_1x(0,j,k,k,j)+msq57_1x(1,j,k,k,j)+msq57_1x(2,j,k,k,j))
         endif
      endif

c--qbar-q
      elseif ((j < 0).and.(k > 0)) then
      if (j==-k) then
      msq(2,j,k)=msq(2,j,k)+(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &  *(msq17_2(0,k)*sub17_2(gq)+msq17_2v(0,k)*sub17_2v)
      msq(5,j,k)=msq(5,j,k)+tr*real(nflav,dp)
     & *(msq17_6(j,k)*sub67_1(gq)-msq67_1v(j,k)*sub67_1v)
      msq(6,j,k)=msq(6,j,k)+tr*real(nflav,dp)
     & *(msq27_6(j,k)*sub67_2(gq)-msq67_2v(j,k)*sub67_2v)
      else 
      msq(2,j,k)=msq(2,j,k)+(xn-1._dp/xn)
     &  *(msq26_1(j,0)*sub26_1(gq)+msq26_1v(j,0)*sub26_1v)
      msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &  *(msq17_2(0,k)*sub17_2(gq)+msq17_2v(0,k)*sub17_2v)

      endif

      if((frag) .and. (ipsgen == 1)) then 
         if (j == -k) then
           do f=1,5
           msq(7,j,k)=msq(7,j,k)+Q(f)**2*sub56_1*(msq56_1x(0,j,k,f,-f)
     &                      +msq56_1x(1,j,k,f,-f)+msq56_1x(2,j,k,f,-f))
           msq(8,j,k)=msq(8,j,k)+Q(f)**2*sub57_1*(msq57_1x(0,j,k,-f,f)
     &                      +msq57_1x(1,j,k,-f,f)+msq57_1x(2,j,k,-f,f))
           enddo
         else
           msq(7,j,k)=Q(abs(k))**2*sub56_1
     &   *(msq56_1x(0,j,k,k,j)+msq56_1x(1,j,k,k,j)+msq56_1x(2,j,k,k,j))
           msq(8,j,k)=Q(abs(j))**2*sub57_1
     &   *(msq57_1x(0,j,k,j,k)+msq57_1x(1,j,k,j,k)+msq57_1x(2,j,k,j,k))
         endif
      endif

      endif


      enddo
      enddo

      return
      end
      
      subroutine fill_z2jet_swap(nd,msq_swap,msq_x_swap)
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
      real(dp):: msq_x_swap(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp):: mqq(0:2,fn:nf,fn:nf),msq0,msq1,msq2
      real(dp):: msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp):: msqx_cs(0:2,-nf:nf,-nf:nf)
      real(dp):: msqd1(0:2,-nf:nf,-nf:nf),msqd2(0:2,-nf:nf,-nf:nf) 
      
      if (incldip(nd)) then
        call getptilde(nd,pdip)
        pdip(5,:)=pdip(5,:)/z_dip(nd)
        call qqb_z2jetx_swap(pdip,msq_swap,msqd1,msq_x_swap,msqd2)
      else
        msq_swap(:,:)=0._dp
      endif
       
      return
      end
