      subroutine qqb_dirgam_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell                                            *
*     October, 2002.                                                   *
*     Modified 2011   By C. Williams                                   *
*    Matrix element SUBTRACTION squared averag'd over init'l colors    *
*    and spins                                                         *
*     f(-p1) + f(-p2) -->  gamma(p3) + parton(p4) + parton(p5)         *
************************************************************************

       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'ewcharge.f'
      include 'frag.f'
      include 'phot_dip.f'
      include 'msqbits.f'
      real(dp):: msqbits34_1(12),msqbits35_1(12)
      real(dp):: msqbits34_1_swap(12),msqbits35_1_swap(12)
      real(dp):: msq34_1_swap(-nf:nf,-nf:nf)
     &,msq35_1_swap(-nf:nf,-nf:nf)
      integer:: j,k,nd
c --- remember: nd will count the dipoles
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp):: 
     & msq14_2(-nf:nf,-nf:nf),msq24_1(-nf:nf,-nf:nf),
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     & msq14_5(-nf:nf,-nf:nf),msq25_4(-nf:nf,-nf:nf),
     & msq45_1v(-nf:nf,-nf:nf),msq45_2v(-nf:nf,-nf:nf),
     & msq25_4v(-nf:nf,-nf:nf),msq25_1v(-nf:nf,-nf:nf),
     & msq14_5v(-nf:nf,-nf:nf),msq15_2v(-nf:nf,-nf:nf),
     & msq45_2(-nf:nf,-nf:nf),msq45_1(-nf:nf,-nf:nf),
     & msq15_4(-nf:nf,-nf:nf),msq15_4v(-nf:nf,-nf:nf),
     & msq24_5(-nf:nf,-nf:nf),msq24_5v(-nf:nf,-nf:nf),
     & msq24_1v(-nf:nf,-nf:nf),
     & msq14_2v(-nf:nf,-nf:nf),
     & msq34_1(-nf:nf,-nf:nf),
     & msq34_2(-nf:nf,-nf:nf),
     & msq35_1(-nf:nf,-nf:nf),
     & msq35_2(-nf:nf,-nf:nf), 
     & msq34_5(-nf:nf,-nf:nf),msq35_4(-nf:nf,-nf:nf),
     & sub14_2(4),sub24_1(4),sub15_2(4),sub25_1(4),
     & sub14_5(4),sub15_4(4),sub25_4(4),sub24_5(4),
     & sub45_1(4),sub45_2(4),sub45_1v,sub45_2v,
     & sub25_4v,sub24_1v,sub25_1v,sub15_4v,sub15_2v,sub14_2v,sub14_5v,
     & sub24_5v,sub35_4,sub34_5,sub35_1,sub34_1,sub35_2,sub34_2 
      external qqb_dirgam,qqb_dirgam_gvec,qqb_2j_t,qqb_2j_s,
     & qqb_2jnoggswap
      external qqb_2jet
c      external qqb_2jnogg,qqb_2j_sqqb_2j_sswap,donothing_gvec
      if (frag) then 
         ndmax=8
      else
         ndmax=6
      endif
      

c--- calculate all the initial-initial dipoles
      call dips(1,p,1,4,2,sub14_2,sub14_2v,msq14_2,msq14_2v,
     & qqb_dirgam,qqb_dirgam_gvec)
      call dips(2,p,2,4,1,sub24_1,sub24_1v,msq24_1,msq24_1v,
     & qqb_dirgam,qqb_dirgam_gvec)
      call dips(3,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & qqb_dirgam,qqb_dirgam_gvec)
      call dips(4,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & qqb_dirgam,qqb_dirgam_gvec)

c--- now the basic initial final ones
      call dips(5,p,1,4,5,sub14_5,sub14_5v,msq14_5,msq14_5v,
     & qqb_dirgam,qqb_dirgam_gvec)
c--- called for final initial the routine only supplies 
c----new values for
c--- sub... and sub...v and msqv
      call dips(5,p,4,5,1,sub45_1,sub45_1v,msq45_1,msq45_1v,
     & qqb_dirgam,qqb_dirgam_gvec)
      call dips(5,p,1,5,4,sub15_4,sub15_4v,msq15_4,msq15_4v,
     & qqb_dirgam,qqb_dirgam_gvec)

      call dips(6,p,2,4,5,sub24_5,sub24_5v,msq24_5,msq24_5v,
     & qqb_dirgam,qqb_dirgam_gvec)
      call dips(6,p,4,5,2,sub45_2,sub45_2v,msq45_2,msq45_2v,
     & qqb_dirgam,qqb_dirgam_gvec)
      call dips(6,p,2,5,4,sub25_4,sub25_4v,msq25_4,msq25_4v,
     & qqb_dirgam,qqb_dirgam_gvec)


     
      if (frag) then
         msqbits(:)=zip
         call dipsfrag(7,p,3,4,1,sub34_1,msq34_1,qqb_2jet) 
         msqbits34_1(:)=msqbits(:) 
         call fill_dirgam_swap(7,msq34_1_swap)
         msqbits34_1_swap(:)=msqbits(:) 
         
         msqbits(:)=zip
         call dipsfrag(8,p,3,5,1,sub35_1,msq35_1,qqb_2jet)
         msqbits35_1(:)=msqbits(:)
         call fill_dirgam_swap(8,msq35_1_swap) 
         msqbits35_1_swap(:)=msqbits(:) 

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


      
      do j=-nf,nf
      do k=-nf,nf
      
      if ((j .ne. 0) .and. (k .ne. 0) .and. (j.ne.-k)) goto 19

c--- do only q-qb and qb-q cases      
      if (  ((j > 0).and.(k < 0))
     & .or. ((j < 0).and.(k > 0))) then
C-----half=statistical factor
      msq(1,j,k)=-half*msq14_2(j,k)*sub14_2(qq)/xn
      msq(2,j,k)=-half*msq24_1(j,k)*sub24_1(qq)/xn
      msq(3,j,k)=-half*msq15_2(j,k)*sub15_2(qq)/xn
      msq(4,j,k)=-half*msq25_1(j,k)*sub25_1(qq)/xn
      msq(5,j,k)=half*xn*(
     &  msq14_5(j,k)*(sub14_5(qq)+0.5_dp*sub45_1(gg))
     & +0.5_dp*msq45_1v(j,k)*sub45_1v
     & +msq14_5(j,k)*(sub15_4(qq)+0.5_dp*sub45_1(gg))
     & +0.5_dp*msq45_1v(j,k)*sub45_1v)
      msq(6,j,k)=half*xn*(
     &  msq25_4(j,k)*(sub25_4(qq)+0.5_dp*sub45_2(gg))
     & +0.5_dp*msq45_2v(j,k)*sub45_2v
     & +msq25_4(j,k)*(sub24_5(qq)+0.5_dp*sub45_2(gg))
     & +0.5_dp*msq45_2v(j,k)*sub45_2v)

      elseif ((k == 0).and.(j.ne.0)) then
c--- q-g and qb-g cases
      msq(2,j,k)=2._dp*tr*msq24_1(j,-j)*sub24_1(qg)
      msq(3,j,k)=xn*msq15_2(j,k)*sub15_2(qq)
      msq(4,j,k)=xn*(msq25_1(j,k)*sub25_1(gg)+msq25_1v(j,k)*sub25_1v)
      msq(5,j,k)=-(msq15_4(j,k)*sub15_4(qq)+msq15_4(j,k)*sub45_1(qq))/xn
      msq(6,j,k)=xn*(msq25_4(j,k)*sub25_4(gg)+msq25_4v(j,k)*sub25_4v
     &              +msq25_4(j,k)*sub45_2(qq))
 
      if(frag) then 
      msq(7,j,k)=Q(j)**2*msq34_1(j,k)*sub34_1
      endif

      elseif ((j == 0).and.(k.ne.0)) then
c--- g-q and g-qb cases
      msq(1,j,k)=2._dp*tr*msq14_2(-k,k)*sub14_2(qg)
      msq(3,j,k)=xn*(msq15_2(j,k)*sub15_2(gg)+msq15_2v(j,k)*sub15_2v)
      msq(4,j,k)=xn*msq25_1(j,k)*sub25_1(qq)
      msq(5,j,k)=xn*(msq15_4(j,k)*sub15_4(gg)+msq15_4v(j,k)*sub15_4v
     &              +msq15_4(j,k)*sub45_1(qq))      
      msq(6,j,k)=-(msq25_4(j,k)*sub25_4(qq)+msq25_4(j,k)*sub45_2(qq))/xn
      
      if(frag) then 
         msq(7,j,k)=Q(k)**2*msq34_1(j,k)*sub34_1
      endif

      elseif ((j == 0).and.(k == 0)) then
c--- g-g case (real process is g(p1)+g(p2) --> gamma(p3)+q(p4)+qb(p5)
c---Hence 14 split multiplies qb(14)+g(p2) --> gamma(p3)+qb(p5)
c---Hence 24 split multiplies g(p1)+qb(p24) --> gamma(p3)+qb(p5)
      msq(1,j,k)=(msq14_2(-1,k)+msq14_2(-2,k)+msq14_2(-3,k)
     &           +msq14_2(-4,k)+msq14_2(-5,k))*sub14_2(qg)*2._dp*tr
      msq(2,j,k)=(msq24_1(k,-1)+msq24_1(k,-2)+msq24_1(k,-3)
     &           +msq24_1(k,-4)+msq24_1(k,-5))*sub24_1(qg)*2._dp*tr
      msq(3,j,k)=(msq15_2(+5,k)+msq15_2(+4,k)+msq15_2(+3,k)
     &           +msq15_2(+2,k)+msq15_2(+1,k))*sub15_2(qg)*2._dp*tr
      msq(4,j,k)=(msq25_1(k,+5)+msq25_1(k,+4)+msq25_1(k,+3)
     &           +msq25_1(k,+2)+msq25_1(k,+1))*sub25_1(qg)*2._dp*tr
      
      if (frag) then 
      msq(7,j,k)=(real(nf-2,dp)*Q(1)**2+2._dp*Q(2)**2)
     &  *msq34_1(j,k)*sub34_1
      msq(8,j,k)=(real(nf-2,dp)*Q(1)**2+2._dp*Q(2)**2)
     &  *msq35_1_swap(j,k)*sub35_1
      endif

      endif

 19   continue
      enddo
      enddo


      do j=-nf,nf
      do k=-nf,nf      

         if (((j > 0).and.(k > 0)) .or. 
     &        ((j < 0).and.(k < 0))) then
c---  q-q or qb-qb
            if (j==k) then
               msq(1,j,k)=msq(1,j,k)+0.5_dp*(xn-1._dp/xn)
     &              *(msq14_2(0,k)*sub14_2(gq)+msq14_2v(0,k)*sub14_2v)
               msq(2,j,k)=msq(2,j,k)+0.5_dp*(xn-1._dp/xn)
     &              *(msq24_1(j,0)*sub24_1(gq)+msq24_1v(j,0)*sub24_1v)
               msq(3,j,k)=msq(3,j,k)+0.5_dp*(xn-1._dp/xn)
     &              *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
               msq(4,j,k)=msq(4,j,k)+0.5_dp*(xn-1._dp/xn)
     &              *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
               
               if(frag) then 
                  msq(7,j,k)=Q(j)**2*msq34_1(j,k)*sub34_1
                  msq(8,j,k)=Q(k)**2*msq35_1_swap(j,k)*sub35_1
               endif

            else
               msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
     &              *(msq14_2(0,k)*sub14_2(gq)+msq14_2v(0,k)*sub14_2v)
               msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &              *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)

               if(frag) then 
                  msq(7,j,k)=Q(j)**2*msq34_1(j,k)*sub34_1
                  msq(8,j,k)=Q(k)**2*msq35_1_swap(j,k)*sub35_1
               endif

               

            endif
         elseif ((j > 0).and.(k < 0)) then

c--- q-qbar
            if (j==-k) then
               msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
     &              *(msq14_2(0,k)*sub14_2(gq)+msq14_2v(0,k)*sub14_2v)
               msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &              *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
               msq(6,j,k)=msq(6,j,k)+2._dp*tr*real(nf,dp)
     &              *(msq25_4(j,k)*sub45_2(gq)-msq45_2v(j,k)*sub45_2v)

                  
               if(frag) then    
                  if  ((abs(j) == 2) .or. (abs(j) == 4)) then
                     msq(7,j,k)=sub34_1*(
     &               Q(2)**2*(msqbits34_1(uub_uub)+msqbits34_1(uub_ccb))
     &               +Q(1)**2*(3._dp*msqbits34_1(uub_ddb)))
                     msq(8,j,k)=sub35_1*(
     &               Q(2)**2*(msqbits35_1_swap(uub_uub)
     &                    +1._dp*msqbits35_1_swap(uub_ccb))
     &               +Q(1)**2*(3._dp*msqbits35_1_swap(uub_ddb)))
                  else
                     msq(7,j,k)=sub34_1*(
     &               Q(1)**2*(msqbits34_1(ddb_ddb)
     &                    +2._dp*msqbits34_1(ddb_ssb))
     &               +Q(2)**2*(2._dp*msqbits34_1(ddb_uub)))
    
                     msq(8,j,k)=sub35_1*(
     &               Q(1)**2*(msqbits35_1_swap(ddb_ddb)
     &                    +2._dp*msqbits35_1_swap(ddb_ssb))
     &               +Q(2)**2*(2._dp*msqbits35_1_swap(ddb_uub)))
                     endif
               endif
               
            else 
               msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
     &              *(msq14_2(0,k)*sub14_2(gq)+msq14_2v(0,k)*sub14_2v)
               msq(4,j,k)=msq(4,j,k)+(xn-1._dp/xn)
     &              *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
               
               if(frag) then 
                  msq(7,j,k)=Q(j)**2*msq34_1(j,k)*sub34_1
                  msq(8,j,k)=Q(k)**2*msq35_1_swap(j,k)*sub35_1
               endif
            endif
            
 

c--- qbar-q
         elseif ((j < 0).and.(k > 0)) then
            if (j==-k) then               
               msq(2,j,k)=msq(2,j,k)+(xn-1._dp/xn)
     &              *(msq24_1(j,0)*sub24_1(gq)+msq24_1v(j,0)*sub24_1v)
               msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &              *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
               msq(6,j,k)=msq(6,j,k)+2._dp*tr*real(nf,dp)
     &              *(msq25_4(j,k)*sub45_2(gq)-msq45_2v(j,k)*sub45_2v)
               
               if(frag) then    
                  if  ((abs(j) == 2) .or. (abs(j) == 4)) then
                     msq(7,j,k)=sub34_1*(
     &               Q(2)**2*(msqbits34_1(ubu_uub)
     &                       +msqbits34_1(ubu_ccb))
     &               +Q(1)**2*(3._dp*msqbits34_1(ubu_ddb)))
                     msq(8,j,k)=sub35_1*(
     &               Q(2)**2*(msqbits35_1_swap(ubu_uub)
     &                       +msqbits35_1_swap(ubu_ccb))
     &               +Q(1)**2*(3._dp*msqbits35_1_swap(ubu_ddb)))
                  else
                     msq(7,j,k)=sub34_1*(
     &               Q(1)**2*(msqbits34_1(dbd_ddb)
     &                   +2._dp*msqbits34_1(dbd_ssb))
     &               +Q(2)**2*(2._dp*msqbits34_1(dbd_uub)))
    
                     msq(8,j,k)=sub35_1*(
     &               Q(1)**2*(msqbits35_1_swap(dbd_ddb)
     &                   +2._dp*msqbits35_1_swap(dbd_ssb))
     &               +Q(2)**2*(2._dp*msqbits35_1_swap(dbd_uub)))
                     endif
               
               
               endif
               
            else 
               msq(2,j,k)=msq(2,j,k)+(xn-1._dp/xn)
     &              *(msq24_1(j,0)*sub24_1(gq)+msq24_1v(j,0)*sub24_1v)
               msq(3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &              *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)               
               if(frag) then      
                  msq(7,j,k)=Q(k)**2*msq34_1(j,k)*sub34_1
                  msq(8,j,k)=Q(j)**2*msq35_1_swap(j,k)*sub35_1
               endif
                              
            endif
         endif


      enddo
      enddo

 
      return
      end
      

      subroutine fill_dirgam_swap(nd,msq_swap)
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
        pdip(3,:)=pdip(3,:)/z_dip(nd)
        call qqb_2jet_swap(pdip,msq_swap)
      else
        msq_swap(:,:)=0._dp
      endif
       
      return
      end
