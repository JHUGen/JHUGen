      subroutine qqb_hflgam_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell                                            *
*     January, 2013.                                                   *
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
      include 'heavyflav.f'
      include 'masses.f'
      include 'phot_dip.f'
      integer:: j,k,nd,n3,n4,n5,n6
c --- remember: nd will count the dipoles
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf),dot
      real(dp):: 
c     & msq14_2(-nf:nf,-nf:nf),msq24_1(-nf:nf,-nf:nf),
     & msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
c     & msq14_5(-nf:nf,-nf:nf),
     & msq25_4(-nf:nf,-nf:nf),
     & msq45_1v(-nf:nf,-nf:nf),msq45_2v(-nf:nf,-nf:nf),
     & msq25_4v(-nf:nf,-nf:nf),msq25_1v(-nf:nf,-nf:nf),
c     & msq14_5v(-nf:nf,-nf:nf),
     & msq15_2v(-nf:nf,-nf:nf),
     & msq45_2(-nf:nf,-nf:nf),msq45_1(-nf:nf,-nf:nf),
     & msq15_4(-nf:nf,-nf:nf),msq15_4v(-nf:nf,-nf:nf),
c     & msq24_5(-nf:nf,-nf:nf),msq24_5v(-nf:nf,-nf:nf),
c     & msq24_1v(-nf:nf,-nf:nf),
c     & msq14_2v(-nf:nf,-nf:nf),
     & msq34_1(-nf:nf,-nf:nf),
     & msq34_2(-nf:nf,-nf:nf),
     & msq35_1(-nf:nf,-nf:nf),
     & msq35_2(-nf:nf,-nf:nf), 
     & msq34_5(-nf:nf,-nf:nf),msq35_4(-nf:nf,-nf:nf),
c     & sub14_2(4),sub24_1(4),
     & sub15_2(4),sub25_1(4),
c     & sub14_5(4),
     & sub15_4(4),sub25_4(4),
c     & sub24_5(4),
     & sub45_1(4),sub45_2(4),sub45_1v,sub45_2v,
     & sub25_4v,sub25_1v,sub15_4v,sub15_2v,
C    & sub14_2v,sub14_5v,sub24_1v,
c     & sub24_5v,
     & sub35_4,sub34_5,sub35_1,sub34_1,sub35_2,sub34_2 
      external qqb_hflgam,qqb_hflgam_gvec,qqb_2j_t,qqb_2j_s,
     & qqb_2jnoggswap
c      external qqb_2jnogg,qqb_2j_sqqb_2j_sswap,donothing_gvec
      if (frag) then 
         ndmax=12
      else
         ndmax=4
      endif

c-- change of names      
      n3=1
      n4=2
      n5=3
      n6=4

c--- calculate all the initial-initial dipoles
c      call dips(1,p,1,4,2,sub14_2,sub14_2v,msq14_2,msq14_2v,
c     & qqb_hflgam,qqb_hflgam_gvec)
c      call dips(2,p,2,4,1,sub24_1,sub24_1v,msq24_1,msq24_1v,
c     & qqb_hflgam,qqb_hflgam_gvec)
      call dips(n3,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & qqb_hflgam,qqb_hflgam_gvec)
      call dips(n4,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & qqb_hflgam,qqb_hflgam_gvec)

c--- now the basic initial final ones
c      call dips(n5,p,1,4,5,sub14_5,sub14_5v,msq14_5,msq14_5v,
c     & qqb_hflgam,qqb_hflgam_gvec)
c--- called for final initial the routine only supplies 
c----new values for
c--- sub... and sub...v and msqv
      call dips(n5,p,4,5,1,sub45_1,sub45_1v,msq45_1,msq45_1v,
     & qqb_hflgam,qqb_hflgam_gvec)
      call dips(n5,p,1,5,4,sub15_4,sub15_4v,msq15_4,msq15_4v,
     & qqb_hflgam,qqb_hflgam_gvec)

c      call dips(n6,p,2,4,5,sub24_5,sub24_5v,msq24_5,msq24_5v,
c     & qqb_hflgam,qqb_hflgam_gvec)
      call dips(n6,p,4,5,2,sub45_2,sub45_2v,msq45_2,msq45_2v,
     & qqb_hflgam,qqb_hflgam_gvec)
      call dips(n6,p,2,5,4,sub25_4,sub25_4v,msq25_4,msq25_4v,
     & qqb_hflgam,qqb_hflgam_gvec)

C------debug
      frag=.false.
C------debug
                 
      if (frag) then
      call dipsfrag(7,p,3,4,1,sub34_1,msq34_1,qqb_2j_t) 
      call dipsfrag(8,p,3,4,2,sub34_2,msq34_2,qqb_2j_t)
      call dipsfrag(9,p,3,4,5,sub34_5,msq34_5,qqb_2j_s)
      call dipsfrag(10,p,3,5,4,sub35_4,msq35_4,qqb_2j_s)
      call dipsfrag(11,p,3,5,1,sub35_1,msq35_1,qqb_2jnoggswap) 
      call dipsfrag(12,p,3,5,2,sub35_2,msq35_2,qqb_2jnoggswap)
      do j=7,12
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


      
      do j=0,flav,flav
      do k=0,flav,flav
      
c      if ((j .ne. 0) .and. (k .ne. 0) .and. (j.ne.-k)) goto 19

c--- do only q-qb and qb-q cases      
c      if (  ((j > 0).and.(k < 0))
c     & .or. ((j < 0).and.(k > 0))) then
C-----half=statistical factor
c      msq(1,j,k)=-half*msq14_2(j,k)*sub14_2(qq)/xn
c      msq(2,j,k)=-half*msq24_1(j,k)*sub24_1(qq)/xn
c      msq(3,j,k)=-half*msq15_2(j,k)*sub15_2(qq)/xn
c      msq(4,j,k)=-half*msq25_1(j,k)*sub25_1(qq)/xn
c      msq(5,j,k)=half*xn*(
c     &  msq14_5(j,k)*(sub14_5(qq)+zip*0.5_dp*sub45_1(gg))
c     & +0.5_dp*msq45_1v(j,k)*sub45_1v
c     & +msq14_5(j,k)*(sub15_4(qq)+0.5_dp*sub45_1(gg))
c     & +0.5_dp*msq45_1v(j,k)*sub45_1v)
c      msq(6,j,k)=half*xn*(
c     &  msq25_4(j,k)*(sub25_4(qq)+zip*0.5_dp*sub45_2(gg))
c     & +0.5_dp*msq45_2v(j,k)*sub45_2v
c     & +msq25_4(j,k)*(sub24_5(qq)+0.5_dp*sub45_2(gg))
c     & +0.5_dp*msq45_2v(j,k)*sub45_2v)

      if ((k == 0).and.(j.ne.0)) then
c--- q-g and qb-g cases
c      msq(2,j,k)=2._dp*tr*msq24_1(j,-j)*sub24_1(qg)
      msq(n3,j,k)=xn*msq15_2(j,k)*sub15_2(qq)
      msq(n4,j,k)=xn*(msq25_1(j,k)*sub25_1(gg)+msq25_1v(j,k)*sub25_1v)
      msq(n5,j,k)=-(msq15_4(j,k)*sub15_4(qq)+msq15_4(j,k)*sub45_1(qq))
     & /xn
      msq(n6,j,k)=xn*(msq25_4(j,k)*sub25_4(gg)+msq25_4v(j,k)*sub25_4v
     &              +msq25_4(j,k)*sub45_2(qq))
 
      if(frag) then 
      msq(7,j,k)=Q(abs(j))**2*msq34_1(j,k)*sub34_1
      endif

      elseif ((j == 0).and.(k.ne.0)) then
c--- g-q and g-qb cases
c      msq(1,j,k)=2._dp*tr*msq14_2(-k,k)*sub14_2(qg)
      msq(n3,j,k)=xn*(msq15_2(j,k)*sub15_2(gg)+msq15_2v(j,k)*sub15_2v)
      msq(n4,j,k)=xn*msq25_1(j,k)*sub25_1(qq)
      msq(n5,j,k)=xn*(msq15_4(j,k)*sub15_4(gg)+msq15_4v(j,k)*sub15_4v
     &              +msq15_4(j,k)*sub45_1(qq))      
      msq(n6,j,k)=-(msq25_4(j,k)*sub25_4(qq)+msq25_4(j,k)*sub45_2(qq))
     & /xn
      
      if(frag) then 
         msq(8,j,k)=Q(abs(k))**2*msq34_2(j,k)*sub34_2
      endif

      elseif ((j == 0).and.(k == 0)) then
c--- g-g case (real process is g(p1)+g(p2) --> gamma(p3)+q(p4)+qb(p5)
c---Hence 14 split multiplies qb(14)+g(p2) --> gamma(p3)+qb(p5)
c---Hence 24 split multiplies g(p1)+qb(p24) --> gamma(p3)+qb(p5)
c      msq(1,j,k)=(msq14_2(-1,k)+msq14_2(-2,k)+msq14_2(-3,k)
c     &           +msq14_2(-4,k)+msq14_2(-5,k))*sub14_2(qg)*2._dp*tr
c      msq(2,j,k)=(msq24_1(k,-1)+msq24_1(k,-2)+msq24_1(k,-3)
c     &           +msq24_1(k,-4)+msq24_1(k,-5))*sub24_1(qg)*2._dp*tr
      msq(n3,j,k)=msq15_2(flav,k)*sub15_2(qg)*2._dp*tr
      msq(n4,j,k)=msq25_1(k,flav)*sub25_1(qg)*2._dp*tr
      
      if (frag) then 
      msq(9,j,k)=(real(nf-2,dp)*Q(1)**2+2._dp*Q(2)**2)
     &  *msq34_5(j,k)*sub34_5
      msq(10,j,k)=(real(nf-2,dp)*Q(1)**2+2._dp*Q(2)**2)
     &  *msq35_4(j,k)*sub35_4
      endif

      endif

c 19   continue
      enddo
      enddo


      do j=-nf,nf
      do k=-nf,nf      

         if (((j > 0).and.(k > 0))
c     & .or. ((j < 0).and.(k < 0))
     &       ) then
c---  q-q or qb-qb
            if ((j == flav) .and. (k. eq. flav)) then
c               msq(1,j,k)=msq(1,j,k)+0.5_dp*(xn-1._dp/xn)
c     &              *(msq14_2(0,k)*sub14_2(gq)+msq14_2v(0,k)*sub14_2v)
c               msq(2,j,k)=msq(2,j,k)+0.5_dp*(xn-1._dp/xn)
c     &              *(msq24_1(j,0)*sub24_1(gq)+msq24_1v(j,0)*sub24_1v)
               msq(n3,j,k)=msq(n3,j,k)+0.5_dp*(xn-1._dp/xn)
     &              *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
               msq(n4,j,k)=msq(n4,j,k)+0.5_dp*(xn-1._dp/xn)
     &              *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
               
               if(frag) then 
                  msq(7,j,k)=Q(abs(j))**2*msq34_1(j,k)*sub34_1
                  msq(12,j,k)=Q(abs(k))**2*msq35_2(j,k)*sub35_2
               endif

            else
c               msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
c     &              *(msq14_2(0,k)*sub14_2(gq)+msq14_2v(0,k)*sub14_2v)
               if (j == flav) then
               msq(n4,j,k)=msq(n4,j,k)+(xn-1._dp/xn)
     &              *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
               endif
               if (k == flav) then
               msq(n3,j,k)=msq(n3,j,k)+(xn-1._dp/xn)
     &              *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
               endif
               
               if(frag) then 
                  msq(7,j,k)=Q(abs(j))**2*msq34_1(j,k)*sub34_1
                  msq(12,j,k)=Q(abs(k))**2*msq35_2(j,k)*sub35_2
               endif

               

            endif
         elseif ((j > 0).and.(k < 0)) then

c--- q-qbar
            if (j==-k) then
c               msq(n1,j,k)=msq(n1,j,k)+(xn-1._dp/xn)
c     &              *(msq14_2(0,k)*sub14_2(gq)+msq14_2v(0,k)*sub14_2v)
               msq(n4,j,k)=msq(n4,j,k)+(xn-1._dp/xn)
     &              *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
c               msq(n6,j,k)=msq(n6,j,k)+2._dp*tr*real(nf,dp)
c     &              *(msq25_4(j,k)*sub45_2(gq)-msq45_2v(j,k)*sub45_2v)

                  
               if(frag) then    
!-----Initial-final dipoles (t channel) 
                  msq(7,j,k)=Q(abs(j))**2*msq34_1(j,k)*sub34_1
!                  msq(8,j,k)=0._dp*Q(abs(k))**2*msq34_2(j,k)*sub34_2
                  msq(11,j,k)=Q(abs(j))**2*msq35_1(j,k)*sub35_1
!                  msq(12,j,k)=0._dp*Q(abs(k))**2*msq35_2(j,k)*sub35_2
!-----Final-final dipoles (nf s-channel) 
                  msq(9,j,k)=(real(nf-2,dp)*Q(1)**2+2._dp*Q(2)**2)
     &                 *msq34_5(j,k)*sub34_5
                  msq(10,j,k)=(real(nf-2,dp)*Q(1)**2+2._dp*Q(2)**2)
     &                 *msq35_4(j,k)*sub35_4                    
               endif
               
            else 
c               msq(1,j,k)=msq(1,j,k)+(xn-1._dp/xn)
c     &              *(msq14_2(0,k)*sub14_2(gq)+msq14_2v(0,k)*sub14_2v)
               msq(n4,j,k)=msq(n4,j,k)+(xn-1._dp/xn)
     &              *(msq25_1(j,0)*sub25_1(gq)+msq25_1v(j,0)*sub25_1v)
               
               if(frag) then 
                  msq(7,j,k)=Q(abs(j))**2*msq34_1(j,k)*sub34_1
                  msq(12,j,k)=Q(abs(k))**2*msq35_2(j,k)*sub35_2
               endif
            endif
            
 

c--- qbar-q
         elseif ((j < 0).and.(k > 0)) then
            if (j==-k) then               
c               msq(n2,j,k)=msq(2,j,k)+(xn-1._dp/xn)
c     &              *(msq24_1(j,0)*sub24_1(gq)+msq24_1v(j,0)*sub24_1v)
               msq(n3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &              *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
c               msq(n6,j,k)=msq(6,j,k)+2._dp*tr*real(nf,dp)
c     &              *(msq25_4(j,k)*sub45_2(gq)-msq45_2v(j,k)*sub45_2v)
               
               if(frag) then    
!----- Initial-final dipoles (t channel) 
                  msq(7,j,k)=Q(abs(j))**2*msq34_1(j,k)*sub34_1
!                  msq(8,j,k)=Q(abs(k))**2*msq34_2(j,k)*sub34_2
                  msq(11,j,k)=Q(abs(j))**2*msq35_1(j,k)*sub35_1
!                  msq(12,j,k)=Q(abs(k))**2*msq35_2(j,k)*sub35_2
!-----Final-final dipoles (nf s-channel) 
                  if(mod(abs(j),2)==1) then 
                     msq(9,j,k)=(real(nf-2,dp)*Q(1)**2+2._dp*Q(2)**2)
     &                    *msq34_5(j,k)*sub34_5
                     msq(10,j,k)=(real(nf-2,dp)*Q(1)**2+2._dp*Q(2)**2)
     &                    *msq35_4(j,k)*sub35_4      
                  else
                     msq(9,j,k)=(real(nf-2,dp)*Q(1)**2+2._dp*Q(2)**2)
     &                    *msq34_5(j,k)*sub34_5
                     msq(10,j,k)=(real(nf-2,dp)*Q(1)**2+2._dp*Q(2)**2)
     &                    *msq35_4(j,k)*sub35_4 
                  endif
               endif
               
            else 
c               msq(n2,j,k)=msq(2,j,k)+(xn-1._dp/xn)
c     &              *(msq24_1(j,0)*sub24_1(gq)+msq24_1v(j,0)*sub24_1v)
               msq(n3,j,k)=msq(3,j,k)+(xn-1._dp/xn)
     &              *(msq15_2(0,k)*sub15_2(gq)+msq15_2v(0,k)*sub15_2v)
               if(frag) then      
                  msq(8,j,k)=Q(abs(k))**2*msq34_2(j,k)*sub34_2
                  msq(11,j,k)=Q(abs(j))**2*msq35_1(j,k)*sub35_1
               endif
                              
            endif
         endif

      enddo
      enddo

c--- remove contributions for qq~ -> photon+QQ~ with m(QQ~)<4*mQ^2  
c--- (to screen collinear divergence in g->QQ~)    
      if     (flav == 4) then
        if (2._dp*dot(p,4,5) < 4._dp*mcsq) then
          msq(:,4,-4)=zip
          msq(:,-4,4)=zip
        endif
      elseif (flav. eq. 5) then
        if (2._dp*dot(p,4,5) < 4._dp*mbsq) then
          msq(:,5,-5)=zip
          msq(:,-5,5)=zip
        endif
      endif
      
 
      return
      end
      
