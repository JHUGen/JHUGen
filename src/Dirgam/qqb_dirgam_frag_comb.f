!========= new style of routine in which integrated frag and dipoles are combined together 
!==== C Williams April 2013
      subroutine qqb_dirgam_frag_combo(p,p_phys,msq) 
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'frag.f'
      include 'msqbits.f'
      include 'sprods_com.f'
      double precision p(mxpart,4),p_phys(mxpart,4)
      double precision msq_qcd(-nf:nf,-nf:nf),msq(-nf:nf,-nf:nf)
      double precision msq34_1(-nf:nf,-nf:nf),msqbits34_1(12)
      double precision msq34_1_swap(-nf:nf,-nf:nf),msqbits34_1_swap(12)
      double precision msq_qcd_swap(-nf:nf,-nf:nf)
      integer j,k,i
      double precision virt_dip,xl,dot,fsq 
      double precision aewo2pi,fi_gaq,ff_gaq
      double precision D(0:5),d_sum_q
      double precision smalla,smallb,smallc,ss,tt,uu,smalld
      double precision ga_ga,ag_ga,gq_gq,qg_gq,gg_gg,qa_gg,aq_gg
      double precision fac
      common/D/D
!$omp threadprivate(/D/)
   
      aewo2pi=esq/(fourpi*twopi)      
      fsq=frag_scale**2
      msq(:,:)=0d0 

!======= integerated dipoles 
      xl=dlog(-two*dot(p_phys,1,3)/fsq)
      virt_dip=+aewo2pi*(fi_gaq(z_frag,p_phys,xl,3,1,2))
!======== fragmentation functions 
c---- Generate array D(j) corresponding to MCFM notation 0=gluon 1=down 2=up ....
      do i=0,5
         D(i)=0d0
         if     (fragset .eq. 'BFGset_I') then
            call get_frag(z_frag,fsq,1,i,D(i))   
         elseif (fragset .eq. 'BFGsetII') then  
            call get_frag(z_frag,fsq,2,i,D(i))   
         elseif (fragset .eq. 'GdRG__LO') then 
            call GGdR_frag(z_frag,i,D(i),0) 
         else
            write(6,*) 'Unrecognized fragmentation set name: ',fragset
            stop        
         endif
      enddo

!============== matrix element functions 
      
      call qqb_2jet(p,msq34_1) 
      msqbits34_1(:)=msqbits(:) 
      call qqb_2jet_swap(p,msq34_1_swap)
!      call fill_dirgam_swap(7,msq34_1_swap)
      msqbits34_1_swap(:)=msqbits(:) 
!===== higher order pieces 
      call dotem(4,p,s)
      fac=gsq**2
     
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)
      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half

      gq_gq=-fac*aveqg*smallc(tt,ss,uu)
      ga_ga=-fac*aveqg*smallc(tt,uu,ss)
      qg_gq=-fac*aveqg*smallc(uu,ss,tt)
      ag_ga=-fac*aveqg*smallc(uu,tt,ss)

      gg_gg=fac*avegg*smalld(ss,tt,uu)*half



      do j=-nf,nf 
         do k=-nf,nf
            
!========= pieces with dipoles and fragmentation 
!========== factors of two arise from pieces with both 7 and 8 (34 and 35) singular regions, pieces with only 0
            if ((k .eq. 0).and.(j.ne.0)) then
               msq(j,k)=(Q(j)**2*virt_dip+D(abs(j)))*msq34_1(j,k)
            elseif((j.eq.0).and.(k.ne.0)) then 
               msq(j,k)=(Q(k)**2*virt_dip+D(abs(k)))*msq34_1(j,k)
            elseif((j.eq.0).and.(k.eq.0)) then 
               msq(j,k)=((3d0*Q(1)**2+2d0*Q(2)**2)*virt_dip
     &              +d_sum_q(0))*msq34_1(j,k)*two
            elseif((j.gt.0).and.(k.gt.0)) then 
               msq(j,k)=(Q(j)**2*virt_dip+D(abs(j)))*msq34_1(j,k)
     &             +(Q(k)**2*virt_dip+D(abs(k)))*msq34_1_swap(j,k)
            elseif((j.lt.0).and.(k.lt.0)) then 
               msq(j,k)=(Q(j)**2*virt_dip+D(abs(j)))*msq34_1(j,k)
     &             +(Q(k)**2*virt_dip+D(abs(k)))*msq34_1_swap(j,k)
            elseif((j.gt.0).and.(k.lt.0)) then 
               if(j.eq.-k) then 
               if  ((abs(j) .eq. 2) .or. (abs(j) .eq. 4)) then
                  msq(j,k)=
     &              ((Q(j)**2*virt_dip+D(j))*msqbits34_1(uub_uub)
     &             + (Q(2)**2*virt_dip)*msqbits34_1(uub_ccb)
     &             + (3d0*Q(1)**2*virt_dip)*msqbits34_1(uub_ddb)
     &               +d_sum_q(j)*msqbits34_1(uub_ddb))           
     &              +((Q(k)**2*virt_dip+D(-k))*msqbits34_1_swap(uub_uub)
     &             + (Q(2)**2*virt_dip)*msqbits34_1_swap(uub_ccb)
     &             + (3d0*Q(1)**2*virt_dip)*msqbits34_1_swap(uub_ddb)
     &               +d_sum_q(-k)*msqbits34_1_swap(uub_ddb))              
               else
                  msq(j,k)=
     &                ((Q(j)**2*virt_dip+D(j))*msqbits34_1(ddb_ddb)
     &               + (2d0*Q(2)**2*virt_dip)*msqbits34_1(ddb_uub)
     &               + (2d0*Q(1)**2*virt_dip)*msqbits34_1(ddb_ssb)
     &               +d_sum_q(j)*msqbits34_1(ddb_ssb))                          
     &              +((Q(k)**2*virt_dip+D(-k))*msqbits34_1_swap(ddb_ddb)
     &               + (2d0*Q(2)**2*virt_dip)*msqbits34_1_swap(ddb_uub)
     &               + (2d0*Q(1)**2*virt_dip)*msqbits34_1_swap(ddb_ssb)
     &               +d_sum_q(k)*msqbits34_1_swap(ddb_ssb))                          
               endif
               else
                  msq(j,k)=(Q(j)**2*virt_dip+D(j))*msq34_1(j,k)   
     &             +(Q(k)**2*virt_dip+D(abs(k)))*msq34_1_swap(j,k)
  
               endif
                  
            elseif((j.lt.0).and.(k.gt.0)) then 
               if(j.eq.-k) then 
               if  ((abs(j) .eq. 2) .or. (abs(j) .eq. 4)) then
                  msq(j,k)=
     &              ((Q(j)**2*virt_dip+D(-j))*msqbits34_1(ubu_uub)
     &             + (Q(2)**2*virt_dip)*msqbits34_1(ubu_ccb)
     &             + (3d0*Q(1)**2*virt_dip)*msqbits34_1(ubu_ddb)
     &               +d_sum_q(-j)*msqbits34_1(ubu_ddb))              
     &              +((Q(k)**2*virt_dip+D(k))*msqbits34_1_swap(ubu_uub)
     &             + (Q(2)**2*virt_dip)*msqbits34_1_swap(ubu_ccb)
     &             + (3d0*Q(1)**2*virt_dip)*msqbits34_1_swap(ubu_ddb)
     &               +d_sum_q(k)*msqbits34_1_swap(ubu_ddb))              
               else
                  msq(j,k)=
     &                ((Q(j)**2*virt_dip+D(-j))*msqbits34_1(dbd_ddb)
     &               + (2d0*Q(2)**2*virt_dip)*msqbits34_1(dbd_uub)
     &               + (2d0*Q(1)**2*virt_dip)*msqbits34_1(dbd_ssb)
     &               +d_sum_q(-j)*msqbits34_1(ddb_ssb))                          
     &             + ((Q(k)**2*virt_dip+D(k))*msqbits34_1_swap(dbd_ddb)
     &               + (2d0*Q(2)**2*virt_dip)*msqbits34_1_swap(dbd_uub)
     &               + (2d0*Q(1)**2*virt_dip)*msqbits34_1_swap(dbd_ssb)
     &               +d_sum_q(k)*msqbits34_1_swap(dbd_ssb))                        
               endif
            else
               msq(j,k)=(Q(j)**2*virt_dip+D(-j))*msq34_1(j,k)    
     &             +(Q(k)**2*virt_dip+D(abs(k)))*msq34_1_swap(j,k)
  
            endif
!============== these are the higher orer pieces which are proportional to the gluon fragmentation
            elseif((j.eq.0).and.(k.lt.0)) then 
               msq(j,k)=msq(j,k)+D(0)*ga_ga
            elseif((j.eq.0).and.(k.gt.0)) then 
               msq(j,k)=msq(j,k)+D(0)*gq_gq
            elseif((j.lt.0).and.(k.eq.0)) then 
               msq(j,k)=msq(j,k)+D(0)*ag_ga
            elseif((j.gt.0).and.(k.eq.0)) then 
               msq(j,k)=msq(j,k)+D(0)*qg_gq
            elseif((j.gt.0).and.(k.lt.0).and.(j.eq.-k)) then 
               msq(j,k)=msq(j,k)+D(0)*qa_gg*two
            elseif((j.lt.0).and.(k.gt.0).and.(j.eq.-k)) then 
               msq(j,k)=msq(j,k)+D(0)*aq_gg*two
            elseif((j.eq.0).and.(k.eq.0)) then 
               msq(j,k)=msq(j,k)+D(0)*gg_gg*two
            endif
         enddo
      enddo
      
      
      return 
      end
      

      double precision function d_sum_q(j)
!-----sum over i=1,5 i!=j   
      implicit none 
      double precision D(0:5) 
      common/D/D
      integer i,j
!$omp threadprivate(/D/)

      d_sum_q=0d0
      do i=1,5 
         if(i.ne.j) then 
            d_sum_q=d_sum_q+D(i) 
         endif
      enddo 
      return 
      end
