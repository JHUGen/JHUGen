      subroutine qqb_zaj_frag(p,msq)
 !------- FRAGMENTATION CONTRIBUTION, Always fragment p5 (=j) , allow phase space integration to symmetrize. 

      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'frag.f'
      integer j,k,i
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      double precision D(0:5),d_sum_q,fsq
      common/D/D
      double precision mqq_z2j(0:2,fn:nf,fn:nf),msq0,msq1,msq2
      double precision msqx_z2j(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      double precision msqx_cs_z2j(0:2,-nf:nf,-nf:nf)
      double precision msq_z2j(-nf:nf,-nf:nf) 
      integer m,n,icol
      double precision qcd_z2j(-nf:nf,-nf:nf,-nf:nf,-nf:nf)
!$omp threadprivate(/D/)

      fsq=frag_scale**2
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

      
!====== NEED TO SETUP equivlient of qqii_ii etc for 
!====== this process
!====== Underlying QCD matrix element 
      call qqb_z2jetx(p,msq_z2j,mqq_z2j,msqx_z2j,msqx_cs_z2j) 
!----- not naming follows z2jet process, we will work with msqx_z2j 

     
!====== We do not need individual colour compenents so define new matrices which are sums over colour 

      do i=-nf,nf
         do j=-nf,nf
            do k=-nf,nf
               do m=-nf,nf 
                  qcd_z2j(i,j,k,m)=0d0 
                  do icol=0,2 
               qcd_z2j(i,j,k,m)=qcd_z2j(i,j,k,m)+msqx_z2j(icol,i,j,k,m)
                   enddo
                enddo
             enddo
          enddo
       enddo



       do j=-nf,nf
       do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
C--qq      
      if ((j .gt. 0) .and. (k .gt. 0)) then
          if (j .eq. k) then

!---- p5 = j 
            msq(j,k)=two*qcd_z2j(j,j,j,j)*D(abs(j))
          else 
            msq(j,k)=qcd_z2j(j,k,j,k)*D(abs(j))+qcd_z2j(j,k,k,j)*D(k)
          endif

C--qa      
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
             
             msq(j,k)=two*qcd_z2j(j,k,j,k)*D(j)
             do m=1,5 
                if(m.ne.j) then 
                   msq(j,k)=msq(j,k)+two*qcd_z2j(j,k,m,-m)*D(m)
                endif
             enddo
             msq(j,k)=msq(j,k)+qcd_z2j(j,k,0,0)*D(0)*two
          else
!-------p3=j 
             msq(j,k)=qcd_z2j(j,k,j,k)*D(j)+qcd_z2j(j,k,k,j)*D(abs(k))
          endif

C--aa      
       elseif ((j .lt. 0) .and. (k .lt. 0)) then
          if (j .eq. k) then
             msq(j,k)=two*qcd_z2j(j,k,j,k)*D(abs(j))
          else
             msq(j,k)=qcd_z2j(j,k,j,k)*D(abs(j))
     &            +qcd_z2j(j,k,k,j)*D(abs(k))
          endif
          
C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=two*(qcd_z2j(j,k,j,k))*D(abs(j))
            do m=1,5 
               if(m.ne.-j) then 
                  msq(j,k)=msq(j,k)+two*qcd_z2j(j,k,m,-m)*D(m)
               endif
            enddo
            msq(j,k)=msq(j,k)+two*(qcd_z2j(j,k,0,0)*D(0))
          else
            msq(j,k)=qcd_z2j(j,k,j,k)*D(abs(j))+qcd_z2j(j,k,k,j)*D(k)
          endif

C--qg_qg      
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=qcd_z2j(j,k,j,0)*D(j)+qcd_z2j(j,k,0,j)*D(0)
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=qcd_z2j(j,k,j,0)*D(abs(j))+qcd_z2j(j,k,0,j)*D(0)
C--gq_gq      
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=qcd_z2j(j,k,k,0)*D(k)+qcd_z2j(j,k,0,k)*D(0)
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=qcd_z2j(j,k,k,0)*D(abs(k))+qcd_z2j(j,k,0,k)*D(0)
C--gg      
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=qcd_z2j(j,k,0,0)*D(0)
            do m=1,5
               msq(j,k)=msq(j,k)+qcd_z2j(0,0,m,-m)*D(m)
            enddo
      endif

      enddo
      enddo


      return
      end



