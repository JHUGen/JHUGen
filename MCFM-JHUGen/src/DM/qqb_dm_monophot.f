
      subroutine qqb_dm_monophot(p,msq) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'dm_params.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'first.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      complex(dp):: qqbg(2,2,2,2),qbqg(2,2,2,2)
      real(dp):: qqbg_sum(nf),qbqg_sum(nf)
      real(dp):: fac 
      integer:: h2,h3,h4 
      integer:: j,k
      complex(dp):: cprop
      real(dp):: propsq
  
      real(dp):: s34
        logical:: check_QED 
      common/check_QED/check_QED 
!$omp threadprivate(/check_QED/)

      check_QED=.false. 

      if(check_QED) then 
!------ make factor into photon propogator for testing 
         s34=2d0*(p(3,4)*p(4,4)-p(3,3)*p(4,3)
     &        -p(3,2)*p(4,2)-p(3,1)*p(4,1))
         dm_lam=sqrt(s34/esq)
         do j=1,nf 
            dmL(j)=Q(j)
            dmR(j)=Q(j) 
         enddo
      endif

      if(effective_th) then 
         fac=8d0/dm_lam**4*esq*xn
      else
!--------full theory => for V,A and PS is simply g_dmx**2*g_dmq**2/prop(34)**2 
!--------for S and G need to be more careful, but G happens elsewhere and S can be done 
!--------in special routine
         s34=(p(3,4)+p(4,4))**2
     &        -(p(3,1)+p(4,1))**2
     &        -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2         
         cprop=cone/cplx2((s34-medmass**2),medmass*medwidth)
         propsq=abs(cprop)**2 
         fac=8d0*esq*xn*propsq*g_dmq**2*g_dmx**2
      endif
      
      

   
      msq(:,:)=0
!------- wrapper for separte operators
      if(dm_mediator=='vector') then 
         call qqb_dm_monojet_Vamps(p,1,5,2,3,4,qqbg)  
         call qqb_dm_monojet_Vamps(p,2,5,1,3,4,qbqg)         
!         call qqb_dm_monojet_vecamps(p,1,5,2,3,4,qqbg)  
!         call qqb_dm_monojet_vecamps(p,2,5,1,3,4,qbqg)   
      elseif(dm_mediator=='axvect') then 
!         call qqb_dm_monojet_Axamps(p,1,2,5,3,4,qgqb)  
!         call qqb_dm_monojet_Axamps(p,5,2,1,3,4,qbgq)
         call qqb_dm_monojet_Axamps(p,1,5,2,3,4,qqbg)  
         call qqb_dm_monojet_Axamps(p,2,5,1,3,4,qbqg)
!        call qqb_dm_monojet_Axamps(p,2,1,5,3,4,gqqb)  
!        call qqb_dm_monojet_Axamps(p,5,1,2,3,4,gqbq)  
      elseif(dm_mediator=='scalar') then 
         call qqb_dm_monojet_Samps(p,1,5,2,3,4,qqbg)  
         call qqb_dm_monojet_Samps(p,2,5,1,3,4,qbqg)
         fac=fac/4d0
      elseif(dm_mediator=='pseudo') then 
         call qqb_dm_monojet_PSamps(p,1,5,2,3,4,qqbg)  
         call qqb_dm_monojet_PSamps(p,2,5,1,3,4,qbqg)
         fac=fac/4d0
      endif

   
      qqbg_sum(:)=0d0
      qbqg_sum(:)=0d0
      do j=1,nf
         do h2=1,2
            do h3=1,2 
               do h4=1,2 
                  qqbg_sum(j)=qqbg_sum(j)
     &                 +abs(dmL(j)*qqbg(1,h2,h3,h4))**2
     &                 +abs(dmR(j)*qqbg(2,h2,h3,h4))**2
                  qbqg_sum(j)=qbqg_sum(j)
     &                 +abs(dmL(j)*qbqg(1,h2,h3,h4))**2
     &                 +abs(dmR(j)*qbqg(2,h2,h3,h4))**2
            enddo
         enddo
      enddo
      enddo


      do j=-nf,nf 
         do k=-nf,nf 
           if((j .ne. -k)) goto 20

           if((j==0).and.(k==0)) then 
              msq(j,k)=0d0 
           elseif((j>0).and.(k<0)) then 
              msq(j,k)=aveqq*qqbg_sum(j)*fac*Q(j)**2 
           elseif((j<0).and.(k>0)) then 
              msq(j,k)=aveqq*qbqg_sum(k)*fac*Q(k)**2
           endif


 20        continue 
        enddo
      enddo

      return 
      end 
