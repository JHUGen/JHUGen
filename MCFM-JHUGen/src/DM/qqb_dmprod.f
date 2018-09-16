      subroutine qqb_dm_prod(p,msq) 
      implicit none
      include 'types.f'
            
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'dm_params.f'
      
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
!------- first index = helicity of fermion line, 
!-------- second index, helicity of xi(p3),3rd xi~(p4) massive allows
!-------- helicity violation.
      complex(dp):: qqb_dmamp(2,2,2)
      complex(dp):: qbq_dmamp(2,2,2)
      real(dp):: qbq_sum,qqb_sum
      complex(dp):: prop
    
      real(dp):: s34,fac
      integer:: h1,h2,h3 
      integer:: j,k,nu

      fac=one/dm_lam**4*xn

      if(dm_mediator=='vector') then 
         call dmprod_vec(p,1,2,3,4,qqb_dmamp)  
         call dmprod_vec(p,2,1,3,4,qbq_dmamp) 
      elseif(dm_mediator=='axvect') then 
!         call dmprod_avec(p,1,2,3,4,qqb_dmamp)
!         call dmprod_avec(p,2,1,3,4,qbq_dmamp)        
      elseif(dm_mediator=='scalar') then 
         call dmprod_scal(p,1,2,3,4,qqb_dmamp)
         call dmprod_scal(p,2,1,3,4,qbq_dmamp)
      elseif(dm_mediator=='pseudo') then 
!         call dmprod_ps(p,1,2,3,4,qqb_dmamp)
!         call dmprod_ps(p,2,1,3,4,qbq_dmamp) 
      endif
      
      s34=0d0 
      do nu=1,3
         s34=s34-(p(3,nu)+p(4,nu))**2 
      enddo
      nu=4
      s34=s34+(p(3,nu)+p(4,nu))**2 
     
      medwidth=1d0
!      prop=cplx2(one/(s34-medmass**2),medwidth*medmass) 
!      prop=one/medmass**2
      prop=one
      qbq_sum=0d0 
      qqb_sum=0d0 
      do h1=1,2 
         do h2=1,2 
            do h3=1,2
            qqb_sum=qqb_sum+abs(prop*qqb_dmamp(h1,h2,h3))**2
            qbq_sum=qbq_sum+abs(prop*qbq_dmamp(h1,h2,h3))**2
         enddo
      enddo
      enddo

     
      do j=-nf,nf 
         do k=-nf,nf
            msq(j,k)=0d0 
         enddo
      enddo


            
      do j=-nf,nf
           if((j<0)) then
              msq(j,-j)=qbq_sum*aveqq*fac
           elseif((j>0)) then 
              msq(j,-j)=qqb_sum*aveqq*fac
           elseif(j==0) then 
              msq(0,0)=0d0 
           endif
        enddo
!      enddo


!        write(6,*) msq
!        pause

      return 
      end 

      subroutine gen_masslessvecs(pin,pout,a,b) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'dm_params.f'
      real(dp):: pin(mxpart,4),pout(mxpart,4) 
      integer:: i,nu
      integer:: a,b
      real(dp):: beta,sab,opb,omb

      sab=0d0 
      do nu=1,3
         sab=sab-(pin(a,nu)+pin(b,nu))**2 
      enddo
      nu=4
      sab=sab+(pin(a,nu)+pin(b,nu))**2 
      
!      write(6,*) ' in genmas', sab

      beta=1d0-4d0*xmass**2/sab
      beta=sqrt(beta)
      opb=(one+beta)
      omb=one-beta
      
      do i=1,mxpart
         do nu=1,4
            if(i==a) then 
               pout(a,nu)=0.5d0*(opb/beta*pin(a,nu)-omb/beta*pin(b,nu))
            elseif(i==b) then 
               pout(b,nu)=0.5d0*(opb/beta*pin(b,nu)-omb/beta*pin(a,nu))
            else
               pout(i,nu)=pin(i,nu)
            endif
      enddo
      enddo
      return 
      end 
 

      subroutine dmprod_scal(p,i1,i2,i3,i4,amp) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4),q(mxpart,4) 
      complex(dp):: amp(2,2,2),amp_dec(2,2)
      integer:: i1,i2,i3,i4 
      integer:: h1,h2
!--------- generate phase space with massless vectors
      call gen_masslessvecs(p,q,i3,i4) 
!--------- generate spinors 
      call spinoru(4,q,za,zb) 
!--------- call generic dm decay 
      call scalar_dm(i3,i4,za,zb,amp_dec)
    

      do h1=1,2
         do h2=1,2 
            amp(1,h1,h2)=-za(i1,i2)*amp_dec(h1,h2) 
            amp(2,h1,h2)=zb(i2,i1)*amp_dec(h1,h2) 
         enddo
      enddo

      return 
      end
      


      subroutine dmprod_vec(p,i1,i2,i3,i4,amp) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f'
      include 'dm_params.f'
      real(dp):: p(mxpart,4),q(mxpart,4)
      complex(dp):: amp(2,2,2) 
!      complex(dp):: jpp(4),jpm(4),jmp(4),jmm(4) 
      integer:: i1,i2,i3,i4 
      complex(dp):: Nmm,Npp
!----- generate phase space with massless vectors       
      call gen_masslessvecs(p,q,i3,i4)
!----- generate spinors 
      call spinoru(4,q,za,zb)


      Nmm=xmass/zb(i3,i4)
      Npp=xmass/za(i3,i4)

      amp(1,2,1)=2d0*za(i1,i4)*zb(i2,i3)
      amp(1,1,2)=2d0*za(i1,i3)*zb(i2,i4)
      amp(1,1,1)=-2d0*za(i1,i4)*zb(i4,i2)*Nmm
      amp(1,2,2)=-2d0*za(i1,i4)*zb(i4,i2)*Npp
      amp(2,2,1)=2d0*za(i2,i4)*zb(i1,i3)
      amp(2,1,2)=2d0*za(i2,i3)*zb(i1,i4)
      amp(2,1,1)=-2d0*za(i2,i4)*zb(i4,i1)
      amp(2,1,2)=2d0*za(i2,i3)*zb(i1,i4)
      amp(2,1,1)=-2d0*za(i2,i4)*zb(i4,i1)*Nmm
      amp(2,2,2)=-2d0*za(i2,i4)*zb(i4,i1)*Npp
      

      return 
      end 
      
