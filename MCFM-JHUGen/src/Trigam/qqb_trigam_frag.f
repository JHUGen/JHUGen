!====== fragmentation routine for trigam 
!==== based of gmgmjt tree C.Williams March 2013
      subroutine qqb_trigam_frag(p,msq) 
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'ewcharge.f' 
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'frag.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      complex(dp):: qqbg(2,2,2,2),qbqg(2,2,2,2)    
      complex(dp):: qgqb(2,2,2,2),qbgq(2,2,2,2)
      complex(dp):: gqqb(2,2,2,2),gqbq(2,2,2,2)
      real(dp):: qqbg_sum,qbqg_sum
      real(dp):: qgqb_sum,qbgq_sum
      real(dp):: gqbq_sum,gqqb_sum
      integer:: h1,h2,h3,h4,j,k 
      real(dp):: fac,statfac
      parameter(statfac=0.5_dp)
      real(dp):: fsq,D(0:5)
      common/D/D
      integer:: i 
      real(dp):: conv_P0qqP0qgam,conv_P0qqDqgam
!$omp threadprivate(/D/)

      qbqg_sum=0._dp 
      qqbg_sum=0._dp 
      qbgq_sum=0._dp 
      gqbq_sum=0._dp 
      qgqb_sum=0._dp 
      gqqb_sum=0._dp 
      msq(:,:)=0._dp 
      
      
      fsq=frag_scale**2
c---- Generate array D(j) corresponding to MCFM notation 0=gluon 1=up 2=down ....
        do i=0,5
           D(i)=0._dp
           if     (fragset == 'BFGset_I') then
              call get_frag(z_frag,fsq,1,i,D(i))   
           elseif (fragset == 'BFGsetII') then  
              call get_frag(z_frag,fsq,2,i,D(i))   
           elseif (fragset == 'GdRG__LO') then 
            call GGdR_frag(z_frag,i,D(i),0)
           elseif (fragset == 'GdRG_NLO') then 
            call GGdR_frag(z_frag,i,D(i),1)
           else
              write(6,*) 'Unrecognized fragmentation set name: ',fragset
              stop        
           endif
        enddo

      fac=8._dp*cf*xn*gsq*esq**2*statfac

      call spinoru(5,p,za,zb)
      call amp_lord_gmgmjt(1,2,5,3,4,za,zb,qqbg)
      call amp_lord_gmgmjt(2,1,5,3,4,za,zb,qbqg)
      call amp_lord_gmgmjt(1,5,2,3,4,za,zb,qgqb)
      call amp_lord_gmgmjt(2,5,1,3,4,za,zb,gqqb)
      call amp_lord_gmgmjt(5,1,2,3,4,za,zb,qbgq)
      call amp_lord_gmgmjt(5,2,1,3,4,za,zb,gqbq)

      do h1 =1,2
         do h2 =1,2
         do h3 =1,2
         do h4 =1,2

            qqbg_sum=qqbg_sum+abs(qqbg(h1,h2,h3,h4))**2
            qbqg_sum=qbqg_sum+abs(qbqg(h1,h2,h3,h4))**2
            qgqb_sum=qgqb_sum+abs(qgqb(h1,h2,h3,h4))**2
            gqqb_sum=gqqb_sum+abs(gqqb(h1,h2,h3,h4))**2
            qbgq_sum=qbgq_sum+abs(qbgq(h1,h2,h3,h4))**2
            gqbq_sum=gqbq_sum+abs(gqbq(h1,h2,h3,h4))**2
            
         enddo
      enddo
      enddo
      enddo

      do j=-nf,nf
         do k=-nf,nf
            msq(j,k)=0._dp
            if((j.ne.0).and.(k.ne.0).and.(j.ne.-k)) goto 20

            if((j<0).and.(k>0)) then 
               msq(j,k)=fac*aveqq*qbqg_sum*Q(k)**4*D(0)
            elseif((j==0).and.(k>0)) then 
               msq(j,k)=fac*aveqg*gqqb_sum*Q(k)**4*D(k)
            elseif((j>0).and.(k==0)) then 
               msq(j,k)=fac*aveqg*qgqb_sum*Q(j)**4*D(j)
            elseif((j<0).and.(k==0)) then 
               msq(j,k)=fac*aveqg*qbgq_sum*Q(j)**4*D(-j)
            elseif((j==0).and.(k<0)) then 
               msq(j,k)=fac*aveqg*gqbq_sum*Q(k)**4*D(-k)
            elseif((j>0).and.(k<0)) then 
               msq(j,k)=fac*aveqq*qqbg_sum*Q(j)**4*D(0)
            endif

 20         continue 
         enddo
      enddo
      
      return 
      end
