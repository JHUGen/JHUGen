!====== fragmentation routine for trigam 
!==== based of gmgmjt tree C.Williams March 2013
      subroutine qqb_trigam_frag(p,msq) 
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'ewcharge.f' 
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'frag.f'
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf)
      double complex qqbg(2,2,2,2),qbqg(2,2,2,2)    
      double complex qgqb(2,2,2,2),qbgq(2,2,2,2)
      double complex gqqb(2,2,2,2),gqbq(2,2,2,2)
      double precision qqbg_sum,qbqg_sum
      double precision qgqb_sum,qbgq_sum
      double precision gqbq_sum,gqqb_sum
      integer h1,h2,h3,h4,j,k 
      double precision fac,statfac
      parameter(statfac=0.5d0)
      double precision fsq,D(0:5)
      common/D/D
      integer i 
      double precision conv_P0qqP0qgam,conv_P0qqDqgam
!$omp threadprivate(/D/)

      qbqg_sum=0d0 
      qqbg_sum=0d0 
      qbgq_sum=0d0 
      gqbq_sum=0d0 
      qgqb_sum=0d0 
      gqqb_sum=0d0 
      msq(:,:)=0d0 
      
      
      fsq=frag_scale**2
c---- Generate array D(j) corresponding to MCFM notation 0=gluon 1=up 2=down ....
        do i=0,5
           D(i)=0d0
           if     (fragset .eq. 'BFGset_I') then
              call get_frag(z_frag,fsq,1,i,D(i))   
           elseif (fragset .eq. 'BFGsetII') then  
              call get_frag(z_frag,fsq,2,i,D(i))   
           elseif (fragset .eq. 'GdRG__LO') then 
            call GGdR_frag(z_frag,i,D(i),0)
           elseif (fragset .eq. 'GdRG_NLO') then 
            call GGdR_frag(z_frag,i,D(i),1)
           else
              write(6,*) 'Unrecognized fragmentation set name: ',fragset
              stop        
           endif
        enddo

      fac=8d0*cf*xn*gsq*esq**2*statfac

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

            qqbg_sum=qqbg_sum+cdabs(qqbg(h1,h2,h3,h4))**2
            qbqg_sum=qbqg_sum+cdabs(qbqg(h1,h2,h3,h4))**2
            qgqb_sum=qgqb_sum+cdabs(qgqb(h1,h2,h3,h4))**2
            gqqb_sum=gqqb_sum+cdabs(gqqb(h1,h2,h3,h4))**2
            qbgq_sum=qbgq_sum+cdabs(qbgq(h1,h2,h3,h4))**2
            gqbq_sum=gqbq_sum+cdabs(gqbq(h1,h2,h3,h4))**2
            
         enddo
      enddo
      enddo
      enddo

      do j=-nf,nf
         do k=-nf,nf
            msq(j,k)=0d0
            if((j.ne.0).and.(k.ne.0).and.(j.ne.-k)) goto 20

            if((j.lt.0).and.(k.gt.0)) then 
               msq(j,k)=fac*aveqq*qbqg_sum*Q(k)**4*D(0)
            elseif((j.eq.0).and.(k.gt.0)) then 
               msq(j,k)=fac*aveqg*gqqb_sum*Q(k)**4*D(k)
            elseif((j.gt.0).and.(k.eq.0)) then 
               msq(j,k)=fac*aveqg*qgqb_sum*Q(j)**4*D(j)
            elseif((j.lt.0).and.(k.eq.0)) then 
               msq(j,k)=fac*aveqg*qbgq_sum*Q(j)**4*D(-j)
            elseif((j.eq.0).and.(k.lt.0)) then 
               msq(j,k)=fac*aveqg*gqbq_sum*Q(k)**4*D(-k)
            elseif((j.gt.0).and.(k.lt.0)) then 
               msq(j,k)=fac*aveqq*qqbg_sum*Q(j)**4*D(0)
            endif

 20         continue 
         enddo
      enddo
      
      return 
      end
