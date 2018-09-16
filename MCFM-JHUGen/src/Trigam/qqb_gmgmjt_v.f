!====== LO ROUTINE FOR GAM GAM JET using amplitudes

      subroutine qqb_gmgmjt_v(p,msq) 
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'ewcharge.f' 
      include 'ewcouple.f'
      include 'qcdcouple.f' 
      include 'epinv.f'
      include 'nflav.f'
      include 'scheme.f'
      double precision subuv
      double precision p(mxpart,4),msq(-nf:nf,-nf:nf)
      double precision msq0(-nf:nf,-nf:nf)
      double complex qqbg_lo(2,2,2,2),qbqg_lo(2,2,2,2)    
      double complex qgqb_lo(2,2,2,2),qbgq_lo(2,2,2,2)
      double complex gqqb_lo(2,2,2,2),gqbq_lo(2,2,2,2)
      double complex qqbg_v(2,2,2,2),qbqg_v(2,2,2,2)    
      double complex qgqb_v(2,2,2,2),qbgq_v(2,2,2,2)
      double complex gqqb_v(2,2,2,2),gqbq_v(2,2,2,2)
      double complex qqbg_nf(2,2,2,2),qbqg_nf(2,2,2,2)    
      double complex qgqb_nf(2,2,2,2),qbgq_nf(2,2,2,2)
      double complex gqqb_nf(2,2,2,2),gqbq_nf(2,2,2,2)
      double precision qqbg_sum,qbqg_sum
      double precision qgqb_sum,qbgq_sum
      double precision gqbq_sum,gqqb_sum
      double precision qqbg_nfs,qbqg_nfs
      double precision qgqb_nfs,qbgq_nfs
      double precision gqbq_nfs,gqqb_nfs
      integer h1,h2,h3,h4,j,k 
      double precision fac,statfac
      parameter(statfac=0.5d0)
      double complex test,virt_gmgmjt_GaMHV
      double precision faclo,nf_fac

      qbqg_sum=0d0 
      qqbg_sum=0d0 
      qbgq_sum=0d0 
      gqbq_sum=0d0 
      qgqb_sum=0d0 
      gqqb_sum=0d0 

      qbqg_nfs=0d0 
      qqbg_nfs=0d0 
      qbgq_nfs=0d0 
      gqbq_nfs=0d0 
      qgqb_nfs=0d0 
      gqqb_nfs=0d0 

      msq(:,:)=0d0 
      scheme='dred'

      faclo=8d0*cf*xn*gsq*esq**2*statfac
      fac=faclo*xn*ason2pi/2d0

!----UV counterterm contains the finite renormalization to arrive
!----at MS bar scheme.      
!====== for now call born separately. probably want to adjust this later 
!====== to make code more efficient
      call qqb_gmgmjt(p,msq0)
      subuv=ason2pi*xn*(epinv*(11d0-2d0*dble(nflav)/xn)-1d0)/6d0
      call spinoru(5,p,za,zb)
      call amp_virt_gmgmjt(1,2,3,4,5,za,zb,qqbg_v)
      call amp_lord_gmgmjt(1,2,5,3,4,za,zb,qqbg_lo)
      call amp_lord_gmgmjt(2,1,5,3,4,za,zb,qbqg_lo)
      call amp_lord_gmgmjt(1,5,2,3,4,za,zb,qgqb_lo)
      call amp_lord_gmgmjt(2,5,1,3,4,za,zb,gqqb_lo)
      call amp_lord_gmgmjt(5,1,2,3,4,za,zb,qbgq_lo)
      call amp_lord_gmgmjt(5,2,1,3,4,za,zb,gqbq_lo)

      call amp_virt_gmgmjt(1,2,5,3,4,za,zb,qqbg_v)
      call amp_virt_gmgmjt(2,1,5,3,4,za,zb,qbqg_v)
      call amp_virt_gmgmjt(1,5,2,3,4,za,zb,qgqb_v)
      call amp_virt_gmgmjt(2,5,1,3,4,za,zb,gqqb_v)
      call amp_virt_gmgmjt(5,1,2,3,4,za,zb,qbgq_v)
      call amp_virt_gmgmjt(5,2,1,3,4,za,zb,gqbq_v)

!======= nf loops
      call amp_virt_nf_gmgmjt(1,2,5,3,4,za,zb,qqbg_nf)
      call amp_virt_nf_gmgmjt(2,1,5,3,4,za,zb,qbqg_nf)
      call amp_virt_nf_gmgmjt(1,5,2,3,4,za,zb,qgqb_nf)
      call amp_virt_nf_gmgmjt(2,5,1,3,4,za,zb,gqqb_nf)
      call amp_virt_nf_gmgmjt(5,1,2,3,4,za,zb,qbgq_nf)
      call amp_virt_nf_gmgmjt(5,2,1,3,4,za,zb,gqbq_nf)
!====== prefactor for nf loops
      nf_fac=-fac/xn
     

      do h1 =1,2
         do h2 =1,2
         do h3 =1,2
         do h4 =1,2

            qqbg_sum=qqbg_sum+
     & Dble(Dconjg(qqbg_v(h1,h2,h3,h4))*qqbg_lo(h1,h2,h3,h4))
     & +Dble(Dconjg(qqbg_lo(h1,h2,h3,h4))*qqbg_v(h1,h2,h3,h4))
            qqbg_nfs=qqbg_nfs+     
     & Dble(Dconjg(qqbg_nf(h1,h2,h3,h4))*qqbg_lo(h1,h2,h3,h4))
     & +Dble(Dconjg(qqbg_lo(h1,h2,h3,h4))*qqbg_nf(h1,h2,h3,h4))

            qbqg_sum=qbqg_sum+
     & Dble(Dconjg(qbqg_v(h1,h2,h3,h4))*qbqg_lo(h1,h2,h3,h4))
     & +Dble(Dconjg(qbqg_lo(h1,h2,h3,h4))*qbqg_v(h1,h2,h3,h4))
            qbqg_nfs=qbqg_nfs+
     & Dble(Dconjg(qbqg_nf(h1,h2,h3,h4))*qbqg_lo(h1,h2,h3,h4))
     & +Dble(Dconjg(qbqg_lo(h1,h2,h3,h4))*qbqg_nf(h1,h2,h3,h4))
            
            qgqb_sum=qgqb_sum+
     & Dble(Dconjg(qgqb_v(h1,h2,h3,h4))*qgqb_lo(h1,h2,h3,h4))
     & +Dble(Dconjg(qgqb_lo(h1,h2,h3,h4))*qgqb_v(h1,h2,h3,h4))
            qgqb_nfs=qgqb_nfs+
     & Dble(Dconjg(qgqb_nf(h1,h2,h3,h4))*qgqb_lo(h1,h2,h3,h4))
     & +Dble(Dconjg(qgqb_lo(h1,h2,h3,h4))*qgqb_nf(h1,h2,h3,h4))

            gqqb_sum=gqqb_sum+
     & Dble(Dconjg(gqqb_v(h1,h2,h3,h4))*gqqb_lo(h1,h2,h3,h4))
     & +Dble(Dconjg(gqqb_lo(h1,h2,h3,h4))*gqqb_v(h1,h2,h3,h4))
            gqqb_nfs=gqqb_nfs+
     & Dble(Dconjg(gqqb_nf(h1,h2,h3,h4))*gqqb_lo(h1,h2,h3,h4))
     & +Dble(Dconjg(gqqb_lo(h1,h2,h3,h4))*gqqb_nf(h1,h2,h3,h4))
            

            qbgq_sum=qbgq_sum+
     & Dble(Dconjg(qbgq_v(h1,h2,h3,h4))*qbgq_lo(h1,h2,h3,h4))
     & +Dble(Dconjg(qbgq_lo(h1,h2,h3,h4))*qbgq_v(h1,h2,h3,h4))
            qbgq_nfs=qbgq_nfs+
     & Dble(Dconjg(qbgq_nf(h1,h2,h3,h4))*qbgq_lo(h1,h2,h3,h4))
     & +Dble(Dconjg(qbgq_lo(h1,h2,h3,h4))*qbgq_nf(h1,h2,h3,h4))

            gqbq_sum=gqbq_sum+
     & Dble(Dconjg(gqbq_v(h1,h2,h3,h4))*gqbq_lo(h1,h2,h3,h4))
     & +Dble(Dconjg(gqbq_lo(h1,h2,h3,h4))*gqbq_v(h1,h2,h3,h4))
            gqbq_nfs=gqbq_nfs+
     & Dble(Dconjg(gqbq_nf(h1,h2,h3,h4))*gqbq_lo(h1,h2,h3,h4))
     & +Dble(Dconjg(gqbq_lo(h1,h2,h3,h4))*gqbq_nf(h1,h2,h3,h4))

         enddo
      enddo
      enddo
      enddo

      do j=-nf,nf
         do k=-nf,nf
            msq(j,k)=0d0
            if((j.ne.0).and.(k.ne.0).and.(j.ne.-k)) goto 20

            if((j.lt.0).and.(k.gt.0)) then 
               msq(j,k)=fac*aveqq*qbqg_sum*Q(k)**4-subuv*msq0(j,k)
     &              +nf_fac*aveqq*qbqg_nfs*(2d0*Q(2)**4+3d0*Q(1)**4)
            elseif((j.eq.0).and.(k.gt.0)) then 
               msq(j,k)=fac*aveqg*gqqb_sum*Q(k)**4-subuv*msq0(j,k)
     &              +nf_fac*aveqg*gqqb_nfs*(2d0*Q(2)**4+3d0*Q(1)**4)
            elseif((j.gt.0).and.(k.eq.0)) then 
               msq(j,k)=fac*aveqg*qgqb_sum*Q(j)**4-subuv*msq0(j,k)
     &              +nf_fac*aveqg*qgqb_nfs*(2d0*Q(2)**4+3d0*Q(1)**4)
            elseif((j.lt.0).and.(k.eq.0)) then 
               msq(j,k)=fac*aveqg*qbgq_sum*Q(j)**4-subuv*msq0(j,k)
     &              +nf_fac*aveqg*qbgq_nfs*(2d0*Q(2)**4+3d0*Q(1)**4) 
            elseif((j.eq.0).and.(k.lt.0)) then 
               msq(j,k)=fac*aveqg*gqbq_sum*Q(k)**4-subuv*msq0(j,k)
     &              +nf_fac*aveqg*gqbq_nfs*(2d0*Q(2)**4+3d0*Q(1)**4)
            elseif((j.gt.0).and.(k.lt.0)) then 
               msq(j,k)=fac*aveqq*qqbg_sum*Q(j)**4-subuv*msq0(j,k)
     &              +nf_fac*aveqq*qqbg_nfs*(2d0*Q(2)**4+3d0*Q(1)**4)
            endif

 20         continue 
         enddo
      enddo
      return 
      end


      subroutine amp_virt_gmgmjt(i1,i2,i3,i4,i5,za,zb,amp)
      implicit none 
      include 'constants.f' 
      include 'zprods_decl.f' 
      integer i1,i2,i3,i4,i5
      double complex amp(2,2,2,2),amp_lc(2,2,2,2)
      double complex amp_slc(2,2,2,2),amp_nf(2,2,2,2)
      integer h1,h2,h3,h4 
!==== call lc amplitudes 
      call amp_virt_lc_gmgmjt(i1,i2,i3,i4,i5,za,zb,amp_lc)
!====== subleading 
      call amp_virt_slc_gmgmjt(i1,i2,i3,i4,i5,za,zb,amp_slc) 
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  amp(h1,h2,h3,h4)=-amp_lc(h1,h2,h3,h4)
     &                 -one/xnsq*amp_slc(h1,h2,h3,h4)
               enddo
            enddo
         enddo
      enddo

      return 
      end



      subroutine amp_virt_lc_gmgmjt(i1,i2,i3,i4,i5,za,zb,amp)
      implicit none 
      include 'constants.f' 
      include 'zprods_decl.f' 
      integer i1,i2,i3,i4,i5
      double complex amp(2,2,2,2) 
      double complex amp_2gam1g,virt_gmgmjt_GMHV,virt_gmgmjt_GaMHV
      double complex virt_gmgmjt_gammaMHV,virt_gmgmjt_gluonMHV
!===== default is for gluon MHV amplitude 

c--- DEBUG ONLY
c      write(6,*)
c      amp(1,2,2,1)=virt_gmgmjt_GMHV(i1,i2,i3,i4,i5,za,zb)
c      amp(1,2,2,1)=virt_gmgmjt_gluonMHV(i1,i2,i3,i4,i5,za,zb)
c      write(6,*)
c      pause
c--- DEBUG ONLY

      amp(:,:,:,:)=czip
      amp(1,1,1,1)=czip
      amp(2,2,2,2)=czip
      amp(2,1,1,1)=czip
      amp(1,2,2,2)=czip

!===== MHV amplitudes 
      amp(1,1,2,2)=virt_gmgmjt_gluonMHV(i1,i2,i3,i4,i5,za,zb)
      amp(1,2,1,2)=virt_gmgmjt_gammaMHV(i1,i2,i3,i5,i4,za,zb)
      amp(1,2,2,1)=virt_gmgmjt_gammaMHV(i1,i2,i3,i4,i5,za,zb)

!===== line reversal 
      amp(2,1,2,2)=virt_gmgmjt_gluonMHV(i2,i1,i3,i4,i5,za,zb)
      amp(2,2,1,2)=virt_gmgmjt_gammaMHV(i2,i1,i3,i5,i4,za,zb)
      amp(2,2,2,1)=virt_gmgmjt_gammaMHV(i2,i1,i3,i4,i5,za,zb)

!====== conjugation 
      amp(2,2,1,1)=-virt_gmgmjt_gluonMHV(i1,i2,i3,i4,i5,zb,za)
      amp(2,1,2,1)=-virt_gmgmjt_gammaMHV(i1,i2,i3,i5,i4,zb,za)
      amp(2,1,1,2)=-virt_gmgmjt_gammaMHV(i1,i2,i3,i4,i5,zb,za)
!====== conjugation + line reversal 
      amp(1,2,1,1)=-virt_gmgmjt_gluonMHV(i2,i1,i3,i4,i5,zb,za)
      amp(1,1,2,1)=-virt_gmgmjt_gammaMHV(i2,i1,i3,i5,i4,zb,za)
      amp(1,1,1,2)=-virt_gmgmjt_gammaMHV(i2,i1,i3,i4,i5,zb,za)
      

      return 
      end 

      subroutine amp_virt_nf_gmgmjt(i1,i2,i3,i4,i5,za,zb,amp)
      implicit none 
      include 'constants.f' 
      include 'zprods_decl.f' 
      integer i1,i2,i3,i4,i5
      double complex amp(2,2,2,2) 
      double complex amp_2gam1g
      double complex virt_gmgmjt_nfGaMHV,virt_gmgmjt_nfgammaMHV
!===== default is for gluon MHV amplitude 

c--- DEBUG ONLY
c      write(6,*)
c      amp(1,2,2,1)=virt_gmgmjt_nfGaMHV(i1,i2,i3,i4,i5,za,zb)
c      amp(1,2,2,1)=virt_gmgmjt_nfgammaMHV(i1,i2,i3,i4,i5,za,zb)
c      write(6,*)
c      pause
c--- DEBUG ONLY

      amp(:,:,:,:)=czip
      amp(1,1,1,1)=czip
      amp(2,2,2,2)=czip
      amp(2,1,1,1)=czip
      amp(1,2,2,2)=czip
!===== MHV amplitudes 
      amp(1,2,2,1)=virt_gmgmjt_nfgammaMHV(i1,i2,i3,i4,i5,za,zb)
      amp(1,2,1,2)=virt_gmgmjt_nfgammaMHV(i1,i2,i3,i5,i4,za,zb)
      amp(1,1,2,2)=virt_gmgmjt_nfgammaMHV(i1,i2,i4,i5,i3,za,zb)

!===== line reversal 
      amp(2,2,2,1)=virt_gmgmjt_nfgammaMHV(i2,i1,i3,i4,i5,za,zb)
      amp(2,2,1,2)=virt_gmgmjt_nfgammaMHV(i2,i1,i3,i5,i4,za,zb)
      amp(2,1,2,2)=virt_gmgmjt_nfgammaMHV(i2,i1,i4,i5,i3,za,zb)

!====== conjugation 
      amp(2,1,1,2)=virt_gmgmjt_nfgammaMHV(i1,i2,i3,i4,i5,zb,za)
      amp(2,1,2,1)=virt_gmgmjt_nfgammaMHV(i1,i2,i3,i5,i4,zb,za)
      amp(2,2,1,1)=virt_gmgmjt_nfgammaMHV(i1,i2,i4,i5,i3,zb,za)
!====== conjugation + line reversal 
      amp(1,1,1,2)=virt_gmgmjt_nfgammaMHV(i2,i1,i3,i4,i5,zb,za)
      amp(1,1,2,1)=virt_gmgmjt_nfgammaMHV(i2,i1,i3,i5,i4,zb,za)
      amp(1,2,1,1)=virt_gmgmjt_nfgammaMHV(i2,i1,i4,i5,i3,zb,za)

      return 
      end 

      subroutine amp_virt_slc_gmgmjt(i1,i2,i3,i4,i5,za,zb,amp)
      implicit none 
      include 'constants.f' 
      include 'zprods_decl.f' 
      integer i1,i2,i3,i4,i5
      double complex amp(2,2,2,2) 
      double complex amp_2gam1g

      call amp_virt_3gam(i1,i2,i3,i4,i5,za,zb,amp)

      return 
      end 

