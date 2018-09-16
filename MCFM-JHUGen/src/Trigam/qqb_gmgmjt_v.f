!====== LO ROUTINE FOR GAM GAM JET using amplitudes

      subroutine qqb_gmgmjt_v(p,msq)
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
      include 'epinv.f'
      include 'nflav.f'
      include 'scheme.f'
      include 'scale.f'
      real(dp):: subuv
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: msq0(-nf:nf,-nf:nf)
      complex(dp):: qqbg_lo(2,2,2,2),qbqg_lo(2,2,2,2)
      complex(dp):: qgqb_lo(2,2,2,2),qbgq_lo(2,2,2,2)
      complex(dp):: gqqb_lo(2,2,2,2),gqbq_lo(2,2,2,2)
      complex(dp):: qqbg_v(2,2,2,2),qbqg_v(2,2,2,2)
      complex(dp):: qgqb_v(2,2,2,2),qbgq_v(2,2,2,2)
      complex(dp):: gqqb_v(2,2,2,2),gqbq_v(2,2,2,2)
      complex(dp):: qqbg_nf(2,2,2,2),qbqg_nf(2,2,2,2)
      complex(dp):: qgqb_nf(2,2,2,2),qbgq_nf(2,2,2,2)
      complex(dp):: gqqb_nf(2,2,2,2),gqbq_nf(2,2,2,2)
      real(dp):: qqbg_sum,qbqg_sum
      real(dp):: qgqb_sum,qbgq_sum
      real(dp):: gqbq_sum,gqqb_sum
      real(dp):: qqbg_nfs,qbqg_nfs
      real(dp):: qgqb_nfs,qbgq_nfs
      real(dp):: gqbq_nfs,gqqb_nfs
      integer:: h1,h2,h3,h4,j,k
      real(dp):: fac,statfac
      parameter(statfac=0.5_dp)
      complex(dp):: test,virt_gmgmjt_GaMHV
      real(dp):: faclo,nf_fac
      real(dp):: qsum

      qbqg_sum=0._dp
      qqbg_sum=0._dp
      qbgq_sum=0._dp
      gqbq_sum=0._dp
      qgqb_sum=0._dp
      gqqb_sum=0._dp

      qbqg_nfs=0._dp
      qqbg_nfs=0._dp
      qbgq_nfs=0._dp
      gqbq_nfs=0._dp
      qgqb_nfs=0._dp
      gqqb_nfs=0._dp
      qsum=0._dp
      msq(:,:)=0._dp
      scheme='dred'


      faclo=8._dp*cf*xn*gsq*esq**2*statfac
      fac=faclo*xn*ason2pi/2._dp
!----UV counterterm contains the finite renormalization to arrive
!----at MS bar scheme.
!====== for now call born separately. probably want to adjust this later
!====== to make code more efficient
      call qqb_gmgmjt(p,msq0)
      subuv=ason2pi*xn*(epinv*(11._dp-2._dp*real(nflav)/xn)-1._dp)/6._dp
      call spinoru(5,p,za,zb)
!===== DEBUG KC CHECK
!      include 'kinpoint_mcfm.f'
!      call spinorz(5,p,za,zb)
!      epinv=0_dp
!      scale=1_dp
!      musq=1_dp
!      call amp_virt_gmgmjt(1,2,3,4,5,za,zb,qqbg_v)
!      call amp_lord_gmgmjt(1,2,3,4,5,za,zb,qqbg_v)
!      write(6,*) 'lord = ',qqbg_v(1,2,2,1)
!      stop
!=====DEBUG END



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
     & real(conjg(qqbg_v(h1,h2,h3,h4))*qqbg_lo(h1,h2,h3,h4))
     & +real(conjg(qqbg_lo(h1,h2,h3,h4))*qqbg_v(h1,h2,h3,h4))
            qqbg_nfs=qqbg_nfs+
     & real(conjg(qqbg_nf(h1,h2,h3,h4))*qqbg_lo(h1,h2,h3,h4))
     & +real(conjg(qqbg_lo(h1,h2,h3,h4))*qqbg_nf(h1,h2,h3,h4))

            qbqg_sum=qbqg_sum+
     & real(conjg(qbqg_v(h1,h2,h3,h4))*qbqg_lo(h1,h2,h3,h4))
     & +real(conjg(qbqg_lo(h1,h2,h3,h4))*qbqg_v(h1,h2,h3,h4))
            qbqg_nfs=qbqg_nfs+
     & real(conjg(qbqg_nf(h1,h2,h3,h4))*qbqg_lo(h1,h2,h3,h4))
     & +real(conjg(qbqg_lo(h1,h2,h3,h4))*qbqg_nf(h1,h2,h3,h4))

            qgqb_sum=qgqb_sum+
     & real(conjg(qgqb_v(h1,h2,h3,h4))*qgqb_lo(h1,h2,h3,h4))
     & +real(conjg(qgqb_lo(h1,h2,h3,h4))*qgqb_v(h1,h2,h3,h4))
           qgqb_nfs=qgqb_nfs+
     & real(conjg(qgqb_nf(h1,h2,h3,h4))*qgqb_lo(h1,h2,h3,h4))
     & +real(conjg(qgqb_lo(h1,h2,h3,h4))*qgqb_nf(h1,h2,h3,h4))

            gqqb_sum=gqqb_sum+
     & real(conjg(gqqb_v(h1,h2,h3,h4))*gqqb_lo(h1,h2,h3,h4))
     & +real(conjg(gqqb_lo(h1,h2,h3,h4))*gqqb_v(h1,h2,h3,h4))
            gqqb_nfs=gqqb_nfs+
     & real(conjg(gqqb_nf(h1,h2,h3,h4))*gqqb_lo(h1,h2,h3,h4))
     & +real(conjg(gqqb_lo(h1,h2,h3,h4))*gqqb_nf(h1,h2,h3,h4))


            qbgq_sum=qbgq_sum+
     & real(conjg(qbgq_v(h1,h2,h3,h4))*qbgq_lo(h1,h2,h3,h4))
     & +real(conjg(qbgq_lo(h1,h2,h3,h4))*qbgq_v(h1,h2,h3,h4))
            qbgq_nfs=qbgq_nfs+
     & real(conjg(qbgq_nf(h1,h2,h3,h4))*qbgq_lo(h1,h2,h3,h4))
     & +real(conjg(qbgq_lo(h1,h2,h3,h4))*qbgq_nf(h1,h2,h3,h4))

            gqbq_sum=gqbq_sum+
     & real(conjg(gqbq_v(h1,h2,h3,h4))*gqbq_lo(h1,h2,h3,h4))
     & +real(conjg(gqbq_lo(h1,h2,h3,h4))*gqbq_v(h1,h2,h3,h4))
           gqbq_nfs=gqbq_nfs+
     & real(conjg(gqbq_nf(h1,h2,h3,h4))*gqbq_lo(h1,h2,h3,h4))
     & +real(conjg(gqbq_lo(h1,h2,h3,h4))*gqbq_nf(h1,h2,h3,h4))


         enddo
      enddo
      enddo
      enddo

      qsum=(2._dp*Q(2)**2+3._dp*Q(1)**2)
      do j=-nf,nf
         do k=-nf,nf
            msq(j,k)=0._dp
            if((j.ne.0).and.(k.ne.0).and.(j.ne.-k)) goto 20

            if((j<0).and.(k>0)) then
               msq(j,k)=fac*aveqq*qbqg_sum*Q(k)**4-subuv*msq0(j,k)
     &              +nf_fac*aveqq*qbqg_nfs*qsum*Q(k)**2
            elseif((j==0).and.(k>0)) then
               msq(j,k)=fac*aveqg*gqqb_sum*Q(k)**4-subuv*msq0(j,k)
     &              +nf_fac*aveqg*gqqb_nfs*qsum*Q(k)**2
            elseif((j>0).and.(k==0)) then
               msq(j,k)=fac*aveqg*qgqb_sum*Q(j)**4-subuv*msq0(j,k)
     &              +nf_fac*aveqg*qgqb_nfs*qsum*Q(j)**2
            elseif((j<0).and.(k==0)) then
               msq(j,k)=fac*aveqg*qbgq_sum*Q(j)**4-subuv*msq0(j,k)
     &              +nf_fac*aveqg*qbgq_nfs*qsum*Q(j)**2
            elseif((j==0).and.(k<0)) then
               msq(j,k)=fac*aveqg*gqbq_sum*Q(k)**4-subuv*msq0(j,k)
     &              +nf_fac*aveqg*gqbq_nfs*qsum*Q(k)**2
            elseif((j>0).and.(k<0)) then
               msq(j,k)=fac*aveqq*qqbg_sum*Q(j)**4-subuv*msq0(j,k)
     &              +nf_fac*aveqq*qqbg_nfs*qsum*Q(j)**2
            endif

 20         continue
         enddo
      enddo

      return
      end


      subroutine amp_virt_gmgmjt(i1,i2,i3,i4,i5,za,zb,amp)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5
      complex(dp):: amp(2,2,2,2),amp_lc(2,2,2,2)
      complex(dp):: amp_slc(2,2,2,2),amp_nf(2,2,2,2)
      integer:: h1,h2,h3,h4
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
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5
      complex(dp):: amp(2,2,2,2)
      complex(dp):: amp_2gam1g,virt_gmgmjt_GMHV,virt_gmgmjt_GaMHV
      complex(dp):: virt_gmgmjt_gammaMHV,virt_gmgmjt_gluonMHV
      integer h1,h2,h3,h4
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
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5
      complex(dp):: amp(2,2,2,2)
      complex(dp):: amp_2gam1g
      complex(dp):: virt_gmgmjt_nfGaMHV,virt_gmgmjt_nfgammaMHV
      complex(dp):: virt_gmgmjt_nfallp
!===== default is for gluon MHV amplitude

c--- DEBUG ONLY
c      write(6,*)
c      amp(1,2,2,1)=virt_gmgmjt_nfGaMHV(i1,i2,i3,i4,i5,za,zb)
c      amp(1,2,2,1)=virt_gmgmjt_nfgammaMHV(i1,i2,i3,i4,i5,za,zb)
c      write(6,*)
c      pause
c--- DEBUG ONLY

      amp(:,:,:,:)=czip
      amp(1,1,1,1)=virt_gmgmjt_nfallp(i1,i2,i3,i4,i5,zb,za)
      amp(2,2,2,2)=virt_gmgmjt_nfallp(i1,i2,i3,i4,i5,za,zb)
      amp(2,1,1,1)=virt_gmgmjt_nfallp(i2,i1,i3,i4,i5,zb,za)
      amp(1,2,2,2)=virt_gmgmjt_nfallp(i2,i1,i3,i4,i5,za,zb)
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
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5
      complex(dp):: amp(2,2,2,2)
      complex(dp):: amp_2gam1g
      integer h1,h2,h3,h4

      call amp_virt_3gam(i1,i2,i3,i4,i5,za,zb,amp)
      return
      end

