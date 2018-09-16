!===== C. Williams Feb 2015
!----- routine for calculating gg=>ZH
!----- note that this routine contains two types of pieces
!------
!_---- 1) Amp box => Feynman diagrams which couple the H directly to the top
!----- 2) Amp tri => Feynman diagrams which couple the H directly to the Z

!==== the routines are called in sepearate files using QCDLoop, this is just
!==== the assemble routine
      subroutine gg_zh(p,msqgg)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'hdecaymode.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'first.f'
      include 'cutoff.f'
      real(kind=dp):: p(mxpart,4),msqgg
      complex(kind=dp):: ampbox_mt(2,2,2),amptri_mt(2,2,2)
      complex(kind=dp):: ampbox_mb(2,2,2),amptri_mb(2,2,2),amptot(2,2,2)
      real(kind=dp):: fac,amp_sq,hdecay,pttwo
      real(kind=dp):: fac_Htop,fac_HZ,fac_Hbot,fac_HZ_top,fac_HZ_bot
      integer:: h1,h2,h3,j,k
      real(kind=dp):: s56,msqgamgam
      real(kind=dp):: v1(2),ggHZ_eft,cutoff_orig
      integer:: nproc
      common/nproc/nproc

      msqgg=zip

c--- numerical safety cut
      if (pttwo(3,4,p) < 0.1_dp) return

      cutoff_orig=cutoff

!==== if Z=>bb this routine is not correct so exit gracefully
      if((nproc==108).or.(nproc==103)) then
         if(first) then
            first=.false.
            write(6,*) 'Z=>bb not implemented for gg=>HZ check coupling'
            write(6,*) 'and decay in gg_zh.f'
            write(6,*) 'Setting msq(g,g)=0'
         endif
         return
      endif

      if(first) then
         call qlinit
         first=.false.
      endif


      ampbox_mt(:,:,:)=czip
      amptri_mt(:,:,:)=czip
      ampbox_mb(:,:,:)=czip
      amptri_mb(:,:,:)=czip
      amptot(:,:,:)=czip

      amp_sq=zip
      call spinoru(4,p,za,zb)

!=====this routine runs in NNLO mode, and can have very small
!=====values of s, since its non-singular lets ensure stability
!==== by doing a smalls cut.
      cutoff=1.E-3_dp
      call smalls(s,4,*999)
      cutoff=cutoff_orig

!------ top loops
      call fill_basis_int_ggzh(1,2,3,4,mt**2)
      call gg_HZ_box(p,ampbox_mt,mt**2)
      call gg_HZ_tri(p,amptri_mt,mt**2)

!------ bottom loops
      call fill_basis_int_ggzh(1,2,3,4,mb**2)
      call gg_HZ_box(p,ampbox_mb,mb**2)
      call gg_HZ_tri(p,amptri_mb,mb**2)

!====== Higgs decay
      if (hdecaymode == 'tlta') then
         s56=s(5,6)+two*mtau**2
         call htautaudecay(p,5,6,hdecay)
      elseif (hdecaymode == 'bqba') then
         s56=s(5,6)+2*mb**2
         call hbbdecay(p,5,6,hdecay)
      elseif (hdecaymode == 'gaga') then
         s56=s(5,6)
         hdecay=msqgamgam(hmass)
      elseif (hdecaymode == 'wpwm') then
         s56=s(5,6)+s(5,7)+s(5,8)+s(6,7)+s(6,8)+s(7,8)
         call hwwdecay(p,5,6,7,8,hdecay)
      else
      write(6,*) 'Unimplemented process in gg_hgg_v'
      stop
      endif
      hdecay=hdecay/((s56-hmass**2)**2+(hmass*hwidth)**2)


!---- couplings for EW side of things
!      write(6,*) 'left handed top coupling ',l(2)
!      write(6,*) 'right handed top coupling ',r(2)
!====== this is couplings for leptonic decay
      v1(1)=l1
      v1(2)=r1
!----- prefactor extracted from KC
! =i_/2/gw/gs^2/ee/Mw*costhw^2* for Z and
!==== gw = ee/sinthw
      fac_HZ=esq/(one-xw)/sqrt(xw)*wmass*gsq
!===== the amplitude we calculated was for the left handed coupling
!===== of Z to top, also need R
      fac_HZ_top=fac_HZ*(l(2)-r(2))
      fac_HZ_bot=fac_HZ*(l(5)-r(5))
!----- prefactor extracted from KC
!     =i_/2/gw/gs^2/ee*Mw/mt^2*
!==== gw=ee/sinthw
      fac_Htop=esq/sqrt(xw)*mt**2/wmass*gsq
!===== the amplitude we calculated was for the left handed coupling
!===== of Z to top, also need R
      fac_Htop=fac_Htop*(l(2)-r(2))
      fac_Hbot=esq/sqrt(xw)*mb**2/wmass*gsq
!===== the amplitude we calculated was for the left handed coupling
!===== of Z to top, also need R
      fac_Hbot=fac_Hbot*(l(5)-r(5))



      do h1=1,2
         do h2=1,2
            do h3=1,2
   !===== dress with  Higgs couplings here
            amptot(h1,h2,h3)=ampbox_mt(h1,h2,h3)*fac_Htop*v1(h3)
     &              +amptri_mt(h1,h2,h3)*fac_HZ_top*v1(h3)
     &              +ampbox_mb(h1,h2,h3)*fac_Hbot*v1(h3)
     &              +amptri_mb(h1,h2,h3)*fac_HZ_bot*v1(h3)
            enddo
         enddo
      enddo


      amp_sq=zip
      do h1=1,2
         do h2=1,2
            do h3=1,2
               amp_sq=amp_sq+abs(amptot(h1,h2,h3))**2
            enddo
         enddo
      enddo

!----- color and initial state and Z=>ee decay and loop factor (squared)
      fac=avegg*V*esq/(16*pisq)**2
      msqgg=fac*hdecay*amp_sq

      return

 999  continue
      cutoff=cutoff_orig
      return

      end



      subroutine fill_basis_int_ggzh(i1,i2,i3,i4,mt2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'scale.f'
      integer:: i1,i2,i3,i4
      real(kind=dp):: p(mxpart,4),mt2
      real(kind=dp):: s12,mZsq,s25,s15,t,mHsq
      integer:: Nbox,Ntri
      parameter(Nbox=3,Ntri=6)
      complex(kind=dp):: D0(Nbox),C0(Ntri)
      integer:: d25_12,d15_12,d15_25
      parameter(d25_12=1,d15_12=2,d15_25=3)
      integer:: c25_Z,cH_25,c12,c15_H,cZ_15,c12_Z_H
      parameter(c25_Z=1,cH_25=2,c12=3,c15_H=4,cZ_15=5,c12_Z_H=6)
      complex(kind=dp):: qlI4,qlI3
      common/ggZH_basisint/D0,C0
!$omp threadprivate(/ggZH_basisint/)


      t(i1,i2,i3)=s(i1,i2)+s(i2,i3)+s(i1,i3)
      s12=s(i1,i2)
      mZsq=s(i3,i4)
      s25=t(i1,i3,i4)
      s15=t(i2,i3,i4)
      mHsq=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

!----- note that we only need one sort of routine, since the
!----- worst that can happen is that 3<->4 to swap leptons, but
!----- this doesnt change any basis integral

      D0(d25_12)=qlI4(mZsq,zip,zip,mHsq,s25,s12,mt2,mt2,mt2,mt2,musq,0)
      D0(d15_12)=qlI4(mZsq,zip,zip,mHsq,s15,s12,mt2,mt2,mt2,mt2,musq,0)
      D0(d15_25)=qlI4(mHsq,zip,mZsq,zip,s15,s25,mt2,mt2,mt2,mt2,musq,0)

!----- fill integrals from QCD loop : triangles
      C0(c25_Z)=qlI3(s25,zip,mZsq,mt2,mt2,mt2,musq,0)
      C0(cH_25)=qlI3(mHsq,zip,s25,mt2,mt2,mt2,musq,0)
      C0(c12)=qlI3(s12,zip,zip,mt2,mt2,mt2,musq,0)
      C0(c15_H)=qlI3(s15,zip,mHsq,mt2,mt2,mt2,musq,0)
      C0(cZ_15)=qlI3(mZsq,zip,s15,mt2,mt2,mt2,musq,0)
      C0(c12_Z_H)=qlI3(s12,mZsq,mHsq,mt2,mt2,mt2,musq,0)

      return
      end

