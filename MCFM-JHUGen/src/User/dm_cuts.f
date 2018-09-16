

!------------ cutting routine for DM
      subroutine dm_cuts(p,failed)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
!      include 'dm_params.f'
      include 'jetlabel.f'
      include 'leptcuts.f'
      real(dp):: p(mxpart,4)
      logical:: failed,first
      real(dp):: met_min,met_calc
      data first /.true./
      save first
      real(dp):: dphi_jj
      real(dp):: phi_jjmax

!======= default is to pass

      failed=.false.

      met_min=misspt
      phi_jjmax=2.5_dp
!      pt_jmin=100._dp


!====== missing ET cuts
!====== met p_T(3,4) always

      met_calc=(p(3,1)+p(4,1))**2+(p(3,2)+p(4,2))**2
      met_calc=sqrt(met_calc)

      if(met_calc<met_min) failed=.true.

      if(first) then
         first=.false.
         write(6,*) '************ DM CUTS ****************'
         write(6,55) 'cutting MET at ',met_min
!         write(6,55) 'cutting jet pt at ',pt_jmin
         write(6,55) 'Delta Phi <  ',phi_jjmax
         write(6,*) '*************************************'
      endif

      if(jets==2) then
!--- cut on phi_jj
        dphi_jj=
     &   (p(5,1)*p(6,1)+p(5,2)*p(6,2))
     &   /sqrt((p(5,1)**2+p(5,2)**2)
     &         *(p(6,1)**2+p(6,2)**2))
      if (dphi_jj < -0.999999999_dp) dphi_jj=-1._dp
      dphi_jj=acos(dphi_jj)
      if(dphi_jj>phi_jjmax) failed=.true.
      endif

 55   format(1x,a20,f8.3)


      return
      end

