

!------------ cutting routine for DM 
      subroutine dm_cuts(p,failed) 
      implicit none 
      include 'constants.f' 
!      include 'dm_params.f' 
      include 'jetlabel.f' 
      include 'leptcuts.f' 
      double precision p(mxpart,4)
      logical failed,first 
      double precision met_min,met_calc
      data first /.true./
      save first
      double precision dphi_jj 
      double precision phi_jjmax

!======= default is to pass 

      failed=.false. 
      
      met_min=misspt
      phi_jjmax=2.5d0
!      pt_jmin=100d0 

    
!====== missing ET cuts
!====== met p_T(3,4) always 

      met_calc=(p(3,1)+p(4,1))**2+(p(3,2)+p(4,2))**2 
      met_calc=dsqrt(met_calc) 

      if(met_calc.lt.met_min) failed=.true. 

      if(first) then 
         first=.false.
         write(6,*) '************ DM CUTS ****************'
         write(6,55) 'cutting MET at ',met_min 
!         write(6,55) 'cutting jet pt at ',pt_jmin 
         write(6,55) 'Delta Phi <  ',phi_jjmax
         write(6,*) '*************************************'
      endif

      if(jets.eq.2) then 
!--- cut on phi_jj 
        dphi_jj=
     .   (p(5,1)*p(6,1)+p(5,2)*p(6,2))
     .   /dsqrt((p(5,1)**2+p(5,2)**2)
     .         *(p(6,1)**2+p(6,2)**2))
      if (dphi_jj .lt. -0.999999999D0) dphi_jj=-1d0
      dphi_jj=dacos(dphi_jj) 
      if(dphi_jj.gt.phi_jjmax) failed=.true. 
      endif
         
 55   format(1x,a20,f8.3)

      
      return 
      end 
      
