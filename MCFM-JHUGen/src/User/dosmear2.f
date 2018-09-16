      function dosmear2(Et,etarap,xran)
      implicit none
      include 'types.f'
      real(dp):: dosmear2
c--- Improved jet pT smearing routine, based on jet rapidity
c--- given values of Et,rapidity for a jet and a random number distributed
c--- according to a unit Gaussian, returns a smeared value for the Et
c--- corresponding to a standard deviation of 
c---  sigma=lambda*(-0.0317+0.578*Sqrt(Et)+0.047*Et)
c--- where lambda depends on the jet rapidity
      
      real(dp):: Et,etarap,xran,lambda,sigma
      
      if     ((etarap >= -2.4_dp) .and. (etarap < -2.0_dp)) then
        lambda=1.56_dp
      elseif ((etarap >= -2.0_dp) .and. (etarap < -1.5_dp)) then
        lambda=1.08_dp-0.0022*Et
      elseif ((etarap >= -1.5_dp) .and. (etarap < -0.9_dp)) then
        lambda=1.22_dp+0.0014*Et
      elseif ((etarap >= -0.9_dp) .and. (etarap < -0.2_dp)) then
        lambda=1.00_dp
      elseif ((etarap >= -0.2_dp) .and. (etarap <  0.2_dp)) then
        lambda=1.10_dp+0.0023*Et
      elseif ((etarap >=  0.2_dp) .and. (etarap <  0.9_dp)) then
        lambda=1.00_dp
      elseif ((etarap >=  0.9_dp) .and. (etarap <  1.5_dp)) then
        lambda=1.29_dp+0.0013*Et
      elseif ((etarap >=  1.5_dp) .and. (etarap <  2.0_dp)) then
        lambda=1.27_dp-0.0050*Et
      elseif ((etarap >=  2.0_dp) .and. (etarap <  2.4_dp)) then
        lambda=1.46_dp
      else
c--- average value beyond |eta|=2.4, shouldn't matter anyway
        lambda=1.5_dp
      endif
      
      sigma=lambda*(-0.0317_dp+0.578_dp*sqrt(Et)+0.047_dp*Et)
      if (sigma < 0._dp) then
c        write(6,*) 'WARNING: smearing standard deviation < 0'
        sigma=0._dp
      endif
      
      dosmear2=Et+sigma*xran
            
      return
      end
