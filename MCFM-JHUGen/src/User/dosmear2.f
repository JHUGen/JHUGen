      double precision function dosmear2(Et,etarap,xran)
c--- Improved jet pT smearing routine, based on jet rapidity
c--- given values of Et,rapidity for a jet and a random number distributed
c--- according to a unit Gaussian, returns a smeared value for the Et
c--- corresponding to a standard deviation of 
c---  sigma=lambda*(-0.0317+0.578*Sqrt(Et)+0.047*Et)
c--- where lambda depends on the jet rapidity
      implicit none
      double precision Et,etarap,xran,lambda,sigma
      
      if     ((etarap .ge. -2.4d0) .and. (etarap .lt. -2.0d0)) then
        lambda=1.56d0
      elseif ((etarap .ge. -2.0d0) .and. (etarap .lt. -1.5d0)) then
        lambda=1.08d0-0.0022*Et
      elseif ((etarap .ge. -1.5d0) .and. (etarap .lt. -0.9d0)) then
        lambda=1.22d0+0.0014*Et
      elseif ((etarap .ge. -0.9d0) .and. (etarap .lt. -0.2d0)) then
        lambda=1.00d0
      elseif ((etarap .ge. -0.2d0) .and. (etarap .lt.  0.2d0)) then
        lambda=1.10d0+0.0023*Et
      elseif ((etarap .ge.  0.2d0) .and. (etarap .lt.  0.9d0)) then
        lambda=1.00d0
      elseif ((etarap .ge.  0.9d0) .and. (etarap .lt.  1.5d0)) then
        lambda=1.29d0+0.0013*Et
      elseif ((etarap .ge.  1.5d0) .and. (etarap .lt.  2.0d0)) then
        lambda=1.27d0-0.0050*Et
      elseif ((etarap .ge.  2.0d0) .and. (etarap .lt.  2.4d0)) then
        lambda=1.46d0
      else
c--- average value beyond |eta|=2.4, shouldn't matter anyway
        lambda=1.5d0
      endif
      
      sigma=lambda*(-0.0317d0+0.578d0*dsqrt(Et)+0.047d0*Et)
      if (sigma .lt. 0d0) then
c        write(6,*) 'WARNING: smearing standard deviation < 0'
        sigma=0d0
      endif
      
      dosmear2=Et+sigma*xran
            
      return
      end
