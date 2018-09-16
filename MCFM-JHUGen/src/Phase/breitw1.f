      subroutine breitw1(msq,mminsq,mmaxsq,rmass,rwidth,wt)       
      implicit none
c---- Given a mass-squared msq generate a weight wt 
c---- breit-wigner should still be included in the matrix element
c      double precision bw
      include 'constants.f'
      include 'zerowidth.f'
      double precision mminsq,mmaxsq,rmass,rwidth,msq,wt
      double precision almin,almax,tanal

      if (zerowidth) then
          tanal=0d0
          almax=+pi/two
          almin=-pi/two
      else
          almin=atan((mminsq-rmass**2)/rmass/rwidth)
          almax=atan((mmaxsq-rmass**2)/rmass/rwidth)
          tanal=(rmass**2-msq)/(rmass*rwidth)
      endif

      wt=(almax-almin)*rmass*rwidth*(1d0+tanal**2)
      return
      end

