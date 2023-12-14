      subroutine scaleset_randrange(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c--- invariant mass of particles 3, 4, 5 and 6
      implicit none
      include 'constants.f'
      include 'process.f'
      include 'facscale.f'
      include 'facscale_range.f'
      double precision p(mxpart,4),mu0
      if(facscale_high.lt.facscale_low) then
        write(6,*)'high scale bound must be < low!'
        write(6,*)'Please edit factorization_range.f'
        stop
      endif
      if(
     &   (case .eq. 'HZZ_4l') .or.
     &   (case .eq. 'HmZZ4l') .or.
     &   (case .eq. 'ggZZ4l') .or.
     &   (case .eq. 'ZZlept') .or.
     &   (case .eq. 'HZZ_tb')
     & ) then
        call random_number(mu0)
        mu0 = (facscale_high - facscale_low)*mu0 + facscale_low
      else
        write(6,*)'dynamicscale rand not supported for this process.'
        stop
      endif

      return
      end

