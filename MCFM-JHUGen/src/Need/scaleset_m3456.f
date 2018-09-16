      subroutine scaleset_m3456(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c--- invariant mass of particles 3, 4, 5 and 6
      implicit none
      include 'constants.f'
      include 'process.f'
      double precision p(mxpart,4),mu0

      if(
     &   (case .eq. 'WWqqbr') .or.
     &   (case .eq. 'WWnpol') .or.
     &   (case .eq. 'WW_jet') .or.
     &   (case .eq. 'WZbbar') .or.
     &   (case .eq. 'ZZlept') .or.
     &   (case .eq. 'WHbbar') .or.
     &   (case .eq. 'WHgaga') .or.
     &   (case .eq. 'ZHbbar') .or.
     &   (case .eq. 'ZHgaga') .or.
     &   (case .eq. 'HWW_4l') .or.
     &   (case .eq. 'HWW_tb') .or.
     &   (case .eq. 'HWWint') .or.
     &   (case .eq. 'HWWH+i') .or.
     &   (case .eq. 'ggWW4l') .or.
     &   (case .eq. 'ggWWbx') .or.
     &   (case .eq. 'HZZ_4l') .or.
     &   (case .eq. 'HZZ_tb') .or.
     &   (case .eq. 'HZZint') .or.
     &   (case .eq. 'HZZH+i') .or.
     &   (case .eq. 'HZZqgI') .or.
     &   (case .eq. 'HZZ+jt') .or.
     &   (case .eq. 'HWWjet') .or.
     &   (case .eq. 'HZZjet') .or.
     &   (case .eq. 'HWW2jt') .or.
     &   (case .eq. 'HZZ2jt') .or.
     &   (case .eq. 'qq_HZZ') .or.
     &   (case .eq. 'fourga') .or.
     &   (case .eq. 'W_2gam') .or.
     &   (case .eq. 'Z_2gam') .or.
     &   (case .eq. 'HmZZ4l') .or.
     &   (case .eq. 'HVV_tb') .or.
     &   (case .eq. 'ggVV4l') .or.
     &   (case .eq. 'ggVVbx') .or.
     &   (case .eq. 'ggZZ4l') .or.
     &   (case .eq. 'ggZZbx') .or.
     &   (case .eq. 'HmZZ4l') .or.
     &   (case .eq. 'qqZZqq') .or.
     &   (case .eq. 'qqWWqq') .or.
     &   (case .eq. 'qqVVqq')
     & ) then
        mu0=(p(3,4)+p(4,4)+p(5,4)+p(6,4))**2
     &     -(p(3,1)+p(4,1)+p(5,1)+p(6,1))**2
     &     -(p(3,2)+p(4,2)+p(5,2)+p(6,2))**2
     &     -(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2
        mu0=dsqrt(dabs(mu0))
      else
        write(6,*)'dynamicscale m(3456) not supported for this process.'
        stop
      endif

      return
      end

