      subroutine scaleset_ddis(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c--- "double deep-inelastic scale" defined for the single top t-channel by:
c---        Q^2=-(pt-pb)^2 on the light quark line
c---        Q^2+mt^2 on the heavy quark line
c---  (c.f. Z. Sullivan, hep-ph/0408049)
c--- AVAILABLE AT LEADING ORDER ONLY
      implicit none
      include 'constants.f'
      include 'process.f'
      include 'masses.f'
      include 'part.f'
      double precision p(mxpart,4),mu0
      double precision b1scale,q2scale,q1scale,b2scale
      common/bqscale/b1scale,q2scale,q1scale,b2scale
!$omp threadprivate(/bqscale/)
      if((case .eq. 'bq_tpq') .or.
     &   (case .eq. 'qg_tbq')) then
        if (case .eq. 'bq_tpq') then
          q1scale=-(
     .         +(p(3,4)+p(4,4)+p(5,4)+p(2,4))**2
     .         -(p(3,1)+p(4,1)+p(5,1)+p(2,1))**2
     .         -(p(3,2)+p(4,2)+p(5,2)+p(2,2))**2
     .         -(p(3,3)+p(4,3)+p(5,3)+p(2,3))**2)
          b2scale=-(
     .         +(p(3,4)+p(4,4)+p(5,4)+p(2,4))**2
     .         -(p(3,1)+p(4,1)+p(5,1)+p(2,1))**2
     .         -(p(3,2)+p(4,2)+p(5,2)+p(2,2))**2
     .         -(p(3,3)+p(4,3)+p(5,3)+p(2,3))**2)+mt**2
          q2scale=-(
     .         +(p(3,4)+p(4,4)+p(5,4)+p(1,4))**2
     .         -(p(3,1)+p(4,1)+p(5,1)+p(1,1))**2
     .         -(p(3,2)+p(4,2)+p(5,2)+p(1,2))**2
     .         -(p(3,3)+p(4,3)+p(5,3)+p(1,3))**2)
          b1scale=-(
     .         +(p(3,4)+p(4,4)+p(5,4)+p(1,4))**2
     .         -(p(3,1)+p(4,1)+p(5,1)+p(1,1))**2
     .         -(p(3,2)+p(4,2)+p(5,2)+p(1,2))**2
     .         -(p(3,3)+p(4,3)+p(5,3)+p(1,3))**2)+mt**2
        else
          q1scale=-(
     .         +(p(5,4)+p(1,4))**2
     .         -(p(5,1)+p(1,1))**2
     .         -(p(5,2)+p(1,2))**2
     .         -(p(5,3)+p(1,3))**2)
          b2scale=-(
     .         +(p(5,4)+p(1,4))**2
     .         -(p(5,1)+p(1,1))**2
     .         -(p(5,2)+p(1,2))**2
     .         -(p(5,3)+p(1,3))**2)+mt**2
          q2scale=-(
     .         +(p(5,4)+p(2,4))**2
     .         -(p(5,1)+p(2,1))**2
     .         -(p(5,2)+p(2,2))**2
     .         -(p(5,3)+p(2,3))**2)
          b1scale=-(
     .         +(p(5,4)+p(2,4))**2
     .         -(p(5,1)+p(2,1))**2
     .         -(p(5,2)+p(2,2))**2
     .         -(p(5,3)+p(2,3))**2)+mt**2
      endif
      q1scale=max(dsqrt(q1scale),1d0) ! min. of 1 GeV for safety
      b2scale=max(dsqrt(b2scale),1d0)
      q2scale=max(dsqrt(q2scale),1d0)
      b1scale=max(dsqrt(b1scale),1d0)
      mu0=q1scale   ! for safety
      else
        write(6,*) 'dynamicscale DDIS not supported for this process.'
        stop
      endif
      
      if (part .ne. 'lord') then
        write(6,*) 'dynamicscale DDIS only available at LO.'
        stop      
      endif
      
      return
      end
      
