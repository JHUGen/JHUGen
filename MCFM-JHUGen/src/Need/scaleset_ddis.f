      subroutine scaleset_ddis(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c--- "double deep-inelastic scale" defined for the single top t-channel by:
c---        Q^2=-(pt-pb)^2 on the light quark line
c---        Q^2+mt^2 on the heavy quark line
c---  (c.f. Z. Sullivan, hep-ph/0408049)
c--- AVAILABLE AT LEADING ORDER ONLY
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      include 'masses.f'
      include 'kpart.f'
      real(dp):: p(mxpart,4),mu0
      real(dp):: b1scale,q2scale,q1scale,b2scale
      common/bqscale/b1scale,q2scale,q1scale,b2scale
!$omp threadprivate(/bqscale/)
      if((kcase==kbq_tpq) .or.
     &   (kcase==kqg_tbq)) then
        if (kcase==kbq_tpq) then
          q1scale=-(
     &         +(p(3,4)+p(4,4)+p(5,4)+p(2,4))**2
     &         -(p(3,1)+p(4,1)+p(5,1)+p(2,1))**2
     &         -(p(3,2)+p(4,2)+p(5,2)+p(2,2))**2
     &         -(p(3,3)+p(4,3)+p(5,3)+p(2,3))**2)
          b2scale=-(
     &         +(p(3,4)+p(4,4)+p(5,4)+p(2,4))**2
     &         -(p(3,1)+p(4,1)+p(5,1)+p(2,1))**2
     &         -(p(3,2)+p(4,2)+p(5,2)+p(2,2))**2
     &         -(p(3,3)+p(4,3)+p(5,3)+p(2,3))**2)+mt**2
          q2scale=-(
     &         +(p(3,4)+p(4,4)+p(5,4)+p(1,4))**2
     &         -(p(3,1)+p(4,1)+p(5,1)+p(1,1))**2
     &         -(p(3,2)+p(4,2)+p(5,2)+p(1,2))**2
     &         -(p(3,3)+p(4,3)+p(5,3)+p(1,3))**2)
          b1scale=-(
     &         +(p(3,4)+p(4,4)+p(5,4)+p(1,4))**2
     &         -(p(3,1)+p(4,1)+p(5,1)+p(1,1))**2
     &         -(p(3,2)+p(4,2)+p(5,2)+p(1,2))**2
     &         -(p(3,3)+p(4,3)+p(5,3)+p(1,3))**2)+mt**2
        else
          q1scale=-(
     &         +(p(5,4)+p(1,4))**2
     &         -(p(5,1)+p(1,1))**2
     &         -(p(5,2)+p(1,2))**2
     &         -(p(5,3)+p(1,3))**2)
          b2scale=-(
     &         +(p(5,4)+p(1,4))**2
     &         -(p(5,1)+p(1,1))**2
     &         -(p(5,2)+p(1,2))**2
     &         -(p(5,3)+p(1,3))**2)+mt**2
          q2scale=-(
     &         +(p(5,4)+p(2,4))**2
     &         -(p(5,1)+p(2,1))**2
     &         -(p(5,2)+p(2,2))**2
     &         -(p(5,3)+p(2,3))**2)
          b1scale=-(
     &         +(p(5,4)+p(2,4))**2
     &         -(p(5,1)+p(2,1))**2
     &         -(p(5,2)+p(2,2))**2
     &         -(p(5,3)+p(2,3))**2)+mt**2
      endif
      q1scale=max(sqrt(q1scale),one) ! min. of 1 GeV for safety
      b2scale=max(sqrt(b2scale),one)
      q2scale=max(sqrt(q2scale),one)
      b1scale=max(sqrt(b1scale),one)
      mu0=q1scale   ! for safety
      else
        write(6,*) 'dynamicscale DDIS not supported for this process.'
        stop
      endif
      
      if (kpart.ne.klord) then
        write(6,*) 'dynamicscale DDIS only available at LO.'
        stop      
      endif
      
      return
      end
      
