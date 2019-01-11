      subroutine scaleset_m34sqsumptjsq(p,mu0)
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(m(34)^2+sum(ptj^2)), where m(34) is the invariant mass of (3+4)
c--- and the sum is over the pt-squared of each jet in the event
      implicit none
      include 'constants.f'
      include 'process.f'
      include 'jetlabel.f'
      include 'npart.f'
      include 'breit.f'
      double precision dot
      integer isub,oldjets
      double precision p(mxpart,4),pjet(mxpart,4),mu0,pt,rcut,sumptsq
      common/rcut/rcut

      if((case .eq. 'W_1jet') .or.
     &   (case .eq. 'Z_1jet') .or.
     &   (case .eq. 'ggfus1') .or.
     &   (case .eq. 'W_2jet') .or.
     &   (case .eq. 'Z_2jet') .or.
     &   (case .eq. 'ggfus2') .or.
     &   (case .eq. 'gmgmjt') .or.
     &   (case .eq. 'gmgmjj')) then

c--- first work out whether this point is real radiation or not
        if (abs(p(npart+2,4)) .gt. 1d-8) then
          isub=0  ! real term
        else
          isub=1  ! subtraction term
        endif
      
c-- cluster jets but make sure recorded number of jets is not changed
        oldjets=jets

        if ((case .eq. 'gmgmjt') .or.(case .eq. 'gmgmjj')) then
c--- if a photon process, call specific jet clustering routine
          call genclustphotons(p,rcut,pjet,isub)
        else
c--- otherwise, generic routine
          call genclust2(p,rcut,pjet,isub)
        endif
      
c        write(6,*) 'partons:'
c        if (p(5,4) .ge. 1d-8) write(6,*) 'pt5',pt(5,p)
c        if (p(6,4) .ge. 1d-8) write(6,*) 'pt6',pt(6,p)  
c      write(6,*) 'jets:'
c        if (pjet(5,4) .ge. 1d-8) write(6,*) 'pt5',pt(5,pjet)
c        if (pjet(6,4) .ge. 1d-8) write(6,*) 'pt6',pt(6,pjet)
c        write(6,*)

c--- work out sum of pt-squared
        if     (jets .eq. 1) then
          sumptsq=pt(5,pjet)**2
        elseif (jets .eq. 2) then
          sumptsq=pt(5,pjet)**2+pt(6,pjet)**2
        elseif (jets .eq. 3) then
          sumptsq=pt(5,pjet)**2+pt(6,pjet)**2+pt(7,pjet)**2
        else
          sumptsq=0d0
        endif
      
c--- restore old value of jets
        jets=oldjets

        mu0=2d0*dot(p,3,4)+sumptsq
        mu0=dsqrt(dabs(mu0))
      
      else
        write(6,*) 'dynamicscale sqrt(m(34)^2+sumptj^2)'//
     &             ' not supported for this process.'
        stop
      endif
      
      return
      end
      
