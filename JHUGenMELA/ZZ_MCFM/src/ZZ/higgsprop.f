      double complex function higgsprop(s)
c--- computes Higgs propagator for Higgs boson four-momentum squared s
c--- if CPscheme = .true. then it is computed in the complex pole
c---   scheme (Goria, Passarino, Rosco, arXiv:1112.5517,
c---           and c.f. Eq. (2.11) of arXiv:1206.4803)
c--- otherwise it takes the usual Breit-Wigner form
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'cpscheme.f'
      include 'first.f'
      double precision s,mhbarsq,mhbar,gammahbar
      double complex cfac
      save mhbarsq,cfac

      if (CPscheme) then
c--- complex pole scheme propagator      
        if (first) then
          mhbarsq=hmass**2+hwidth**2
          mhbar=sqrt(mhbarsq)
          gammahbar=mhbar/hmass*hwidth
          cfac=dcmplx(1d0,gammahbar/mhbar)
          first=.false.
        write(6,*)
        write(6,*)'****************************************************'
        write(6,*)'*  Using complex pole scheme for Higgs propagator  *'
        write(6,99) mhbar,gammahbar
        write(6,*)'****************************************************'
        write(6,*)
        endif
        higgsprop=cfac/(s*cfac-dcmplx(mhbarsq,0d0))
      else
c--- Breit Wigner propagator      
        higgsprop=1d0/dcmplx(s-hmass**2,hmass*hwidth)
      endif
      
      return
   
   99 format(' *    MHB = ',f9.4,' GeV    GHB = ',f9.4,' GeV    *')
      
      end
      
      
      
      
      
      
      double complex function SMhiggsprop(s)
c--- computes Higgs propagator for Higgs boson four-momentum squared s
c--- if CPscheme = .true. then it is computed in the complex pole
c---   scheme (Goria, Passarino, Rosco, arXiv:1112.5517,
c---           and c.f. Eq. (2.11) of arXiv:1206.4803)
c--- otherwise it takes the usual Breit-Wigner form
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'cpscheme.f'
      include 'first.f'
      double precision s,mhbarsq,mhbar,gammahbar
      double complex cfac
      double precision, parameter :: SMhmass=125.0d0, SMhwidth=0.00407d0
      save mhbarsq,cfac
      
      
!SMhmass, SMhwidth
!       if (CPscheme) then
! c--- complex pole scheme propagator      
!         if (first) then
!           mhbarsq=SMhmass**2+SMhwidth**2
!           mhbar=sqrt(mhbarsq)
!           gammahbar=mhbar/SMhmass*hwidth
!           cfac=dcmplx(1d0,gammahbar/mhbar)
!           first=.false.
!         write(6,*)
!         write(6,*)'****************************************************'
!         write(6,*)'* Using complex pole scheme for SM Higgs propagator*'
!         write(6,99) mhbar,gammahbar
!         write(6,*)'****************************************************'
!         write(6,*)
!         endif
!         SMhiggsprop=cfac/(s*cfac-dcmplx(mhbarsq,0d0))
!       else
! c--- Breit Wigner propagator      
        SMhiggsprop=1d0/dcmplx(s-SMhmass**2,SMhmass*SMhwidth)
!       endif
      
      return
   
   99 format(' *    MHB = ',f9.4,' GeV    GHB = ',f9.4,' GeV    *')
      
      end
      