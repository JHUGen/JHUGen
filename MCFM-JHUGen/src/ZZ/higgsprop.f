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
      !include 'first.f'
      double precision s,mhbarsq,mhbar,gammahbar
      double complex cfac
      save mhbarsq,cfac

      if (CPscheme) then
c--- complex pole scheme propagator      
        mhbarsq=hmass**2+hwidth**2
        mhbar=sqrt(mhbarsq)
        gammahbar=mhbar/hmass*hwidth
        cfac=dcmplx(1d0,gammahbar/mhbar)
        higgsprop=cfac/(s*cfac-dcmplx(mhbarsq,0d0))
      else
c--- Breit Wigner propagator      
        higgsprop=1d0/dcmplx(s-hmass**2,hmass*hwidth)
      endif
      
      return
   
   99 format(' *    MHB = ',f9.4,' GeV    GHB = ',f9.4,' GeV    *')
      
      end
      

      double complex function higgs2prop(s)
c--- computes HM Higgs propagator for Higgs boson four-momentum squared s
c--- if CPscheme = .true. then it is computed in the complex pole
c---   scheme (Goria, Passarino, Rosco, arXiv:1112.5517,
c---           and c.f. Eq. (2.11) of arXiv:1206.4803)
c--- otherwise it takes the usual Breit-Wigner form
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'cpscheme.f'
      !include 'first.f'
      include 'spinzerohiggs_anomcoupl.f'
      double precision s,mhbarsq,mhbar,gammahbar
      double complex cfac
      save mhbarsq,cfac

      if(h2mass.lt.zip) then
         higgs2prop=czip
         return
      endif

      if (CPscheme) then
c--- complex pole scheme propagator      
        mhbarsq=h2mass**2+h2width**2
        mhbar=sqrt(mhbarsq)
        gammahbar=mhbar/h2mass*h2width
        cfac=dcmplx(1d0,gammahbar/mhbar)
        higgs2prop=cfac/(s*cfac-dcmplx(mhbarsq,0d0))
      else
c--- Breit Wigner propagator      
        higgs2prop=1d0/dcmplx(s-h2mass**2,h2mass*h2width)
      endif
      
      return
   
   99 format(' *    MHB = ',f9.4,' GeV    GHB = ',f9.4,' GeV    *')
      
      end
      
      
