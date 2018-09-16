      function higgsprop(s)
      implicit none
      include 'types.f'
      complex(dp):: higgsprop
c--- computes Higgs propagator for Higgs boson four-momentum squared s
c--- if CPscheme = .true. then it is computed in the complex pole
c---   scheme (Goria, Passarino, Rosco, arXiv:1112.5517,
c---           and c.f. Eq. (2.11) of arXiv:1206.4803)
c--- otherwise it takes the usual Breit-Wigner form
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'cpscheme.f'
      include 'first.f'
      real(dp):: s,mhbarsq,mhbar,gammahbar
      complex(dp):: cfac
      save mhbarsq,cfac

      if (CPscheme) then
c--- complex pole scheme propagator      
        if (first) then
          mhbarsq=hmass**2+hwidth**2
          mhbar=sqrt(mhbarsq)
          gammahbar=mhbar/hmass*hwidth
          cfac=cplx2(one,gammahbar/mhbar)
          first=.false.
        write(6,*)
        write(6,*)'****************************************************'
        write(6,*)'*  Using complex pole scheme for Higgs propagator  *'
        write(6,99) mhbar,gammahbar
        write(6,*)'****************************************************'
        write(6,*)
        endif
        higgsprop=cfac/(s*cfac-cplx2(mhbarsq,zip))
      else
c--- Breit Wigner propagator      
        higgsprop=one/cplx2(s-hmass**2,hmass*hwidth)
      endif
      
      return
   
   99 format(' *    MHB = ',f9.4,' GeV    GHB = ',f9.4,' GeV    *')
      
      end
      
      
      
