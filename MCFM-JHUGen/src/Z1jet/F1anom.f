      function F1anom(s12,s45,mt2,musq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'qlfirst.f'
      complex(dp)::F1anom
      complex(dp)::qlI3,qlI2,lnrat
      real(dp)::s12,s45,mt2,musq

      if (qlfirst) then
         call qlinit
         qlfirst=.false.
      endif

      if (mt2 .eq. zero) then
      F1anom=two/(s45-s12)*(one+s45/(s45-s12)*lnrat(-s12,-s45))
      return
      else
      F1anom=one/(s45-s12)
     & *(two+four*mt2*qlI3(s12,zero,s45,mt2,mt2,mt2,musq,0)
     & +(two+two*s12/(s45-s12))
     & *(qlI2(s45,mt2,mt2,musq,0)-qlI2(s12,mt2,mt2,musq,0)))
      endif
            
      return
      end

