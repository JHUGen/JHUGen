      subroutine setuptau(ipermille)
      implicit none
      include 'types.f'
c---- Routine to setup appropriate value of tau for performing SCET
c---- calculation at NNLO at "ipermille" per mille level of precision
      include 'taucut.f'
      include 'nproc.f'
      integer ipermille

      if ((ipermille == 10) .or. (ipermille == 2)) then
        continue
      else
        write(6,*) 'Invalid per mille accuracy in setuptau: ',ipermille
        stop
      endif
      
      if     ((nproc == 1) .or. (nproc == 6)) then
        if (ipermille == 10) taucut=5.e-3_dp
        if (ipermille ==  2) taucut=1.e-3_dp
      elseif ((nproc == 31) .or. (nproc == 32)) then
        if (ipermille == 10) taucut=1.e-2_dp
        if (ipermille ==  2) taucut=2.e-3_dp
      elseif ((nproc >= 91) .and. (nproc <= 100)) then
        if (ipermille == 10) taucut=2.e-1_dp
        if (ipermille ==  2) taucut=1.e-2_dp
      elseif ((nproc >= 101) .and. (nproc <= 110)) then
        if (ipermille == 10) taucut=3.e-1_dp
        if (ipermille ==  2) taucut=2.e-2_dp
      elseif ((nproc == 111) .or. (nproc == 112) 
     &   .or. (nproc == 119)) then
        if (ipermille == 10) taucut=3.e-2_dp
        if (ipermille ==  2) taucut=2.e-3_dp
      elseif ((nproc == 285)) then
        if (ipermille == 10) taucut=1.e-2_dp
        if (ipermille ==  2) taucut=1.e-3_dp
      else
        write(6,*) 'Pre-determined value of taucut not available'
        write(6,*) 'for nproc = ',nproc
        stop
      endif
      
      return
      end
      
