      subroutine gen4_intf(r,p,wt,*)
      implicit none
      include 'types.f'
c--- Generates phase space for gg -> VV | gg-> H -> VV interference processes.
c--- If the mass of the Higgs boson is above (or close) to the minimum
c--- invariant mass required by the cuts (m34min and m56min),
c--- the phase space is generated using a B.W. for the Higgs boson (gen4h.f);
c--- if the Higgs mass is below this threshold then the phase space is
c--- generated using the standard VV routine gen4.f

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'mxdim.f'
      include 'limits.f'
      include 'first.f'
      logical, save::useHiggsBW
      real(dp):: r(mxdim),p(mxpart,4),wt,threshold

c--- initialization on first call
      if (first) then
        threshold=wsqmin+bbsqmin
c--- generate using B.W. if Higgs mass above or close to threshold
        if (hmass > threshold-5._dp*hwidth) then
c--- otherwise generate without B.W.
          useHiggsBW=.true.
        else
          useHiggsBW=.false.
        endif
        write(6,*)
        write(6,*) 'gen4_intf: useHiggsBW = ',useHiggsBW
        first=.false.
      endif

c--- switch between gen4h and gen4
      if (useHiggsBW) then
        call gen4h(r,p,wt,*999)
      else
        call gen4(r,p,wt,*999)
      endif

      return

 999  wt=0._dp
      return 1

      end

