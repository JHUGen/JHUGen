      subroutine breitw_mod(x1,mminsq,mmaxsq,mflatsq,rmass,rwidth,
     & msq,wt)
      implicit none
      include 'constants.f'
c---- Given a number 0<x<1 generate a mass-squared msq and a weight wt 
c---- such that mminsq<msq<mmaxsq
c---- points are generated around resonance position rmass, but 
c---- breit-wigner should still be included in the matrix element
c     wt is the jacobian between integration in msq and integration in x1
c---- MODIFIED FROM BREITW:
c---- extra argument mflatsq that ensures that generation occurs
c---- according to BW around resonance mass for mflatsq<msq<mmaxsq
c---- but then flat in the region  mminsq<msq<mflatsq
      double precision x1,mminsq,mmaxsq,mflatsq,rmass,rwidth,msq,wt,z1
      double precision almin,almax,al,tanal,flatfrac,mmin,mflat
      include 'zerowidth.f'
      include 'leptcuts.f'

c--- fraction of events generated flat in invariant mass:
c--- this form is chosen such that as the photon cut is decreased,
c--- more events are generated this way due to copious photon radiation
c--- from leptons in Z decay
      flatfrac=(1d0-sqrt(gammpt/rmass))/2d0
      if (flatfrac .lt. 0d0) flatfrac=0d0
      if (mflatsq .lt. mminsq) then
        mflatsq=mminsq
        flatfrac=0d0
      endif

c--- sanity check
      if ((mflatsq .lt. mminsq) .or. (mflatsq .gt. rmass**2)) then
        write(6,*) 'Problem in breitw_mod: need mminsq<mflatsq<rmass**2'
        write(6,*) '   mminsq =',mminsq
        write(6,*) '  mflatsq =',mflatsq
        write(6,*) ' rmass**2 =',rmass**2
        stop
      endif

c--- in case the maximum msq is very small, just generate linearly for safety
      if ((mmaxsq .lt. rmass*1d-3) .and. (zerowidth .eqv. .false.)) then
        msq=mminsq+x1*(mmaxsq-mminsq)
        wt=mmaxsq-mminsq
        return
      endif

c--- if zerowidth, flat parameters make no difference
      if (zerowidth) then
          tanal=0d0
          almax=+pi/two
          almin=-pi/two
          msq=rmass**2+rmass*rwidth*tanal
c---- bw=(1d0+tanal**2)*rmass**2*rwidth**2
          wt=(almax-almin)*rmass*rwidth*(1d0+tanal**2)
          return
      endif
      
      if (x1 .gt. flatfrac) then
c--- generate according to BW
          z1=(x1-flatfrac)/(1d0-flatfrac)
          almin=datan((mflatsq-rmass**2)/rmass/rwidth)
          almax=datan((mmaxsq-rmass**2)/rmass/rwidth)
          al=(almax-almin)*z1+almin
          tanal=dtan(al)
          msq=rmass**2+rmass*rwidth*tanal
c---- bw=(1d0+tanal**2)*rmass**2*rwidth**2
          wt=(almax-almin)*rmass*rwidth*(1d0+tanal**2)
          wt=wt/(1d0-flatfrac)
      else
c---- generate flat in mass
        z1=x1/flatfrac
        mmin=sqrt(max(mminsq,1d-8))
        mflat=sqrt(max(mflatsq,1d-8))
        msq=(mmin+z1*(mflat-mmin))**2
        wt=2d0*(mflat-mmin)*sqrt(msq)
        wt=wt/flatfrac
      endif

      return
      end

