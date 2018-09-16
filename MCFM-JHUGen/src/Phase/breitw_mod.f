      subroutine breitw_mod(x1,mminsq,mmaxsq,mflatsq,rmass,rwidth,
     & msq,wt)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
c---- Given a number 0<x<1 generate a mass-squared msq and a weight wt 
c---- such that mminsq<msq<mmaxsq
c---- points are generated around resonance position rmass, but 
c---- breit-wigner should still be included in the matrix element
c     wt is the jacobian between integration in msq and integration in x1
c---- MODIFIED FROM BREITW:
c---- extra argument mflatsq that ensures that generation occurs
c---- according to BW around resonance mass for mflatsq<msq<mmaxsq
c---- but then flat in the region  mminsq<msq<mflatsq
      real(dp):: x1,mminsq,mmaxsq,mflatsq,rmass,rwidth,msq,wt,z1
      real(dp):: almin,almax,al,tanal,flatfrac,mmin,mflat
      include 'zerowidth.f'
      include 'leptcuts.f'

c--- fraction of events generated flat in invariant mass:
c--- this form is chosen such that as the photon cut is decreased,
c--- more events are generated this way due to copious photon radiation
c--- from leptons in Z decay
      flatfrac=(1._dp-sqrt(gammpt/rmass))/2._dp
      if (flatfrac < 0._dp) flatfrac=0._dp
      if (mflatsq < mminsq) then
        mflatsq=mminsq
        flatfrac=0._dp
      endif

c--- sanity check
      if ((mflatsq < mminsq) .or. (mflatsq > rmass**2)) then
        write(6,*) 'Problem in breitw_mod: need mminsq<mflatsq<rmass**2'
        write(6,*) '   mminsq =',mminsq
        write(6,*) '  mflatsq =',mflatsq
        write(6,*) ' rmass**2 =',rmass**2
        stop
      endif

c--- in case the maximum msq is very small, just generate linearly for safety
      if ((mmaxsq < rmass*1.e-3_dp) .and. (zerowidth .eqv. .false.)) then
        msq=mminsq+x1*(mmaxsq-mminsq)
        wt=mmaxsq-mminsq
        return
      endif

c--- if zerowidth, flat parameters make no difference
      if (zerowidth) then
          tanal=0._dp
          almax=+pi/two
          almin=-pi/two
          msq=rmass**2+rmass*rwidth*tanal
c---- bw=(1._dp+tanal**2)*rmass**2*rwidth**2
          wt=(almax-almin)*rmass*rwidth*(1._dp+tanal**2)
          return
      endif
      
      if (x1 > flatfrac) then
c--- generate according to BW
          z1=(x1-flatfrac)/(1._dp-flatfrac)
          almin=atan((mflatsq-rmass**2)/rmass/rwidth)
          almax=atan((mmaxsq-rmass**2)/rmass/rwidth)
          al=(almax-almin)*z1+almin
          tanal=tan(al)
          msq=rmass**2+rmass*rwidth*tanal
c---- bw=(1._dp+tanal**2)*rmass**2*rwidth**2
          wt=(almax-almin)*rmass*rwidth*(1._dp+tanal**2)
          wt=wt/(1._dp-flatfrac)
      else
c---- generate flat in mass
        z1=x1/flatfrac
        mmin=sqrt(max(mminsq,1.e-8_dp))
        mflat=sqrt(max(mflatsq,1.e-8_dp))
        msq=(mmin+z1*(mflat-mmin))**2
        wt=2._dp*(mflat-mmin)*sqrt(msq)
        wt=wt/flatfrac
      endif

      return
      end

