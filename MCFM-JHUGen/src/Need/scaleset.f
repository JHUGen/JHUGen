      subroutine scaleset(rscalestart,fscalestart,p)
      implicit none
      include 'types.f'
c--- wrapper routine to set a dynamic scale; please refer to individual
c--- routines for exact definitions of the scales;
c--- upgraded 3/2016 to convert string to integer on first call, for speed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'frag.f'
      include 'nlooprun.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'stopscales.f'
      include 'masses.f'
      include 'first.f'
      include 'mpicommon.f'
      real(dp):: rscalestart,fscalestart,p(mxpart,4),mu0,
     & alphas
      integer, parameter ::
     & khmass=1,kwmass=2,kzmass=3,kmt=4,km34=5,km345=6,km3456=7,
     & kMsqpt34sq=8,kMsqpt345sq=9,kMsqpt5sq=10,kMsqptj1sq=11,kMsqsumptjsq=12,
     & km34sqsumptjsq=13,kptphoton=14,kHT=15,kddis=16,kmVpmH=17,kshat=18
      integer, save :: scaleindex
!$omp threadprivate(scaleindex)

! Initialization and write-out
      if (first) then
        if (rank == 0) then
        write(6,*)
        write(6,*)'************** Dynamic scale choice ****************'
        write(6,*)'*                                                  *'
        write(6,*)'*                 RENORMALIZATION                  *'
        write(6,45) ' mu_ren  =',rscalestart,dynstring
        write(6,*)'*                                                  *'
        write(6,*)'*                  FACTORIZATION                   *'
        write(6,45) ' mu_fac  =',fscalestart,dynstring
        if (frag) then
        write(6,*)'*                                                  *'
        write(6,*)'*                  FRAGMENTATION                   *'
        write(6,45) ' mu_frag =',frag_scalestart,dynstring
        endif
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        endif
        first=.false.
        if     ((dynstring == 'mh') .or. (dynstring == 'mH')
     &     .or. (dynstring == 'Mh') .or. (dynstring == 'MH')) then
          scaleindex=khmass
        elseif ((dynstring == 'mw') .or. (dynstring == 'mW')
     &     .or. (dynstring == 'Mw') .or. (dynstring == 'MW')) then
          scaleindex=kwmass
        elseif ((dynstring == 'mz') .or. (dynstring == 'mZ')
     &     .or. (dynstring == 'Mz') .or. (dynstring == 'MZ')) then
          scaleindex=kzmass
        elseif ((dynstring == 'mt') .or. (dynstring == 'mT')
     &     .or. (dynstring == 'Mt') .or. (dynstring == 'MT')) then
          scaleindex=kmt
        elseif (dynstring == 'm(34)') then
          scaleindex=km34
        elseif (dynstring == 'm(345)') then
          scaleindex=km345
        elseif (dynstring == 'm(3456)') then
          scaleindex=km3456
        elseif (dynstring == 'sqrt(M^2+pt34^2)') then
          scaleindex=kMsqpt34sq
        elseif (dynstring == 'sqrt(M^2+pt345^2)') then
          scaleindex=kMsqpt345sq
        elseif (dynstring == 'sqrt(M^2+pt5^2)') then
          scaleindex=kMsqpt5sq
        elseif (dynstring == 'sqrt(M^2+ptj1^2)') then
          scaleindex=kMsqptj1sq
        elseif (dynstring == 'sqrt(M^2+sumptj^2)') then
          scaleindex=kMsqsumptjsq
        elseif (dynstring == 'sqrt(m(34)^2+sumptj^2)') then
          scaleindex=km34sqsumptjsq
        elseif (dynstring == 'pt(photon)') then
          scaleindex=kptphoton
        elseif (dynstring == 'HT') then
          scaleindex=kHT
        elseif (dynstring == 'DDIS') then
          scaleindex=kddis
        elseif (dynstring == 'mV+mH') then
          scaleindex=kmVpmH
        elseif (dynstring == 's-hat') then
          scaleindex=kshat
        else
          write(6,*) 'Dynamic scale choice not recognized'
          write(6,*) '   dynamicscale = ',dynstring
          stop
        endif
      endif

! Set scale
      if     (scaleindex == khmass) then
        mu0=hmass
      elseif (scaleindex == kwmass) then
        mu0=wmass
      elseif (scaleindex == kzmass) then
        mu0=zmass
      elseif (scaleindex == kmt) then
        mu0=mt
      elseif (scaleindex == km34) then
        call scaleset_m34(p,mu0)
      elseif (scaleindex == km345) then
        call scaleset_m345(p,mu0)
      elseif (scaleindex == km3456) then
        call scaleset_m3456(p,mu0)
      elseif (scaleindex == kMsqpt34sq) then
        call scaleset_Msqpt34sq(p,mu0)
      elseif (scaleindex == kMsqpt345sq) then
        call scaleset_Msqpt345sq(p,mu0)
      elseif (scaleindex == kMsqpt5sq) then
        call scaleset_Msqpt5sq(p,mu0)
      elseif (scaleindex == kMsqptj1sq) then
        call scaleset_Msqptj1sq(p,mu0)
      elseif (scaleindex == kMsqsumptjsq) then
        call scaleset_Msqsumptjsq(p,mu0)
      elseif (scaleindex == km34sqsumptjsq) then
        call scaleset_m34sqsumptjsq(p,mu0)
      elseif (scaleindex == kptphoton) then
        call scaleset_ptphoton(p,mu0)
      elseif (scaleindex == kHT) then
        call scaleset_HT(p,mu0)
      elseif (scaleindex == kDDIS) then
        call scaleset_ddis(p,mu0)
      elseif (scaleindex == kmVpmH) then
        call scaleset_mVpmH(p,mu0)
      elseif (scaleindex == kshat) then
        call scaleset_shat(p,mu0)
      else
        write(6,*) 'Dynamic scale choice not recognized'
        write(6,*) '   scaleindex = ',scaleindex
        stop
      endif
      
      scale=rscalestart*mu0
      facscale=fscalestart*mu0
      frag_scale=frag_scalestart*mu0
          
c--- catch absurdly large and small scales      
      if  (scale > 100000._dp) scale=60000._dp
      if  (facscale > 100000._dp) facscale=60000._dp
      if  (frag_scale > 990._dp) frag_scale=900._dp
      if  (scale < 1._dp) scale=1._dp
      if  (facscale < 1._dp) facscale=1._dp
      if  (frag_scale < 1._dp) frag_scale=1._dp

c--- run alpha_s
      as=alphas(scale,amz,nlooprun)
      
      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as
      musq=scale**2

c--- these are additional scales used in the t-channel single top + b routines
      renscale_H=scale
      renscale_L=scale
      facscale_H=facscale
      facscale_L=facscale
      as_H=as
      as_L=as
            
      return

 45   format(1x,'* ',a15,f6.2,' x ',a24,' *')

      end
      
