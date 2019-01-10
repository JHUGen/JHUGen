      subroutine coupling2
c--- this routine calculates alpha-s using the now-determined
c--- value of nflav and makes CKM matrix diagonal if necessary
      implicit none
      include 'constants.f'
      include 'masses.f'
c      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'verbose.f'
      include 'nlooprun.f'
c      include 'process.f'
c      include 'ewinput.f'
      include 'nflav.f'
      include 'b0.f'
      include 'dynamicscale.f'
      include 'stopscales.f'
      include 'fourthgen.f'
      include 'couple.f'
      include 'part.f'
      include 'nproc.f'
      double precision alphas,cmass,bmass
      double precision Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,
     &             Vcd,Vcs,Vcb
      common/qmass/cmass,bmass
c      common/em/aemmz
c      common/mypart/mypart

c--- set up the beta-function
      b0=(xn*11d0-2d0*nflav)/6d0

c--- initialize the pdf set
      nlooprun=0
      call pdfwrap      

      musq=scale**2
c--- set up masses used in running of alphas (alfamz.f)
      cmass=dsqrt(mcsq)
      if (fourthgen) then
        continue ! normal b mass already set in chooser.f
      else
        bmass=dsqrt(mbsq) ! use the mass specified in the input file
      endif
      
      if (nflav .lt. 5) then
        bmass=1001d0
      endif
      if (nflav .lt. 4) then
        cmass=1000d0
      endif
 
c--- set the number of loops to use in the running of alpha_s
c--- if it hasn't been set by pdfwrap already
      if (nlooprun .eq. 0) then
        if (part .eq. 'lord') then
          nlooprun=1
        else
          nlooprun=2
        endif
      endif

c--- initialize alpha_s
      as=alphas(abs(scale),amz,nlooprun)

      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as

c--- if we're doing W + jets, automatically make the CKM matrix
c--- diagonal since we're not interested in these small effects   
      if ((nproc .eq. 11) .or. (nproc .eq. 16) .or.
     .    (nproc .eq. 22) .or. (nproc .eq. 27)) then
        Vud=1d0
        Vus=0d0
        Vub=0d0
        Vcd=0d0
        Vcs=1d0
        Vcb=0d0
      endif

      if (verbose) then
      write(6,*)
      write(6,*) '***************** CKM mixing matrix ****************'
      write(6,*) '*                                                  *'
      write(6,47) Vud,Vus,Vub
      write(6,48) Vcd,Vcs,Vcb
      write(6,*) '****************************************************'
      if ((nproc .eq. 11) .or. (nproc .eq. 16) .or.
     .    (nproc .eq. 22) .or. (nproc .eq. 27)) then
      write(6,*) '* Forced to be diagonal for simplicity in W + jets *'
      write(6,*) '****************************************************'
      endif
 47   format(' *      Vud=',g10.5,'Vus=',g10.5,'Vub=',g10.5,'  *')
 48   format(' *      Vcd=',g10.5,'Vcs=',g10.5,'Vcb=',g10.5,'  *')
      endif      

c--- special write-out for stop+b case
      if ((verbose) .and. (initrenscale_L .gt. 0d0)) then
      write(6,*)
      write(6,*) '************* Strong coupling, alpha_s  ************'
      write(6,*) '*                                                  *'
      write(6,49) 'alpha_s (zmass)    ',amz
      write(6,49) 'alpha_s (hvy scale)',as_H
      write(6,49) 'alpha_s (lgt scale)',as_L
      write(6,50) ' (using ',nlooprun,'-loop running of alpha_s)'  
      write(6,*) '****************************************************'
      return
      endif

      if ((verbose) .and. (scale .gt. 0d0)) then      
      write(6,*)
      write(6,*) '************* Strong coupling, alpha_s  ************'
      write(6,*) '*                                                  *'
      if (dynamicscale .eqv. .false.) then
      write(6,49) 'alpha_s (scale)',gsq/fourpi
      write(6,49) 'alpha_s (zmass)',amz
      else
      write(6,*) '*  Dynamic scale - alpha_s changed event-by-event  *'
      write(6,49) 'alpha_s (zmass)',amz
      endif
      write(6,50) ' (using ',nlooprun,'-loop running of alpha_s)'  
      write(6,*) '****************************************************'
 49   format(' *  ',a20,f12.8,16x,'*')
 50   format(' *  ',6x,a8,i1,a25,8x,'*')
      endif
      
      return
      end
