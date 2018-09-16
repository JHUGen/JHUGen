      subroutine coupling
      implicit none
      include 'types.f'
c--- initialize electroweak couplings and calculate alpha-s; this
c--- must be called at the beginning of "chooser" to enable it to
c--- set up all the variables (e.g. twidth). Once nflav is set,
c--- alpha-s will be determined again (in coupling2).
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'scale.f'
      include 'nlooprun.f'
      include 'ewinput.f'
      include 'nflav.f'
      include 'b0.f'
      include 'kpart.f'
      include 'verbose.f'
      include 'mpicommon.f'
      integer:: i
      real(dp):: aemmz,alphas,cmass,bmass,lotopdecaywidth
      character*3 inlabel(10)
      character*5 tworder
      common/qmass/cmass,bmass
      common/em/aemmz

c--- blank out labels that indicate input parameters
      do i=1,10
        inlabel(i)='   '
      enddo
      inlabel(3)='(+)'
      inlabel(4)='(+)'
      inlabel(8)='(+)'

      if (ewscheme == -1) then
c--- This is the MCFM default, corresponding to an effective
c--- field theory approach valid for scales below the top-mass
C--- (see Georgi, Nucl. Phys. B 363 (1991) 301).
c--- There are 4 inputs here instead of the usual 3 ...
         Gf = Gf_inp
         aemmz  = aemmz_inp
         wmass  = wmass_inp
         zmass  = zmass_inp
         inlabel(5)='(+)'
         inlabel(6)='(+)'
         inlabel(1)='(+)'
         inlabel(2)='(+)'
         inlabel(8)='   '
c--- ... and as result, both xw and mtop are derived
c--- (using wmass=zmass*sqrt(rho)*cos(theta_w)
c---    and rho=1+3._dp*aemmz/16._dp/pi/xw*(mt/wmass)**2 )
         xw  = fourpi*aemmz/(8._dp*wmass**2*Gf/rt2)
         mt  = sqrt(16._dp*pisq/3._dp/rt2/Gf*(
     &          wmass**2/zmass**2/(1._dp-xw)-1._dp))

      elseif (ewscheme == 0) then
c------------------------------------------------------------
c     option=0 : MadEvent default (= AlpGen with iewopt=2)
c------------------------------------------------------------

c-- equal to the input values
         xw  = xw_inp
         aemmz  = aemmz_inp
         zmass  = zmass_inp
         inlabel(7)='(+)'
         inlabel(6)='(+)'
         inlabel(1)='(+)'
c-- derived
         wmass  = zmass * sqrt( One - xw )
         Gf = aemmz * Fourpi/xw/(8._dp*wmass**2/Rt2)

      elseif (ewscheme == 1) then
c-----------------------------------------------------
c     option=1 : LUSIFER and AlpGen (iewopt=3) default
c-----------------------------------------------------

c-- equal to the input values
         zmass  = zmass_inp
         wmass  = wmass_inp
         Gf = Gf_inp
         inlabel(1)='(+)'
         inlabel(2)='(+)'
         inlabel(5)='(+)'
c-- derived
         xw  = One-(wmass/zmass)**2
         aemmz  = Rt2*Gf*wmass**2*xw/pi

      elseif (ewscheme == 2) then
c-------------------------------------------------------------------
c     option=2 : W and Z mass are derived from couplings
c-------------------------------------------------------------------

c-- equal to the input values
         Gf = Gf_inp
         aemmz  = aemmz_inp
         xw  = xw_inp
         inlabel(5)='(+)'
         inlabel(6)='(+)'
         inlabel(7)='(+)'
c-- derived
         wmass  = sqrt(aemmz*pi/xw/Gf/Rt2)
         zmass  = wmass/sqrt(One-xw)

      elseif (ewscheme == 3) then
c-----------------------------------------------------------------
c     option=3 : USER choice : you should know what you're doing!!
c-----------------------------------------------------------------
         Gf = Gf_inp
         aemmz  = aemmz_inp
         xw  = xw_inp
         wmass  = wmass_inp
         zmass  = zmass_inp
         inlabel(5)='(+)'
         inlabel(6)='(+)'
         inlabel(7)='(+)'
         inlabel(1)='(+)'
         inlabel(2)='(+)'

      else
         write(6,*) 'ewscheme=',ewscheme,' is not a valid input.'
         stop
      endif

c--- Now set up the other derived parameters
      gwsq=fourpi*aemmz/xw
      esq=gwsq*xw
      gw=sqrt(gwsq)
      call couplz(xw)

c--- Calculate the appropriate Higgs vacuum expectation value.
c--- This vevsq is defined so that gwsq/(4*wmass**2)=Gf*rt2=1/vevsq
c--- (ie differs from definition in ESW)
      vevsq=1._dp/rt2/Gf

c--- set up the beta-function
      b0=(xn*11._dp-2._dp*nflav)/6._dp

c--- initialize the pdf set
      nlooprun=0
      call pdfwrap      

      cmass=sqrt(mcsq)
      bmass=sqrt(mbsq)
      musq=scale**2
 
c--- set the number of loops to use in the running of alpha_s
c--- if it hasn't been set by pdfwrap already
      if (nlooprun == 0) then
        if (kpart==klord) then
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

c--- Set-up twidth, using LO formula everywhere
      twidth=lotopdecaywidth(mt,mb,wmass,wwidth)     
      tworder='(LO) '
c      if ( (kpart==ktodk) .or. (mykpart==ktodk) ) then
c        twidth=twidth*nlotopdecaywidth(mt,mb,wmass)
c        tworder='(NLO) '
c      endif
c      write(6,*) 'twidth LO',twidth
c      write(6,*) 'twidth NLO',twidth*topwidth(mt,wmass)
c      write(6,*) 'twidth (NLO-LO)/LO = ',topwidth(mt,wmass)-1._dp

      if (rank.eq.0) then
      if (verbose) then
      write(6,*) '************** Electroweak parameters **************'
      write(6,*) '*                                                  *'
      write(6,75) 'zmass',inlabel(1),zmass,'wmass',inlabel(2),wmass
      write(6,75) 'zwidth',inlabel(3),zwidth,'wwidth',inlabel(4),wwidth
      write(6,76) 'Gf',inlabel(5),gf,'1/aemmz',inlabel(6),1._dp/aemmz
      write(6,75) 'xw',inlabel(7),xw,'mtop',inlabel(8),mt
      write(6,75) 'gwsq',inlabel(9),gwsq,'esq',inlabel(10),esq
      write(6,77) 'top width',twidth,tworder
      write(6,78) 'mb',mb,'mc',mc
      write(6,*) '*                                                  *'
      write(6,*) '* Parameters marked (+) are input, others derived  *'
      write(6,*) '****************************************************'
      endif
      endif

   75 format(' * ',a6,a3,f13.7,3x,a7,a3,f12.7,'  *')
   76 format(' * ',a6,a3,d13.6,3x,a7,a3,f12.7,'  *')
   77 format(' * ',a9,f13.7,1x,a5,19x,'  *')
   78 format(' * ',a5,4x,f13.7,6x,a4,2x,f13.7,'  *')
      
      return
      end


      block data wsalam1
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ewcharge.f'
      data Q(-5)/+0.333333333333333_dp/
      data Q(-4)/-0.666666666666667_dp/
      data Q(-3)/+0.333333333333333_dp/
      data Q(-2)/-0.666666666666667_dp/
      data Q(-1)/+0.333333333333333_dp/
      data Q(0)/+0._dp/
      data Q(+1)/-0.333333333333333_dp/
      data Q(+2)/+0.666666666666667_dp/
      data Q(+3)/-0.333333333333333_dp/
      data Q(+4)/+0.666666666666667_dp/
      data Q(+5)/-0.333333333333333_dp/
      data tau/1._dp,-1._dp,1._dp,-1._dp,1._dp,0._dp,-1._dp,1._dp,-1._dp,1._dp,-1._dp/
      end 



