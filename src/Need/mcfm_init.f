      subroutine mcfm_init(inputfile,workdir)
************************************************************************
*                                                                      *
*  This routine should initialize any necessary variables and          *
*  perform the usual print-outs                                        *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'cutoff.f'
      include 'limits.f'
      include 'npart.f'
      include 'phasemin.f'
      include 'facscale.f'
      include 'scale.f'
      include 'verbose.f'
      include 'phot_dip.f'
      include 'includect.f'
      include 'TRtensorcontrol.f'
      include 'qlfirst.f'
      include 'frag.f'
      include 'histo.f'
      include 'irregbins_incl.f'
      include 'energy.f'
C -- GZ
      include 'first_time.f'
      include 'TRmaxindex.f'
      include 'masses.f'
      include 'nflav.f'
      include 'nores.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'ehisto.f'

c--- for APPLgrid
c      include 'ptilde.f'
c      include 'APPLinclude.f'
      double precision rtsmin,p1ext(4),p2ext(4),p(mxpart,4),val
      integer j,k
      character*72 inputfile,workdir
      common/rtsmin/rtsmin
      common/pext/p1ext,p2ext
      data p/mxpart*3d0,mxpart*4d0,mxpart*0d0,mxpart*5d0/
!$omp threadprivate(/pext/)

* Welcome banner
      call banner

*Initialise data in commonblock /pvmaxindex/
      maxcindex=3
      maxdindex=4
      maxeindex=5
* Flag to control whether QCDLoop needs to be called yet, if appropriate
      qlfirst=.true.
* Set all plots to linear scale by default
      linlog(:)='lin'
* Set up ebook initial settings
      IHISTOMATCH(:)=0
      ICOUNTHISTO=0
* Set up mbook initial settings
      book(:)='NO'

************************************************************************
*     Masses, widths and initial-state flavour information             *
************************************************************************
c--- if true, nores removes all of the gg contribution
      nores=.false.
c--- Masses: note that "mtausq" is typically used throughout the
c--- program to calculate couplings that depend on the mass, while
c--- "mtau" is the mass that appears in the rest of the matrix
c--- elements and phase space (and may be set to zero in the program,
c--- depending on the process number)

      mtau=1.777d0
      mtausq=3.157729d0
c----   Note: after v5.6, the masses for top, bottom and charm quarks
c----         are set in the input file

c---  Widths: note that the top width is calculated in the program
c---  The W width of 2.1054 is derived using the measured BR of
c---    10.80 +/- 0.09 % (PDG) and the LO partial width calculation
c---    for Mw=80.398 GeV
      wwidth=2.1054d0
      zwidth=2.4952d0
      tauwidth=2.269d-12
c--- Number of active flavours in the initial state: this parameter
c--- may be changed in the program for some processes
      nflav=5
c--- Masses below here are currently unused
      md=5d-3
      mu=5d-3
      ms=1d-1
      mel=0.510997d-3
      mmu=0.105658389d0
*
************************************************************************
*     Dim. Reg. parameter epsilon, used for checking the proper        *
*      operation of the NLO code in the program                        *
************************************************************************
*
      epinv=1d3
      epinv2=1d3


      call reader_input(inputfile,workdir)

      first_time = .true.

      if (verbose) then
      write(6,*)
      write(6,*) '****************************************'
      write(6,*) '*     Cross section in femtobarns      *'
      write(6,*) '****************************************'
      write(6,*)
      endif

* Counter-terms for radiation in top decay should be included
      includect=.true.

* Set-up incoming beams and PS integration cut-offs
c--- Note: version 6.4 onwards, scale cutoff with c.o.m. energy
      cutoff=cutoff*(sqrts/2000d0)**2
      rtsmin=min(rtsmin,dsqrt(wsqmin+cutoff))
      rtsmin=min(rtsmin,dsqrt(bbsqmin+cutoff))
      taumin=(rtsmin/sqrts)**2
      xmin=1d-8

      p1ext(4)=-half*sqrts
      p1ext(1)=0d0
      p1ext(2)=0d0
      p1ext(3)=-half*sqrts

      p2ext(4)=-half*sqrts
      p2ext(1)=0d0
      p2ext(2)=0d0
      p2ext(3)=+half*sqrts

* Set-up run name
      call setrunname(scale,facscale)

* Initialize all histograms
* Setup for histograms with irregular bins
      nirreg=0
      irregbin= (/ (.false.,j=1,maxhisto) /)

* npart=9 is a dummy value, to ensure that all histograms are included
      npart=9
      val=1d-15
      call nplotter(p,val,val**2,1)

      do j=1,mxpart
      do k=1,4
      p(j,k)=0d0
      enddo
      enddo

* Initialize flag for photon fragmentation dipoles
      phot_dip(:)=.false.
      fragint_mode=.false.
* Initialize integer used in TensorReduction to zero
      TRtensorcontrol=0

      return
      end

