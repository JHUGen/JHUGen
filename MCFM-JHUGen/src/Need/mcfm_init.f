      subroutine mcfm_init(inputfile,workdir)
      implicit none
      include 'types.f'
************************************************************************
*                                                                      *
*  This routine should initialize any necessary variables and          *
*  perform the usual print-outs                                        *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
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
      include 'mpicommon.f'

c--- for APPLgrid
c      include 'ptilde.f'
c      include 'APPLinclude.f'
      real(dp):: rtsmin,p1ext(4),p2ext(4),p(mxpart,4),val
      integer:: j,k
      character*72 inputfile,workdir
      common/rtsmin/rtsmin
      common/pext/p1ext,p2ext
      data p/mxpart*3._dp,mxpart*4._dp,mxpart*0._dp,mxpart*5._dp/

!$omp threadprivate(/pext/)

* Welcome banner
      if (rank.eq.0) call banner

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

* Initialize parameter settings
      call mdata

      epinv=1.e1_dp
      epinv2=1.e1_dp

      call reader_input(inputfile,workdir)
      if (rank.ge.1) verbose=.false.

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
      cutoff=cutoff*(sqrts/2000._dp)**2
      rtsmin=min(rtsmin,sqrt(wsqmin+cutoff))
      rtsmin=min(rtsmin,sqrt(bbsqmin+cutoff))
      taumin=(rtsmin/sqrts)**2
      xmin=1.e-8_dp

      p1ext(4)=-half*sqrts
      p1ext(1)=0._dp
      p1ext(2)=0._dp
      p1ext(3)=-half*sqrts

      p2ext(4)=-half*sqrts
      p2ext(1)=0._dp
      p2ext(2)=0._dp
      p2ext(3)=+half*sqrts

* Set-up run name
      call setrunname(scale,facscale)

* Initialize all histograms
* Setup for histograms with irregular bins
      nirreg=0
      irregbin= (/ (.false.,j=1,maxhisto) /)

* npart=9 is a dummy value, to ensure that all histograms are included
      npart=9
      val=1.e-15_dp
      call nplotter(p,val,val**2,1)

      do j=1,mxpart
      do k=1,4
      p(j,k)=0._dp
      enddo
      enddo

* Initialize flag for photon fragmentation dipoles
      phot_dip(:)=.false.
      fragint_mode=.false.
* Initialize integer:: used in TensorReduction to zero
      TRtensorcontrol=0

      return
      end

