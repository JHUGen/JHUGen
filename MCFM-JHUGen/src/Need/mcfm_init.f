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

* Flag to control whether QCDLoop needs to be called yet, if appropriate
      qlfirst=.true.

      return
      end
            
