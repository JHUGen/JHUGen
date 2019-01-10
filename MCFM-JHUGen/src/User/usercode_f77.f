      logical function userincludedipole(nd, ppart, mcfm_result)
************************************************************************
*   "user" routine that gets called to allow the user to veto events
*   that might otherwise pass the MCFM cuts. 
*   
*   It can be used, e.g. to force MCFM to generate only events that are
*   above some large HT threshold, which comes in useful when trying to
*   get precision on the tails of some distributions.
*   
*   Variables passed are
*
*   - nd:           index of the dipole
*   - ppart:        momenta of incoming and outgoing particles
*   - mcfm_result:  the decision that was taken by mcfm about whether to keep
*                   this event
*
*   NB: if mcfm_result is .false., MCFM's jet finding may not have run, 
*       and anything in MCFM that depends on it (e.g. some dynamic scales) 
*       will misbehave. It is safer not to return .true. when
*       mcfm_result=.false.
************************************************************************
      implicit none
      include 'constants.f'
      include 'npart.f'
      integer          nd
      double precision ppart(mxpart,4)
      logical          mcfm_result

c--- take the MCFM result as the default return value      
      userincludedipole = mcfm_result

c--- example of code (using compact f90 notation) that places a cut on HT
c      if (sum(sqrt(sum(ppart(3:2+npart,1:2)**2,dim=2))) < 500d0) then
c         userincludedipole = .false.
c      end if
      end 



      subroutine userplotter(pjet, wt,wt2, nd)
************************************************************************
*   Subroutine that is called to allow the user to bin their own
*   histograms
*
*   Variables passed to this routine:
*   
*        p:  4-momenta of incoming partons(i=1,2), outgoing leptons and
*            jets in the format p(i,4) with the particles numbered
*            according to the input file and components labelled by
*            (px,py,pz,E)?
*   
*       wt:  weight of this event
*   
*      wt2:  weight^2 of this event
*   
*       nd:  an integer specifying the dipole number of this contribution
*            (if applicable), otherwise equal to zero
************************************************************************
      implicit none
      include 'constants.f'
      include 'ptilde.f'
      include 'npart.f'
      include 'nplot.f'
      double precision pjet(mxpart,4)
      double precision wt,wt2
      integer          nd
c--- if you don't need extra histograms, you don't need to touch this code.
c--- Follow the commented example below to see what kind of things
c--- you might want to do.

c      integer          iplot 
c      double precision ht, htjet
c      logical          first
c      data first/.true./
c      save first
c      character*4 tag
c
c      if (first) then
c         tag   = "book"
c         first = .false.
c      else                
c         tag = "plot"
c      end if
c      
c      ! each histogram has an index; nextnplot from the nplot common block
c      ! holds the info on the next available common block
c      iplot = nextnplot
c
c      ! ht is scalar pt sum of all outgoing particles -- accessed from
c      ! the ptilde common block
c      ht    = sum(sqrt(sum(ptilde   (nd,3:2+npart,1:2)**2,dim=2)))
c      
c      call bookplot(iplot,tag,'UserHT',ht,wt,wt2,0d0,1000d0,50d0,'lin'); 
c      iplot = iplot + 1
      
      end 



      subroutine userwriteinfo(unitno,comment_string,xsec,xsec_err,itno)
************************************************************************
*   subroutine that gets called after MCFM has written its comments to 
*   one of the output files. It allows the user to write their own comments
*   to that same file
*
*   Variables passed to this routine:
*
*   - unitno: the unit number to which output is being sent
*   - comment_string: a comment character that precedes each line of output
*   - xsec, xsec: the cross section and its error (in case you care!)
*   - itno: the iteration number (0 at the end of the last iteration)
************************************************************************
      implicit none
      integer          unitno
      character*2      comment_string
      double precision xsec, xsec_err
      integer          itno
  
c      write(6,*) "have reached iteration number", itno
c      write(unitno,"(a,a)") comment_string, "ht cut > 500 GeV"
      end

c---- SSbegin
c----------------------------------------------------------
      subroutine userhistofin(xsec,xsec_err,itno,itmx)
c     This function allows for extra user-defined operations 
c     at the end of each iteration (itno>0) and at the end of 
c     the run of the program (itno=0).
      implicit none
      integer itno,itmx
      double precision xsec,xsec_err
      end
c---- SSend
