c--- J. Campbell, June 25th, 2008.

c--- A collection of routines that interface MCFM with the FROOT
c--- package for writing out ROOT ntuples.


      subroutine bookfill(tag,p,wt)
c--- This is the routine called by nplotter: when the output is to be
c--- initialized (tag="book") or added to (tag="plot")
      implicit none
      include 'constants.f'
      include 'maxwt.f'

      character tag*4
      double precision p(mxpart,4)
      double precision wt 

      if (.not.skipnt) then
        if     (tag .eq. 'book') then
          call FROOT_book
        elseif (tag .eq. 'plot') then
          call FROOT_fill(p,wt)
        endif
      endif

      return
      end

      subroutine NTfinalize
c--- This is the routine called at the end of the program's execution;
c--- it should finalize the output and close opened files, if necessary
      implicit none
      
      call PrintNT()
      call RootNTOutp()
      
      return
      end


      subroutine FROOT_book
      implicit none
      include 'npart.f'
      include 'mxdim.f'
      include 'scale.f'
      include 'facscale.f'
      include 'PDFerrors.f'
c--- Added to keep track of number of momenta entries to be filled
      include 'part.f'
      double precision scale_store,facscale_store
c--- Extra definitions to facilitate dummy call to lowint
      double precision dummy,wgt,r(mxdim),lowint
      integer ifill
      integer imaxmom,ipdf
      character*100 outfile
      character*255 runname
      integer lenocc
      logical first
      common/iarray/imaxmom,ipdf            
      common/runname/runname
      data first/.true./
      save first

c--- Need to ascertain the correct size for momenta n-tuples when this routine
c--- is called for the first time, achieved via a dummy call to lowint
      if (first) then      
        do ifill=1,mxdim
          r(ifill)=0.5d0
        enddo
c--- Be careful that dynamic scale choices aren't ruined
c--- (in versions 5.1 and before, this occured when calling lowint)
        scale_store=scale
        facscale_store=facscale
        dummy=lowint(r,wgt)
        scale=scale_store
        facscale=facscale_store
        
        imaxmom=npart
        if ((part.eq.'real').or.(part.eq.'tota').or.(part.eq.'todk'))
     .    imaxmom=imaxmom+1

        first=.false.
      endif

c--- determine if we need space in array to store PDF weights (ipdf)
      ipdf=0
      if (PDFerrors) then
        ipdf= maxPDFsets
      endif
           
c --- Create an empty ntuple with the usual file name:
      outfile=runname(1:lenocc(runname))//'.root'
      call InitRootNT(outfile,'RECREATE')
      write(6,*) '<----- Initialized FROOT n-tuple ----->' 
      call flush(6)
            
      return
      end



      subroutine FROOT_fill(p,wt)
      implicit none
      include 'constants.f'
      include 'wts_bypart.f'
      include 'PDFerrors.f'      
      double precision p(mxpart,4)
      double precision wt 
c--- Extra common block to carry the information about maximum momenta entries
      integer imaxmom,ipdf
      common/iarray/imaxmom,ipdf
      integer i
c--- assume at most 10 final state particles and 60 PDF sets
      real pfill(105)
      character*3 labelE,labelx,labely,labelz
      character*5 labelPDF
      logical first
c--- force pfill to be allocated statically
      common/pfillcommon/pfill      
      data first/.true./
      save first
      
c--- If the event weight is zero, don't bother to add the n-tuple
      if (wt .eq. 0d0) then
        return
      endif

c--- On the first call, must set-up all branches
      if (first) then
        do i=1,imaxmom
          write(labelx,72) i+2
          write(labely,73) i+2
          write(labelz,74) i+2
          write(labelE,71) i+2
c---    first, momenta
          call AddNTBranch(pfill(4*(i-1)+1),labelx)
          call AddNTBranch(pfill(4*(i-1)+2),labely)
          call AddNTBranch(pfill(4*(i-1)+3),labelz)
          call AddNTBranch(pfill(4*(i-1)+4),labelE)
        enddo
c---    now, the weights
        call AddNTBranch(pfill(imaxmom*4+1),'wt_ALL')
        call AddNTBranch(pfill(imaxmom*4+2),'wt_gg')
        call AddNTBranch(pfill(imaxmom*4+3),'wt_gq')
        call AddNTBranch(pfill(imaxmom*4+4),'wt_qq')
        call AddNTBranch(pfill(imaxmom*4+5),'wt_qqb')
c---    then include PDF errors if necessary
        if (PDFerrors) then
          do i=1,ipdf
          write(labelPDF,75) i
          call AddNTBranch(pfill(imaxmom*4+5+i),labelPDF)
            enddo
        endif
        first=.false.
      endif
      

      do i=1,imaxmom
c---    set up single precision variables for momenta
        pfill(4*(i-1)+1)=sngl(p(i+2,1))
        pfill(4*(i-1)+2)=sngl(p(i+2,2))
        pfill(4*(i-1)+3)=sngl(p(i+2,3))
        pfill(4*(i-1)+4)=sngl(p(i+2,4))
      enddo

c---    set up single precision variables for the event weights
      pfill(imaxmom*4+1)=sngl(wt)
      pfill(imaxmom*4+2)=sngl(wt_gg)
      pfill(imaxmom*4+3)=sngl(wt_gq)
      pfill(imaxmom*4+4)=sngl(wt_qq)
      pfill(imaxmom*4+5)=sngl(wt_qqb)

c--- include PDF errors if necessary
      if (PDFerrors) then
        do i=1,ipdf
        pfill(imaxmom*4+5+i)=sngl(wt*PDFwgt(i)/PDFwgt(0))
        enddo
      endif

c--- now that all branches are set-up, fill all entries
      call FillNTBranch('all')
      
      return
      
c--- format for labels
   71 format('E_',i1)
   72 format('px',i1)
   73 format('py',i1)
   74 format('pz',i1)
   75 format('PDF',i2.2)
      
      end






      

c--- dummy routines, as in dsw_dummy

      subroutine dswhbook(n,titlex,dx,xmin,xmax)
      implicit none
      integer n
      character titlex*8
      real*8 dx,xmin,xmax

      call dsw_error

      return
      end
c

      subroutine dswhfill(n,var,wgt)
      implicit none
      integer n
      real*8 var,wgt

      call dsw_error

      return
      end

c

      subroutine dsw_error
      implicit none

      write(6,*) 'This version of MCFM has not been compiled for'
      write(6,*) 'FROOT output only; DSW-style histograms are not'
      write(6,*) 'available.'
      write(6,*)
      stop
      end
c      

