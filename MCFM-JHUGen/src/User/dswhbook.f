c --- Dave Waters, 28.05.2001
c --- =======================
c --- Provide routines that are called from mcfm and fill
c --- histograms in an hbook file. It is more convenient 
c --- and flexible to obtain hbook histograms
c --- in this way rather than parse the standard mcfm
c --- output files.

      subroutine dswhbook(n,titlex,dx,xmin,xmax)
      implicit none
      include 'types.f'
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      
      integer:: n
      character titlex*8
      real(dp)::dx,xmin,xmax
      logical:: first
      data first /.true./
      save first
      
      integer:: NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR

      integer:: NHISTOMAX
      parameter(NHISTOMAX=200)
      integer:: BOOKED(NHISTOMAX)
      common/BOOKPATTERN/BOOKED

      integer:: ISTAT,i
      character*100 outfile

      character*255 runname
      integer:: lenocc
      real(dp):: savedx(200)
      common/runname/runname
      common/dswsavedx/savedx
c     ------------------------------------------------------------------
      
      if (first) then
c ---   Open the file :
       outfile=runname(1:lenocc(runname))//'.rz'
        call hlimit(NWPAWC)
c ---   A record length of 8192 words will produce a warning but is 
c ---   OK for most systems (see HBOOK manual, p.21) : 
        call hropen(30,'HISTOS',outfile,'N',8192,ISTAT)
        if (ISTAT.NE.0) then
          write(6,*) 'ERROR ReadInOut : hropen (HISTOS), ISTAT =',
     +                                                   ISTAT
          stop
        endif
        do i=1,NHISTOMAX
          BOOKED(i)=0
        enddo
        first = .false.
      endif

c --- Book the histogram :
      call hcdir('//HISTOS',' ')
      call hbook1(n,titlex,int(sngl((xmax-xmin)/dx)),
     +            sngl(xmin),sngl(xmax),0.)
      BOOKED(n)=1
      savedx(n)=dx

      return
      end
c

      subroutine dswhfill(n,var,wgt)
      implicit none
      include 'types.f'
      integer:: n
      real(dp)::var,wgt

      integer:: NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR
      real(dp):: savedx(200)
      common/dswsavedx/savedx
c     ------------------------------------------------------------------

c --- Fill the histogram :
      call hcdir('//HISTOS',' ')
      call hf1(n,sngl(var),sngl(wgt/savedx(n)))

      return
      end

c

      subroutine dswhrout
      implicit none
      include 'types.f'
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      

      integer:: NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR

      integer:: NHISTOMAX
      parameter(NHISTOMAX=200)
      integer:: BOOKED(NHISTOMAX)
      common/BOOKPATTERN/BOOKED

c --- Common block to control output
c --- (histograms or ntuple)
      include 'outputflags.f'

      integer:: ICYCLE,i
c     ------------------------------------------------------------------

      if (.NOT.creatent) then
c ---   Read out the booked histograms :
        do i=1,NHISTOMAX
          if (BOOKED(i).EQ.1) then
            call hcdir('//HISTOS',' ')
c           write(6,*) 'Calling hrout'
            call hrout(i,ICYCLE,' ')
            if (ICYCLE.NE.1) then
              write(6,*) 'ERROR dswhrout : ICYCLE =',ICYCLE
              write(6,*) 'ERROR dswhrout : Dupicated HISTO ID'
            endif
          endif
        enddo

      else

c ---   Read out the ntuple :
c ---   =====================
c       write(6,*) 'Calling hrout'
c ---   For some reason, this call to hrout produces an ICYCLE of 2.
c ---   Possibly because the ntuple is explicitly a disk ntuple.
        call hrout(300,ICYCLE,' ')
c       if (ICYCLE.NE.1) then
c         write(6,*) 'ERROR dswhrout : ICYCLE =',ICYCLE
c         write(6,*) 'ERROR dswhrout : Dupicated NTUPLE ID'
c       endif

      endif

      return
      end

c

      subroutine dswclose
      implicit none
      include 'types.f'
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      

      integer:: NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR

      integer:: ICYCLE
c     ------------------------------------------------------------------

      call hrend('HISTOS')
      close(30)
      
      write(6,*) '<----- Completed a batch of n-tuples ----->' 
      call flush(6)
      
      return
      end

c

      subroutine bookfill(tag,p,wt)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'maxwt.f'
      integer tag
      real(dp):: p(mxpart,4)
      real(dp):: wt 
      integer, parameter:: tagbook=1, tagplot=2

      if (.not.skipnt) then
        if (tag == tagbook) then
          call dswntuplebook
        elseif (tag == tagplot) then
          call dswntuplefill(p,wt)
        endif
      endif

      return
      end

c
      subroutine dswntuplebook
      implicit none
      include 'types.f'
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      
c--- Included for the extra code below
      include 'npart.f'
      include 'mxdim.f'
      include 'scale.f'
      include 'facscale.f'
      include 'PDFerrors.f'
      include 'kpart.f'

      real(dp):: scale_store,facscale_store

      integer:: NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR
      integer:: IQUEST
      common/QUEST/IQUEST(100)

c--- Added to keep track of number of momenta entries to be filled
c--- Extra definitions to facilitate dummy call to lowint
      real(dp):: dummy,wgt,r(mxdim),lowint
      integer:: ifill
      integer:: imaxmom,ipdf
      common/iarray/imaxmom,ipdf
            
      include 'chtags.f'
      
      integer:: ISTAT
      character*100 outfile

      character*255 runname
      integer:: lenocc
      common/runname/runname

      logical:: first
      integer:: batchno
      character*3 batchstr,getstr
      data first/.true./
      save first,batchno
c     ------------------------------------------------------------------

c--- Need to ascertain the correct size for momenta n-tuples when this routine
c--- is called for the first time, achieved via a dummy call to lowint
      if (first) then      
        do ifill=1,mxdim
          r(ifill)=0.5_dp
        enddo
c--- Be careful that dynamic scale choices aren't ruined
c--- (in versions 5.1 and before, this occured when calling lowint)
        scale_store=scale
       facscale_store=facscale
        dummy=lowint(r,wgt)
        scale=scale_store
       facscale=facscale_store
       
        imaxmom=npart
        if ((kpart==kreal).or.(kpart==ktota).or.(kpart==ktodk))
     &    imaxmom=imaxmom+1
         batchno=-1
         first=.false.
      endif

c--- determine if we need space in array to store PDF weights (ipdf)
      ipdf=0
      if (PDFerrors) then
        ipdf= maxPDFsets
      endif
           
c---- Increment the batch number counter and convert to a string    
      batchno=batchno+1
      batchstr=getstr(batchno)
      
c --- Create the output file :
      outfile=runname(1:lenocc(runname))//batchstr//'.rz'
      write(6,*) '<----- Creating batch ',batchno,' of n-tuples ----->' 
      call flush(6)

      call hlimit(NWPAWC)
      IQUEST(10) = 65000
c --- A record length of 8192 words will produce a warning but is 
c --- OK for most systems (see HBOOK manual, p.21) : 
      call hropen(30,'HISTOS',outfile,'N',8192,ISTAT)
      if (ISTAT.NE.0) then
        write(6,*) 'ERROR ReadInOut : hropen (HISTOS), ISTAT =',
     +        ISTAT
      endif

c --- Book an extremely simple row-wise ntuple. Make it explicitly
c --- a disk resident ntuple by specifying the top directory name
c --- of the previously opened RZ file in the 4th argument 
c --- (see HBOOK manual, p.19) :
      if     (ipdf == 0) then
c--- don't need extra space for PDF sets
        if     (imaxmom == 2) then
          call hbookn(300,'MCFM',13,'//HISTOS',4096,CHTAGS4)
        elseif (imaxmom == 3) then
          call hbookn(300,'MCFM',17,'//HISTOS',4096,CHTAGS5)
        elseif (imaxmom == 4) then
          call hbookn(300,'MCFM',21,'//HISTOS',4096,CHTAGS6)
        elseif (imaxmom == 5) then
          call hbookn(300,'MCFM',25,'//HISTOS',4096,CHTAGS7)
        elseif (imaxmom == 6) then
          call hbookn(300,'MCFM',29,'//HISTOS',4096,CHTAGS8)
        elseif (imaxmom == 7) then
          call hbookn(300,'MCFM',33,'//HISTOS',4096,CHTAGS9)
        elseif (imaxmom == 8) then
          call hbookn(300,'MCFM',37,'//HISTOS',4096,CHTAGS10)
        else
          write(6,*) 'Problem in dswntuplebook - value npart=',npart
          write(6,*) 'not anticipated. Program halted.'
          stop
        endif
      elseif (ipdf == 30) then
        if     (imaxmom == 2) then
          call hbookn(300,'MCFM',43,'//HISTOS',4096,CHTAGS4p30)
        elseif (imaxmom == 3) then
          call hbookn(300,'MCFM',47,'//HISTOS',4096,CHTAGS5p30)
        elseif (imaxmom == 4) then
          call hbookn(300,'MCFM',51,'//HISTOS',4096,CHTAGS6p30)
        elseif (imaxmom == 5) then
          call hbookn(300,'MCFM',55,'//HISTOS',4096,CHTAGS7p30)
        elseif (imaxmom == 6) then
          call hbookn(300,'MCFM',59,'//HISTOS',4096,CHTAGS8p30)
        elseif (imaxmom == 7) then
          call hbookn(300,'MCFM',63,'//HISTOS',4096,CHTAGS9p30)
        elseif (imaxmom == 8) then
          call hbookn(300,'MCFM',67,'//HISTOS',4096,CHTAGS10p30)
        else
          write(6,*) 'Problem in dswntuplebook - value npart=',npart
          write(6,*) 'not anticipated. Program halted.'
          stop
        endif
      elseif (ipdf == 40) then
        if     (imaxmom == 2) then
          call hbookn(300,'MCFM',53,'//HISTOS',4096,CHTAGS4p40)
        elseif (imaxmom == 3) then
          call hbookn(300,'MCFM',57,'//HISTOS',4096,CHTAGS5p40)
        elseif (imaxmom == 4) then
          call hbookn(300,'MCFM',61,'//HISTOS',4096,CHTAGS6p40)
        elseif (imaxmom == 5) then
          call hbookn(300,'MCFM',65,'//HISTOS',4096,CHTAGS7p40)
        elseif (imaxmom == 6) then
          call hbookn(300,'MCFM',69,'//HISTOS',4096,CHTAGS8p40)
        elseif (imaxmom == 7) then
          call hbookn(300,'MCFM',73,'//HISTOS',4096,CHTAGS9p40)
        elseif (imaxmom == 8) then
          call hbookn(300,'MCFM',77,'//HISTOS',4096,CHTAGS10p40)
        else
          write(6,*) 'Problem in dswntuplebook - value npart=',npart
          write(6,*) 'not anticipated. Program halted.'
          stop
        endif
      elseif (ipdf == 44) then
        if     (imaxmom == 2) then
          call hbookn(300,'MCFM',57,'//HISTOS',4096,CHTAGS4p44)
        elseif (imaxmom == 3) then
          call hbookn(300,'MCFM',61,'//HISTOS',4096,CHTAGS5p44)
        elseif (imaxmom == 4) then
          call hbookn(300,'MCFM',65,'//HISTOS',4096,CHTAGS6p44)
        elseif (imaxmom == 5) then
          call hbookn(300,'MCFM',69,'//HISTOS',4096,CHTAGS7p44)
        elseif (imaxmom == 6) then
          call hbookn(300,'MCFM',73,'//HISTOS',4096,CHTAGS8p44)
        elseif (imaxmom == 7) then
          call hbookn(300,'MCFM',77,'//HISTOS',4096,CHTAGS9p44)
        elseif (imaxmom == 8) then
          call hbookn(300,'MCFM',81,'//HISTOS',4096,CHTAGS10p44)
        else
          write(6,*) 'Problem in dswntuplebook - value npart=',npart
          write(6,*) 'not anticipated. Program halted.'
          stop
        endif
      else
          write(6,*) 'Number of PDF uncertainty sets, ',ipdf
          write(6,*) 'not anticipated. Program halted.'
          stop        
      endif
      
      return
      end

c
      subroutine dswntuplefill(p,wt)
      implicit none
      include 'types.f'
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'wts_bypart.f'
      include 'PDFerrors.f'
      
      real(dp):: p(mxpart,4)
      real(dp):: wt 

      integer:: NWPAWC
      parameter(NWPAWC=10000000)
      real         HMEMOR(NWPAWC)
      common/PAWC/ HMEMOR

c--- Extra common block to carry the information about maximum momenta entries
      integer:: imaxmom,ipdf
      common/iarray/imaxmom,ipdf

      integer:: i
      real pfill(imaxmom*4+5+ipdf)

c--- Variables to count the number of ntuples filled so far and set the
c--- maximum number filled per batch; the current value is somewhat arbitrary
c--- and may be modified by the user
      integer:: icount,batchlimit
      data icount/0/
      save icount
      
      parameter(batchlimit=1000000)

c     ------------------------------------------------------------------

c--- if the weight is zero, don't bother to add the n-tuple
      if (wt == 0._dp) then
        return
      endif

c--- The n-tuple should be filled, but first check if the current batch
c--- is already full; if so, close and then open a new one
      icount=icount+1
      if (mod(icount,batchlimit) == 0) then
        call dswhrout
        call dswclose
        call dswntuplebook
      endif
      
c --- Fill the ntuple :
      do i=1,imaxmom
        pfill(4*(i-1)+1)=sngl(p(i+2,1))
        pfill(4*(i-1)+2)=sngl(p(i+2,2))
        pfill(4*(i-1)+3)=sngl(p(i+2,3))
        pfill(4*(i-1)+4)=sngl(p(i+2,4))
      enddo
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
      
c     write(6,*) 'Filling ntuple with weight',wt

      call hfn(300,pfill)
      
      return
      end

c


c--- finalize the ntuple by producing output and closing file
      subroutine NTfinalize
      implicit none
      include 'types.f'
      
      
      call dswhrout
      call dswclose
      
      return
      end
      
