*****************
* PDFLIB version*
*****************
      subroutine pdfwrap
      implicit none
      include 'masses.f'
      include 'pdlabel.f'
      include 'couple.f'
* PDFLIB variables
      integer NPTYPE,NGROUP,NSET
      character*20 pdflib_parm(20)
      double precision pdflib_val(20),alphas2
      common/pdflib/NPTYPE,NGROUP,NSET
      

      pdflib_parm(1)='INIT0'
      pdflib_val(1)=0d0
      call pdfset(pdflib_parm,pdflib_val)
      
      NPTYPE=1
      
c      if (newinput .eqv. .false.) then
c        open(unit=21,file='pdflib.DAT',status='old',err=999)
c        call checkversion(21,'pdflib.DAT')
c        read(21,*) NGROUP
c        read(21,*) NSET            
c        close(21)
c      endif
      
      write(6,*)
      write(6,*) '**************************'
      write(6,*) '* MCFM is calling PDFLIB *'
      write(6,*) '*                        *'
      write(6,99) 'NPTYPE',NPTYPE
      write(6,99) 'NGROUP',NGROUP
      write(6,99) 'NSET  ',NSET
      write(6,*) '**************************'

      pdflib_parm(1)='NPTYPE'
      pdflib_parm(2)='NGROUP'
      pdflib_parm(3)='NSET'
      pdflib_val(1)=NPTYPE
      pdflib_val(2)=NGROUP
      pdflib_val(3)=NSET
      call pdfset(pdflib_parm,pdflib_val)
      amz=alphas2(zmass)

c--- rename pdlabel to get sensible output name
      pdlabel='pdflib_'

      return
 
   99 format(' *     ',a8,i3,'        *')

  999 write(6,*) 'Error reading pdflib.DAT'
      call flush(6)
      stop

      end
 

