      subroutine reader_input(inputfile,workdir)
      implicit none
      include 'types.f'
************************************************************************
*   Routine to read in the file input.DAT, which is a consolidated     *
*   form of all the input files and new to version 3.4 of MCFM         *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'debug.f'
      include 'new_pspace.f'
      include 'couple.f'
      include 'kpart.f'
      include 'virtonly.f'
      include 'realonly.f'
      include 'noglue.f'
      include 'realwt.f'
      include 'lc.f'
      include 'cutoff.f'
      include 'maxwt.f'
      include 'masses.f'
      include 'kprocess.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'zerowidth.f'
      include 'removebr.f'
      include 'flags.f'
      include 'clustering.f'
      include 'anomcoup.f'
      include 'gridinfo.f'
      include 'verbose.f'
      include 'limits.f'
      include 'werkdir.f'
      include 'jetcuts.f'
      include 'leptcuts.f'
      include 'lhapdf.f'
      include 'alfacut.f'
      include 'betacut.f'
      include 'pdlabel.f'
      include 'qcdcouple.f'
      include 'nlooprun.f'
      include 'initialscales.f'
      include 'stopscales.f'
      include 'vanillafiles.f'
      include 'frag.f'
      include 'outputoptions.f'
      include 'outputflags.f'
      include 'breit.f'
      include 'anomHiggs.f'
      include 'anom_higgs.f'
      include 'vdecayid.f'
      include 'runstring.f'
      include 'energy.f'
      include 'nproc.f'
      include 'taucut.f'
      include 'iterat.f'
      include 'mpicommon.f'
c--- APPLgrid - flag using grid
c      include 'ptilde.f'
c      include 'APPLinclude.f'
c--- APPLgrid - end
      character*72 workdir,inputfile
      character*90 line
      character*15 part
      logical:: spira,dryrun,makecuts,creategrid,writerefs
      integer:: nmin,nmax,ii
      integer:: ih1,ih2,idum,origij
      integer:: NPTYPE,NGROUP,NSET
      real(dp):: rtsmin,factor
      real(dp):: mbbmin,mbbmax,Mwmin,Mwmax
      real(dp):: Rcut
      logical:: technicalincluded
      real(dp):: ran2,ran2nr,randummy
      real(dp):: alphas
      
      common/writerefs/writerefs
      common/spira/spira
      common/nmin/nmin
      common/nmax/nmax
      common/rtsmin/rtsmin 

      common/density/ih1,ih2
      common/ranno/idum
      common/dryrun/dryrun
      
      common/pdflib/NPTYPE,NGROUP,NSET
      
      common/Rcut/Rcut
      common/makecuts/makecuts

c      common/qmass/cmass,bmass

      common/origij/origij

      save /ranno/

      data cttH/1.0_dp/
      data cWWH/1.0_dp/

      verbose=.false.
      if (rank.eq.0) verbose=.true.

      werkdir=workdir
c--- work out the name of the input file and open it

      if (rank.eq.0) write(6,*) '* Using input file named ',inputfile

      open(unit=20,file=inputfile,status='old',err=999)
      call checkversion(20,inputfile)

c--- rea.e-_dpin the user inputs

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- flags for the mode of MCFM   
      read(20,*) nevtrequested
      if (verbose) call writeinput(6,' * ',' ','nevtrequested')
      read(20,*) creatent
      if (verbose) call writeinput(6,' * ',' ','creatent')
      skipnt=.false.
!      read(20,*) skipnt
!      if (verbose) call writeinput(6,' * ',' ','skipnt')
      read(20,*) dswhisto
      if (verbose) call writeinput(6,' * ',' ','dswhisto')
c--- APPLgrid - added to read flag for grids
       read(20,*) creategrid
       if (verbose) call writeinput(6,' * ',' ','creategrid')
c--- APPLgrid - end
      read(20,*) writerefs
      if (verbose) call writeinput(6,' * ',' ','writerefs')
      read(20,*) writetop
      if (verbose) call writeinput(6,' * ',' ','writetop')
      read(20,*) writedat
      if (verbose) call writeinput(6,' * ',' ','writedat')
      read(20,*) writegnu
      if (verbose) call writeinput(6,' * ',' ','writegnu')
      read(20,*) writeroot
      if (verbose) call writeinput(6,' * ',' ','writeroot')
      read(20,*) writepwg
      if (verbose) call writeinput(6,' * ',' ','writepwg')

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- general options

c--- read in whole line for nproc
      read(20,99) line
      ii=index(line,'.')
      if (ii > 0) then
        vdecayid=.true.
c------ special string present to specify V decays
        read(line(1:ii-1),*) nproc
        v34id=line(ii+1:ii+2)
        v56id=line(ii+3:ii+4)
c        write(6,*) 'special nproc=',nproc
c        write(6,*) 'v34id',v34id
c        write(6,*) 'v56id',v56id
c        stop
c------ normal case
      else
        vdecayid=.false.
        read(line,*) nproc
      endif
      if (verbose) call writeinput(6,' * ',' ','nproc')
      read(20,*) part
      coeffonly=.false.
      kpart=0
      if     ((part == 'lo') .or. (part == 'lord')) then
        kpart=klord
      elseif (part == 'virt') then
        kpart=kvirt
      elseif (part == 'real') then
        kpart=kreal
      elseif ((part == 'nlo') .or. (part == 'tota')
     &   .or. (part == 'nlocoeff') .or. (part == 'totacoeff')) then
        kpart=ktota
      elseif (part == 'frag') then
        kpart=kfrag
      elseif (part == 'todk') then
        kpart=ktodk
      elseif ((part == 'snlo') .or. (part == 'scetnlo')
     &   .or. (part == 'snlocoeff') .or. (part == 'scetnlocoeff')) then
        kpart=ksnlo
      elseif ((part == 'nnlo') .or. (part == 'nnlocoeff')) then
        kpart=knnlo
      endif
      if (index(part,'coeff') > 0) then
        coeffonly=.true.
      endif
      if (kpart == 0) then
        write(6,*) 'Invalid value of part = ',part
        stop
      endif
      if (verbose) call writeinput(6,' * ',' ','part')
      read(20,*) runstring
      if (verbose) call writeinput(6,' * ',' ','runstring')
      read(20,*) sqrts
      if (verbose) call writeinput(6,' * ',' ','sqrts')
      read(20,*) ih1
      if (verbose) call writeinput(6,' * ',' ','ih1')
      read(20,*) ih2
      if (verbose) call writeinput(6,' * ',' ','ih2')
      read(20,*) hmass
      if (verbose) call writeinput(6,' * ',' ','hmass')
      read(20,*) scale
      initscale=scale
      if (verbose) call writeinput(6,' * ',' ','scale')
      read(20,*) facscale
      initfacscale=facscale
      if (verbose) call writeinput(6,' * ',' ','facscale')
      
      initrenscale_L=0._dp
      initfacscale_L=0._dp
      initrenscale_H=0._dp
      initfacscale_H=0._dp
c--- catch special scale choices for stop+b process
      if (((nproc >= 231) .and. (nproc <= 240)) .and.      
     &     (scale == 0._dp) .and. (facscale == 0._dp)) then
        read(20,*) initrenscale_L
      renscale_L=initrenscale_L
        if (verbose) call writeinput(6,' * ',' ','renscale_L')
        read(20,*) initfacscale_L
      facscale_L=initfacscale_L
        if (verbose) call writeinput(6,' * ',' ','facscale_L')
        read(20,*) initrenscale_H
      renscale_H=initrenscale_H
        if (verbose) call writeinput(6,' * ',' ','renscale_H')
        read(20,*) initfacscale_H
      facscale_H=initfacscale_H
        if (verbose) call writeinput(6,' * ',' ','facscale_H')
      scale=initrenscale_H
      facscale=initfacscale_H        
      endif

      read(20,*) dynstring 
      if (verbose) call writeinput(6,' * ',' ','dynamicscale')
      read(20,*) zerowidth
      if (verbose) call writeinput(6,' * ',' ','zerowidth')
      read(20,*) removebr
      if (verbose) call writeinput(6,' * ',' ','removebr')
      read(20,*) itmx1
      if (verbose) call writeinput(6,' * ',' ','itmx1')
      read(20,*) ncall1
      if (verbose) call writeinput(6,' * ',' ','ncall1')
      read(20,*) itmx2
      if (verbose) call writeinput(6,' * ',' ','itmx2')
      read(20,*) ncall2
      if (verbose) call writeinput(6,' * ',' ','ncall2')
      read(20,*) line
      if     (index(line,'1%acc') > 0) then
        if ((kpart == knnlo) .or. (kpart == ksnlo)) then
          call setuptau(10)
        else
          taucut=0.1_dp
        endif
      elseif (index(line,'0.2%acc') > 0) then
        if ((kpart == knnlo) .or. (kpart == ksnlo)) then
          call setuptau(2)
        else
          taucut=0.1_dp
        endif
      else
        read(line,*) taucut
      endif
      if (verbose) call writeinput(6,' * ',' ','taucut')
      read(20,*) origij
      if (verbose) call writeinput(6,' * ',' ','ij')
      read(20,*) dryrun
      if (verbose) call writeinput(6,' * ',' ','dryrun')
      read(20,*) Qflag
      if (verbose) call writeinput(6,' * ',' ','Qflag')
      read(20,*) Gflag
      if (verbose) call writeinput(6,' * ',' ','Gflag')
      
      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- heavy quark masses 
      read(20,*) mt
      if (verbose) call writeinput(6,' * ',' ','top mass')
      read(20,*) mb
      if (verbose) call writeinput(6,' * ',' ','bottom mass')
      read(20,*) mc
      if (verbose) call writeinput(6,' * ',' ','charm mass')

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- pdf options 
      read(20,*) pdlabel
      if (verbose) call writeinput(6,' * ',' ','pdlabel')
      read(20,*) NGROUP
      if (verbose) call writeinput(6,' * ',' ','NGROUP')
      read(20,*) NSET
      if (verbose) call writeinput(6,' * ',' ','NSET')
      read(20,*) PDFname
      if (verbose) call writeinput(6,' * ',' ','LHAPDF group')
      read(20,*) PDFmember
      if (verbose) call writeinput(6,' * ',' ','LHAPDF set')

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- jets and cuts options 
      read(20,*) Mwmin
      wsqmin=Mwmin**2
      if (verbose) call writeinput(6,' * ',' ','m34min')
      read(20,*) Mwmax 
      if (Mwmax > sqrts*0.9999_dp) Mwmax=sqrts*0.9999_dp ! physical cap on m34max
      wsqmax=Mwmax**2
      if (verbose) call writeinput(6,' * ',' ','m34max')
      read(20,*) mbbmin
      bbsqmin=mbbmin**2
      if (verbose) call writeinput(6,' * ',' ','m56min')
      read(20,*) mbbmax 
      if (mbbmax > sqrts*0.9999_dp) Mbbmax=sqrts*09.999_dp ! physical cap on m56max
      bbsqmax=mbbmax**2
      if (verbose) call writeinput(6,' * ',' ','m56max')
      read(20,*) m3456min 
      if (verbose) call writeinput(6,' * ',' ','m3456min')
      read(20,*) m3456max 
      if (m3456max > sqrts) m3456max=sqrts ! physical cap on m3456max
      if (verbose) call writeinput(6,' * ',' ','m3456max')
      read(20,*) inclusive
      if (verbose) call writeinput(6,' * ',' ','inclusive')
      read(20,*) algorithm
      if (verbose) call writeinput(6,' * ',' ','algorithm')
      read(20,*) ptjetmin
      if (verbose) call writeinput(6,' * ',' ','ptjetmin')
      read(20,*) etajetmin 
      if (verbose) call writeinput(6,' * ',' ','etajetmin')
      read(20,*) etajetmax 
      if (verbose) call writeinput(6,' * ',' ','etajetmax')
      read(20,*) Rcut
      if (verbose) call writeinput(6,' * ',' ','Rcut')
      read(20,*) makecuts
      if (verbose) call writeinput(6,' * ',' ','makecuts')
      read(20,*) leptpt
      if (verbose) call writeinput(6,' * ',' ','leptpt')
      read(20,*) leptrap
      if (verbose) call writeinput(6,' * ',' ','leptrap')
      read(20,*) leptveto1min,leptveto1max
      if (verbose) call writeinput(6,' * ',' ','leptveto')
      read(20,*) misspt
      if (verbose) call writeinput(6,' * ',' ','misspt')
      read(20,*) leptpt2
      if (verbose) call writeinput(6,' * ',' ','leptpt2')
      read(20,*) leptrap2
      if (verbose) call writeinput(6,' * ',' ','leptrap2')
      read(20,*) leptveto2min,leptveto2max
      if (verbose) call writeinput(6,' * ',' ','leptveto2')
      read(20,*) mtrans34cut
      if (verbose) call writeinput(6,' * ',' ','mtrans34cut')
      read(20,*) Rjlmin
      if (verbose) call writeinput(6,' * ',' ','Rjlmin')
      read(20,*) Rllmin
      if (verbose) call writeinput(6,' * ',' ','Rllmin')
      read(20,*) delyjjmin
      if (verbose) call writeinput(6,' * ',' ','delyjjmin')
      read(20,*) jetsopphem 
      if (verbose) call writeinput(6,' * ',' ','jetsopphem')
      read(20,*) lbjscheme 
      if (verbose) call writeinput(6,' * ',' ','lbjscheme')
      read(20,*) ptbjetmin
      if (verbose) call writeinput(6,' * ',' ','ptbjetmin')
      read(20,*) etabjetmax
      if (verbose) call writeinput(6,' * ',' ','etabjetmax')

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- settings for photon processes 
      read(20,*) frag
      if (verbose) call writeinput(6,' * ',' ',kfrag)
      read(20,*) fragset
      if (verbose) call writeinput(6,' * ',' ','fragset')
      read(20,*) frag_scale
      frag_scalestart=frag_scale
      if (verbose) call writeinput(6,' * ',' ','frag_scale')
      read(20,*) gammpt
      if (verbose) call writeinput(6,' * ',' ','gammpt')
      read(20,*) gammrap
      if (verbose) call writeinput(6,' * ',' ','gammrap')
      read(20,*) gammpt2
      if (verbose) call writeinput(6,' * ',' ','gammpt2')
      read(20,*) gammpt3
      if (verbose) call writeinput(6,' * ',' ','gammpt3')
      read(20,*) Rgalmin 
      if (verbose) call writeinput(6,' * ',' ','Rgalmin')
      read(20,*) Rgagamin 
      if (verbose) call writeinput(6,' * ',' ','Rgagamin')
      read(20,*) Rgajetmin
      if (verbose) call writeinput(6,' * ',' ','Rgajetmin')
      read(20,*) cone_ang
      if (verbose) call writeinput(6,' * ',' ','cone_ang')
      read(20,*) epsilon_h
      if (verbose) call writeinput(6,' * ',' ','epsilon_h')
      read(20,*) n_pow
      if (verbose) call writeinput(6,' * ',' ','n_pow')

      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- anomalous couplings 
      read(20,*) delg1_z
      if (verbose) call writeinput(6,' * ',' ','delg1_z')
      read(20,*) delk_z
      if (verbose) call writeinput(6,' * ',' ','delk_z')
      read(20,*) delk_g
      if (verbose) call writeinput(6,' * ',' ','delk_g')
      read(20,*) lambda_z
      if (verbose) call writeinput(6,' * ',' ','lambda_z')
      read(20,*) lambda_g
      if (verbose) call writeinput(6,' * ',' ','lambda_g')
      read(20,*) h1Z
      if (verbose) call writeinput(6,' * ',' ','h1Z')
      read(20,*) h1gam
      if (verbose) call writeinput(6,' * ',' ','h1gam')
      read(20,*) h2Z
      if (verbose) call writeinput(6,' * ',' ','h2Z')
      read(20,*) h2gam
      if (verbose) call writeinput(6,' * ',' ','h2gam')
      read(20,*) h3Z
      if (verbose) call writeinput(6,' * ',' ','h3Z')
      read(20,*) h3gam
      if (verbose) call writeinput(6,' * ',' ','h3gam')
      read(20,*) h4Z
      if (verbose) call writeinput(6,' * ',' ','h4Z')
      read(20,*) h4gam
      if (verbose) call writeinput(6,' * ',' ','h4gam')
      read(20,*) tevscale
      if (verbose) call writeinput(6,' * ',' ','tevscale')

      if ( nproc >= 550 .and. nproc <= 557 ) then
         if (verbose) call writeinput(6,' * ',' ','cttH')
         if (verbose) call writeinput(6,' * ',' ','cWWH')
      endif
      
      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- width of the Higgs 
      read(20,*) hwidth_ratio
      if (verbose) call writeinput(6,' * ',' ','hwidth_ratio')
      
      if (verbose) write(6,*)
      read(20,99) line
c--- write-out comment line
      read(20,99) line
      if (verbose) write(6,*) '* ',line
c--- grid information 
      read(20,*) readin
      if (verbose) call writeinput(6,' * ',' ','readin')
      read(20,*) writeout
      if (verbose) call writeinput(6,' * ',' ','writeout')
      read(20,*) ingridfile
      if (verbose) call writeinput(6,' * ',' ','ingridfile')
      read(20,*) outgridfile
      if (verbose) call writeinput(6,' * ',' ','outgridfile')

      if (verbose) write(6,*)

c--- check if the contents of technical.DAT are included here
c--- (the default behaviour going forward)
      technicalincluded=.false.

      read(20,99,end=88) line
      read(20,99,end=88) line
c      write(6,*) line
      if (line(2:10) .ne. 'Technical') goto 88
      technicalincluded=.true.

   88 continue
      
c--- if the contents of technical.DAT are not included, close the input
c--- file and open technical.DAT instead; otherwise continue on
      if (technicalincluded .eqv. .false.) then
        close(20)
        open(unit=20,file='technical.DAT',status='old',err=999)
        call checkversion(20,'technical.DAT')
      endif

      if (verbose) write(6,*) '* [Technical parameters that'//
     &                        ' should not normally be changed]'
      if (verbose) write(6,*)

c---- read in the technical parameters

      read(20,*) debug
      if (verbose) call writeinput(6,' * ',' ','debug')
      read(20,*) verbose
      if (rank.ne.0) verbose=.false.
      if (verbose) call writeinput(6,' * ',' ','verbose')
      read(20,*) new_pspace
      if (verbose) call writeinput(6,' * ',' ','new_pspace')
      read(20,*) virtonly
      if (verbose) call writeinput(6,' * ',' ','virtonly')
      read(20,*) realonly
      if (verbose) call writeinput(6,' * ',' ','realonly')
      read(20,*) spira
      if (verbose) call writeinput(6,' * ',' ','spira')
      read(20,*) noglue
      if (verbose) call writeinput(6,' * ',' ','noglue')
      read(20,*) ggonly
      if (verbose) call writeinput(6,' * ',' ','ggonly')
      read(20,*) gqonly
      if (verbose) call writeinput(6,' * ',' ','gqonly')
      read(20,*) omitgg
      if (verbose) call writeinput(6,' * ',' ','omitgg')
      read(20,*) vanillafiles
      if (verbose) call writeinput(6,' * ',' ','vanillafiles')
! 4/28/16: removed nmin, nmax from input file (no longer used)
!          set values to nmin=1, nmax=2 in case of unforeseen problems
      nmin=1
      nmax=2
!      read(20,*) nmin
!      if (verbose) call writeinput(6,' * ',' ','nmin')
!      read(20,*) nmax
!      if (verbose) call writeinput(6,' * ',' ','nmax')
      read(20,*) clustering
      if (verbose) call writeinput(6,' * ',' ','clustering')
      read(20,*) realwt
      if (verbose) call writeinput(6,' * ',' ','realwt')
      read(20,*) colourchoice
      if (verbose) call writeinput(6,' * ',' ','colourchoice')
      read(20,*) rtsmin
      if (verbose) call writeinput(6,' * ',' ','rtsmin')
      read(20,*) cutoff
      if (verbose) call writeinput(6,' * ',' ','cutoff')
      read(20,*) aii
      if (verbose) call writeinput(6,' * ',' ','aii')
      read(20,*) aif
      if (verbose) call writeinput(6,' * ',' ','aif')
      read(20,*) afi
      if (verbose) call writeinput(6,' * ',' ','afi')
      read(20,*) aff
      if (verbose) call writeinput(6,' * ',' ','aff')
      read(20,*) bfi
      if (verbose) call writeinput(6,' * ',' ','bfi')
      read(20,*) bff
      if (verbose) call writeinput(6,' * ',' ','bff')
      if (verbose) write(6,*)
      close(unit=20)

      if ((etajetmin < 0._dp) .or. (etajetmax < 0._dp)) then
        write(6,*) 'etajetmin and etajetmax are absolute values,'
      write(6,*) ' please reset to a positive value.'
      stop
      endif

c--- for W+2 jet and Z+2 jet processes, set aff equal to afi
!      if ( ((nproc == 22) .or. (nproc == 27)
!     &  .or.(nproc == 44) .or. (nproc == 46)) .or.
!     &     ((kpart == knnlo) .and. ((nproc == 11)
!     &         .or. (nproc == 16) .or. (nproc == 41)
!     &         .or. (nproc == 42) .or. (nproc == 43))) ) then
!        aff=afi
!        if (rank == 0) then
!          write(6,*) '>>> Over-riding input value of aff; aff = afi = ',afi
!          write(6,*)
!        endif
!      endif

c--- determine whether SCET is to be used for calculation
      if ((kpart==knnlo) .or. (kpart==ksnlo)) then
        usescet=.true.
        abovecut=.false.
        if (taucut < 0) then
          write(6,*) 'Must specify taucut > 0 for SCET calculation'
          stop
        endif
      else
        usescet=.false.
        abovecut=.false.
      endif
      
      if     (runstring(1:4) == 'tau0') then
        usescet=.true.
        abovecut=.true.
        ntau=0
        if (rank == 0) then
          write(6,*) 'WARNING: performing SCET calculation with tau0'
        endif
      elseif (runstring(1:4) == 'tau1') then
        usescet=.true.
        abovecut=.true.
        ntau=1
        if (rank == 0) then
          write(6,*) 'WARNING: performing SCET calculation with tau1'
        endif
      endif

      if (index(runstring,'toponly') > 0) then
        toponly=.true.
      else
        toponly=.false.
      endif
      
c      if     (index(runstring,'mc1.3') > 0) then
c        mc=1.3_dp
c      mcsq=mc**2
c      elseif (index(runstring,'mc1.4') > 0) then
c        mc=1.4_dp
c      mcsq=mc**2
c      elseif (index(runstring,'mc1.5') > 0) then
c        mc=1.5_dp
c      mcsq=mc**2
c      endif
      
c      if (runstring(1:3) == 'mlm') then
c        write(6,*) 'WARNING: cross sections divided by Ecm**2'
c      write(6,*)
c      endif
      
c---  create logical:: variable dynamicscale for use in other routines
      if (  (dynstring == 'no') .or. (dynstring == '.false.')
     & .or. (dynstring == 'none') ) then 
         dynamicscale=.false. 
      else
         dynamicscale=.true. 
      endif

c--- print warning messages if some parton fluxes are not included      
      if (noglue) then
        write(6,*) 'WARNING: no gluon contribution included in PDF'
      write(6,*)
      endif
      if (ggonly) then
        write(6,*) 'WARNING: only gluon-gluon flux included'
      write(6,*)
      endif
      if (gqonly) then
        write(6,*) 'WARNING: only gluon-quark flux included'
      write(6,*)
      endif
      if (omitgg) then
        write(6,*) 'WARNING: no gluon-gluon contribution included'
      write(6,*)
      endif
      
c--- assign squared masses for b- and c-quarks
      if (abs(mb) > 1.e-8_dp) then
        mbsq=mb**2
      else
        mbsq=4.75_dp**2
      endif
      if (abs(mc) > 1.e-8_dp) then
        mcsq=mc**2
      else
        mcsq=1.5_dp**2
      endif
      
c--- set-up the variables for the process we wish to consider
      call chooser

c--- set-up the random number generator with a negative seed
      idum=-abs(origij)
      randummy=ran2nr()

c--- initialize masses for alpha_s routine ! this is now done in coupling2.f
c      cmass=sqrt(mcsq)
c      bmass=sqrt(mbsq)

c--- E-M gauge invariance requires that delg1_g=0
      delg1_g=0._dp

c--- check that we have a valid value of 'part'
      if ( (kpart.ne.klord) .and. (kpart.ne.kreal) .and.
     &     (kpart.ne.kvirt) .and. (kpart.ne.ktota) ) then
        if    ( (kpart==ktodk) .and.
     &          ((kcase==kbq_tpq) .or. (kcase==kt_bbar)
     &      .or. (kcase==kW_twdk) .or. (kcase==ktt_bbl)
     &      .or. (kcase==ktt_bbh) .or. (kcase==k4ftwdk)
     &      .or. (kcase==kHWW2lq) .or. (kcase==kqq_ttw)
     &      .or. (kcase==kWWqqbr) .or. (kcase==ktt_bbu)
     &      .or. (kcase==kWHbbar) .or. (kcase==kZHbbar)) ) then
c--- this is an allowed combination
        elseif ( (kpart==kfrag) .and.
     &          ((kcase==kWgamma) .or. (kcase==kZgamma)
     &      .or. (kcase==kgamgam) .or. (kcase==kgg2gam)
     &      .or. (kcase==kdirgam) .or. (kcase==kdm_gam)
     &      .or. (kcase==kgmgmjt) .or .(kcase==ktrigam)
     &      .or. (kcase==kfourga)
     &      .or. (kcase==kZ_2gam) .or. (kcase==kZgajet)
     &      .or. (kcase==kW_2gam)) ) then
c--- this is an allowed combination
        elseif ((kpart==knnlo) .or. (kpart==ksnlo)) then 
c--- this is an allowed combination
        else 
          write(6,*) 'part=',part,' is not a valid option'
          write(6,*) 'for this process number.'
          stop     
        endif
      endif      

c--- check that we are not trying to calculate radiation in decay at LO
      if    ( (kpart==klord) .and.
     &       ((kcase==kWWqqdk) .or. (kcase==kHWWdkW)
     &   .or. (kcase==ktt_ldk) .or. (kcase==ktt_udk)
     &   .or. (kcase==ktt_hdk) .or. (kcase==ktthWdk)
     &   .or. (kcase==kttdkay) .or. (kcase==kWtdkay)
     &   .or. (kcase==kdk_4ft) .or. (kcase==kttwldk)) ) then
          write(6,*) 'This process number cannot be used for'
          write(6,*) 'a LO calculation.'
          stop     
      endif

c--- Assign choice of jet algorithm to integer variable
      if     (algorithm == 'ktal') then
        jetalgorithm=kt
      elseif (algorithm == 'ankt') then
        jetalgorithm=antikt
      elseif (algorithm == 'cone') then
        jetalgorithm=Rsepcone
      elseif (algorithm == 'hqrk') then
        jetalgorithm=hqrk
      elseif (algorithm == 'none') then
        jetalgorithm=noclustering
      else
        write(6,*) 'Invalid choice of jet algorithm: should be one of'
        write(6,*) 'ktal, ankt, cone, hqrk, none'
        stop
      endif

c--- set up the default choices of static scale, if required
      if (scale < 0._dp) then
      if     (scale == -2._dp) then
        factor=0.25_dp
      elseif (scale == -3._dp) then
        factor=0.5_dp
      elseif (scale == -4._dp) then
        factor=0.75_dp
      elseif (scale == -5._dp) then
        factor=1._dp
      elseif (scale == -6._dp) then
        factor=2._dp
      elseif (scale == -7._dp) then
        factor=4._dp
        else
        factor=1._dp
      endif        
        if ((n2+n3 .ne. 0) .or. (kcase==ktt_tot)) then
c--- special case for t-tbar production
        if (kcase==ktt_tot) then
          scale=factor*mt
c--- special cases where Higgs mass is neither mass2 nor mass3
c        elseif ((case(1:1) == 'H') .or. (kcase==kWHbbar)
c     &     .or. (kcase==kZHbbar) .or. (kcase==kqq_Hqq)) then
c          scale=factor*hmass
        else     
          scale=factor*(n2*mass2+n3*mass3)/real(n2+n3,dp)
        endif
        as=alphas(scale,amz,nlooprun)
        ason2pi=as/twopi
        ason4pi=as/fourpi
        gsq=fourpi*as
        musq=scale**2
        write(6,*)
        write(6,*)'************* Strong coupling, alpha_s  ************'
        write(6,*)'*                                                  *'
        write(6,49)'alpha_s (scale)',gsq/fourpi
        write(6,49)'alpha_s (zmass)',amz
        write(6,50)' (using ',nlooprun,'-loop running of alpha_s)'  
        write(6,*)'****************************************************'
        write(6,*)
        write(6,*)'****************************************************'
        write(6,76) scale
        write(6,*)'****************************************************'
        else
        write(6,*) 'Invalid choice of renormalization scale!'
        stop
        endif
      endif
      if (facscale < 0._dp) then
      if     (facscale == -2._dp) then
        factor=0.25_dp
      elseif (facscale == -3._dp) then
        factor=0.5_dp
      elseif (facscale == -4._dp) then
        factor=0.75_dp
      elseif (facscale == -5._dp) then
        factor=1._dp
      elseif (facscale == -6._dp) then
        factor=2._dp
      elseif (facscale == -7._dp) then
        factor=4._dp
        else
        factor=1._dp
      endif        
        if ((n2+n3 .ne. 0) .or. (kcase==ktt_tot)) then
c--- special case for t-tbar production
        if (kcase==ktt_tot) then
          facscale=factor*mt
c--- special cases where Higgs mass is neither mass2 nor mass3
c        elseif ((case(1:1) == 'H') .or. (kcase==kWHbbar)
c     &     .or. (kcase==kZHbbar) .or. (kcase==kqq_Hqq)) then
c          facscale=factor*hmass
        else     
          facscale=factor*(n2*mass2+n3*mass3)/real(n2+n3,dp)
        endif
        write(6,*)
        write(6,*)'****************************************************'
        write(6,77) facscale
        write(6,*)'****************************************************'
       else
        write(6,*) 'Invalid choice of factorization scale!'
        stop
        endif
      endif
      
      return

   49 format(' *  ',a20,f12.8,16x,'*')
   50 format(' *  ',6x,a8,i1,a25,8x,'*')
   76 format(' *      Renormalization scale =',f7.2,'              *')
   77 format(' *        Factorization scale =',f7.2,'              *')
   99 format(a90)

  999 continue
      write(6,*) 'Problem reading ',inputfile
      write(6,*)
      write(6,*) 'Refer to documentation for the format of input.DAT'
      write(6,*)
      stop

      end
      
