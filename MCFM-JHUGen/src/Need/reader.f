      block data reader_data
      include 'impsample.f'
      logical:: msbar
      common/msbar/msbar
      data impsample/.false./
      data msbar/.false./
      end

      subroutine reader
      implicit none
      include 'types.f'
      include 'masses.f'
      include 'scale.f'
      include 'debug.f'
      include 'realonly.f'
      include 'virtonly.f'
      include 'noglue.f'
      include 'realwt.f'
      include 'zerowidth.f'
      include 'removebr.f'
      include 'new_pspace.f'
      include 'impsample.f'
      include 'cutoff.f'
      include 'clustering.f'
      include 'flags.f'
      include 'lc.f'
      include 'gridinfo.f'
      include 'maxwt.f'
      include 'verbose.f'
      include 'kprocess.f'
      include 'limits.f'
      include 'pdlabel.f'
      include 'kpart.f'
      include 'outputflags.f'
      include 'runstring.f'
      include 'energy.f'
      include 'nproc.f'
      include 'iterat.f'
      integer:: ih1,ih2,idum,nmin,nmax
      integer:: nargs
      real(dp):: Rcut,randummy,ran2
      real(dp):: cmass,bmass
      character*72 optionsfile
      logical:: makecuts,dryrun
      logical:: spira
      real(dp):: rtsmin
      real(dp):: mbbmin,mbbmax,Mwmin,Mwmax
      common/spira/spira
      common/ranno/idum
      common/Rcut/Rcut


      common/makecuts/makecuts
      common/dryrun/dryrun
      common/qmass/cmass,bmass

      common/density/ih1,ih2

      common/rtsmin/rtsmin
      common/nmin/nmin
      common/nmax/nmax
      
 
      nargs=iargc()
      if (nargs >= 1) then
      call getarg(1,optionsfile)
      else
      optionsfile='options.DAT'
      endif

      write(6,*) 'Using options file named ',optionsfile
      write(6,*) '****************'
      write(6,*) 'Options file:'
      write(6,*)
      open(unit=20,file=optionsfile,status='old',err=999)
      call checkversion(20,optionsfile)
      read(20,*) nproc
      write(6,*) 'nproc=',nproc
      read(20,*) part
      write(6,*) 'part=',part
      read(20,*) runstring
      write(6,*) 'runstring=',runstring
      read(20,*) verbose
      write(6,*) 'verbose=',verbose
      read(20,*) sqrts
      if (verbose) write(6,*) 'sqrts=',sqrts
      read(20,*) ih1
      if (verbose) write(6,*) 'ih1',ih1
      read(20,*) ih2
      if (verbose) write(6,*) 'ih2',ih2
      read(20,*) pdlabel
      if (verbose) write(6,*) 'pdlabel ',pdlabel
      read(20,*) hmass
      if (verbose) write(6,*) 'hmass',hmass
      read(20,*) scale
      if (verbose) write(6,*) 'scale',scale
        read(20,*) Mwmin
      if (verbose) write(6,*) 'm34min',Mwmin
      read(20,*) Mwmax 
      if (verbose) write(6,*) 'm34max',Mwmax
      read(20,*) mbbmin
      if (verbose) write(6,*) 'm56min',mbbmin
      read(20,*) mbbmax 
      if (verbose) write(6,*) 'm56max',mbbmax
      read(20,*) rtsmin
      if (verbose) write(6,*) 'rtsmin',rtsmin
      read(20,*) zerowidth
      if (verbose) write(6,*) 'zerowidth',zerowidth
      read(20,*) makecuts
      if (verbose) write(6,*) 'makecuts',makecuts
      read(20,*) Rcut
      if (verbose) write(6,*) 'Rcut',Rcut
      read(20,*) itmx1
      if (verbose) write(6,*) 'itmx1',itmx1
      read(20,*) ncall1
      if (verbose) write(6,*) 'ncall1',ncall1
      read(20,*) itmx2
      if (verbose) write(6,*) 'itmx2',itmx2
      read(20,*) ncall2
      if (verbose) write(6,*) 'ncall2',ncall2
      read(20,*) idum
      if (verbose) write(6,*) 'idum',idum
      read(20,*) realwt
      if (verbose) write(6,*) 'realwt',realwt
      read(20,*) cutoff
      if (verbose) write(6,*) 'cutoff',cutoff
      read(20,*) dryrun
      if (verbose) write(6,*) 'dryrun',dryrun
      read(20,*) debug
      if (verbose) write(6,*) 'debug',debug
      read(20,*) Qflag
      if (verbose) write(6,*) 'Qflag',Qflag
      read(20,*) Gflag
      if (verbose) write(6,*) 'Gflag',Gflag
      read(20,*) colourchoice
      if (verbose) write(6,*) 'colourchoice',colourchoice

      close(unit=20)

c      if (debug) then
      open(unit=21,file='debug.DAT',status='old',err=999)
      call checkversion(21,'debug.DAT')
      read(21,*) new_pspace
      if (verbose) write(6,*) 'new_pspace',new_pspace
      read(21,*) virtonly
      if (verbose) write(6,*) 'virtonly',virtonly
      read(21,*) realonly
      if (verbose) write(6,*) 'realonly',realonly
      read(21,*) spira
      if (verbose) write(6,*) 'spira',spira
      read(21,*) noglue
      if (verbose) write(6,*) 'noglue',noglue
      read(21,*) ggonly
      if (verbose) write(6,*) 'ggonly',ggonly
      read(21,*) gqonly
      if (verbose) write(6,*) 'gqonly',gqonly
      read(21,*) nmin
      if (verbose) write(6,*) 'nmin',nmin
      read(21,*) nmax
      if (verbose) write(6,*) 'nmax',nmax
      read(21,*) clustering
      if (verbose) write(6,*) 'clustering',clustering
      close(unit=21)
c      endif
      write(6,*) '****************'
      write(6,*)
 
      if (ggonly .and. gqonly) then
        write(6,*) 'ggonly and gqonly BOTH .true. - pick one!'
        stop
      endif

      call chooser
      
      write(6,*)
      write(6,*) '****************'
      write(6,*) 'Grid file:'
      write(6,*)
c--- read in grid options file
      open(unit=21,file='gridinfo.DAT',status='old',err=998)
      call checkversion(21,'gridinfo.DAT')
      read(21,*) readin
      if (verbose) write(6,*) 'readin',readin
      read(21,*) writeout
      if (verbose) write(6,*) 'writeout',writeout
      read(21,99) ingridfile
      read(21,99) outgridfile

      if (ingridfile == '') then
        ingridfile=case//'_'//part//'_grid'
      endif
      if (outgridfile == '') then
        outgridfile=case//'_'//part//'_grid'
      endif

      if (verbose) write(6,*) 'ingridfile =',ingridfile
      if (verbose) write(6,*) 'outgridfile=',outgridfile
      close(unit=21)
      write(6,*) '****************'

* Fill common block 'outputflags' to control output (ntuples or not,
*  DSW histograms or not) - this should probably be in options.DAT
*      creatent = .false.
* For running MCFM standalone, don't want to skip ntupling :
*      skipnt   = .false.
* This determines whether to use hbook or mbook histograms
* independently of the above option :
*      dswhisto= .false.
* this flag determines if we want to generate weighted events
* (running mcfm in old mode) or the user wants to generate 
* (unweighted) events 
*      evtgen = .false.

c--- read in MCFM mode file
      open(unit=21,file='mcfmmode.DAT',status='old',err=997)
      call checkversion(21,'mcfmmode.DAT')
      read(21,*) evtgen
      read(21,*) creatent
      read(21,*) skipnt
      read(21,*) dswhisto
      close(unit=21)

      write(6,*)
      if ((evtgen .eqv. .false.) .and. (creatent .eqv. .false.)) then
        write(6,*) '****************************************'
        write(6,*) '*   MCFM running in traditional mode   *'
        if (dswhisto .eqv. .false.) then
        write(6,*) '*    (traditional MBOOK histograms)    *'
        else
        write(6,*) '*           (HBOOK histograms)         *'
        endif
        write(6,*) '****************************************'
      else
        write(6,*) '****************'
        write(6,*) 'MCFM mode file:'
        write(6,*)
        if (verbose) write(6,*) 'evtgen',evtgen
        if (verbose) write(6,*) 'creatent',creatent
        if (verbose) write(6,*) 'skipnt',skipnt
        if (verbose) write(6,*) 'dswhisto',dswhisto
        write(6,*) '****************'
      endif
c-----initialize various quantities

c--- set-up the random number generator with a negative seed
      idum=-abs(idum)
      randummy=ran2()

c---initialize masses for alpha_s routine
      cmass=sqrt(mcsq)
      bmass=sqrt(mbsq)


      bbsqmin=mbbmin**2
      bbsqmax=mbbmax**2

      wsqmin=Mwmin**2
      wsqmax=Mwmax**2

c

c-----stange-marciano formula for resolution
c      deltam=sqrt(0.64_dp*hmass+0.03_dp**2*hmass**2)
c      if (verbose) write(6,*) 'delta m',deltam

  99  format(a16)

      return

 997  continue
      write(6,*) 'Problem reading mcfmmode.DAT'
      write(6,*)
      write(6,*) 'Required format is:'
      write(6,*) 'logical::     [evtgen]'
      write(6,*) 'logical::     [creatent]'
      write(6,*) 'logical::     [skipnt]'
      write(6,*) 'logical::     [dswhisto]'
      write(6,*)
      stop

 998  continue
      write(6,*) 'Problem reading gridinfo.DAT'
      write(6,*)
      write(6,*) 'Required format is:'
      write(6,*) 'logical::     [readin]'
      write(6,*) 'logical::     [writeout]'
      write(6,*) 'char*16     [ingridfile]'
      write(6,*) 'char*16     [outgridfile]'
      write(6,*)
      write(6,*) 'READIN/WRITEOUT = True/False specify whether'
      write(6,*) 'a grid should be rea.e-_dpin and/or written-out'
      write(6,*) 'INGRIDFILE/OUTGRIDFILE specify the names of the'
      write(6,*) 'files read/written, but may be left blank for default'
      write(6,*)
      stop

 999  continue
      write(6,*) 'Problem reading ',optionsfile
      write(6,*)
      write(6,*) 'Refer to documentation for the format of options.DAT'
      write(6,*)
      stop

      end

