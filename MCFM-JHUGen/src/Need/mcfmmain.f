      subroutine mcfmmain(inputfile,workdir,r,er)
************************************************************************
*                                                                      *
*  This is the main program for MCFM                                   *
*                                                                      *
*  The sequence of calls should always be:                             *
*   call mcfm_init          : basic variable initialization, print-out *
*   call mcfm_vegas(warmup) : warm-up the Vegas grid                   *
*   call mcfm_vegas(accum)  : accumulate results                       *
*   call mcfm_exit          : final processing and print-out           *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'eventbuffer.f'
      include 'gridinfo.f'
      include 'part.f'
      include 'xs_store_info.f'
      include 'vegas_common.f'
      include 'ipsgen.f'
      integer itmx1,ncall1,itmx2,ncall2,itmxplots
      double precision integ,integ_err,r,er
      logical dryrun
      integer i,pflav,pbarflav
      double precision p(mxpart,4),wt
      character*72 inputfile,workdir
      logical newreadimpl
      parameter (newreadimpl = .true.)
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/dryrun/dryrun
       
* basic variable initialization, print-out
      call mcfm_init(inputfile,workdir)

* tell VEGAS to write out pertinent information
      nprn=0
      
* in initial phases, we don't want any unweighting to take place.
* this will be set true in the first call to getevent.
      unweight = .false.
* if we're reading in a grid, there's no need to do any warming-up
      if (readin) dryrun=.true.

* This is the mcfm_vegas(warmup) call
* The Vegas parameters are those read from options.DAT for
* the warm-up stage (itmx1,ncall1) and binning should only take
* place if dryrun is set to true

* There are now 3 modes of operation:
*   dryrun = .false. : warmup, then freeze grid and accumulate
*   dryrun = .true. , readin = .false. : accumulate during warmup
*   dryrun = .true. , readin = .true.  : accumulate with frozen grid

      njetzero=0
      ncutzero=0
      ntotzero=0
      ntotshot=0
      if (.not.readin .or. .not.dryrun) then
* Initialize efficiency variables      
        call mcfm_vegas(0,itmx1,ncall1,dryrun,integ,integ_err)
        itmxplots=itmx1
      endif
* This is the mcfm_vegas(accum) call
* This takes place only if dryrun is false
* The Vegas parameters are those read from options.DAT for
* the results stage (itmx2,ncall2) and binning takes place (.true.)
* wtmax may have been set during the dry run, so re-set here :
      wtmax = 0d0
      if ((.not.dryrun) .or.
     .    (readin .and. .not.newreadimpl)) then
* Initialize efficiency variables      
        njetzero=0
        ncutzero=0
        ntotzero=0
        ntotshot=0
        call mcfm_vegas(1,itmx2,ncall2,.true.,integ,integ_err)
        itmxplots=itmx2
      endif

      if (newreadimpl) then
        if(.not. readin) then
          open(unit=50, file="CSmax.bin", form="unformatted",
     &      status="replace")
          write(50) wtmax, integ, integ_err
          close(unit=50)
        else
          open(unit=50, file="CSmax.bin", form="unformatted",
     &      status="old")
          read(50) wtmax, integ, integ_err
          close(unit=50)
        endif
      endif

* So far we have not used VEGAS to generate any events.
* Make sure future calls to "getevent" are aware of this :
      numstored = 0

c--- nevtrequested is the number of unweighted events to produce
c--- (so this stage is skipped if nevtrequested <= 0)
      if (nevtrequested .gt. 0) then
        if (part .ne. 'lord') then
          write(6,*) 'LHE events not available beyond LO'
          stop
        endif
        if (doipsgen) then
          write(6,*) 'LHE events not available for this process'
          stop
        endif
        evtgen=.true.
        xs_store_val=integ
        xs_err_store_val=integ_err
        write(6,*) 'Generate events :'
        call mcfm_getunweighted
        call lhefwritefooter(84)
c        do i=1,500
c          call mcfm_getevent(p,wt,pflav,pbarflav)
c          call fill_stdhep(p,0,0,wt)
c         call write_stdhep(6)
c        enddo
      endif

* final processing and print-out
      r=integ
      er=integ_err
      call mcfm_exit(itmxplots,integ,integ_err)
!      stop
      end
       
