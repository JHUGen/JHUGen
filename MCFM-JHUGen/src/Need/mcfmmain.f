      subroutine mcfmmain(inputfile,workdir,r,er)
      implicit none
      include 'types.f'
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
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'efficiency.f'
      include 'maxwt.f'
      include 'eventbuffer.f'
      include 'gridinfo.f'
      include 'kpart.f'
      include 'xs_store_info.f'
      include 'vegas_common.f'
      include 'ipsgen.f'
      include 'iterat.f'
      include 'mpicommon.f'
      real(dp):: integ,integ_err,r,er
      logical:: dryrun
      integer:: itmxplots
      character*72 inputfile,workdir
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

      if ((dryrun .eqv. .false.) .or. 
     &    ((dryrun) .and. (readin .eqv. .false.))) then
* Initialize efficiency variables      
        njetzero=0
        ncutzero=0
        ntotzero=0
        ntotshot=0
        call mcfm_vegas(0,itmx1,ncall1,dryrun,integ,integ_err)
        itmxplots=itmx1
      endif
* This is the mcfm_vegas(accum) call
* This takes place only if dryrun is false
* The Vegas parameters are those read from options.DAT for
* the results stage (itmx2,ncall2) and binning takes place (.true.)
* wtmax may have been set during the dry run, so re-set here :
      wtmax = 0._dp
      if ((dryrun .eqv. .false.) .or. 
     &    ((dryrun) .and. (readin .eqv. .true.))) then
* Initialize efficiency variables      
        njetzero=0
        ncutzero=0
        ntotzero=0
        ntotshot=0
        call mcfm_vegas(1,itmx2,ncall2,.true.,integ,integ_err)
        itmxplots=itmx2
      endif
      
* So far we have not used VEGAS to generate any events.
* Make sure future calls to "getevent" are aware of this :
      numstored = 0

c--- nevtrequested is the number of unweighted events to produce
c--- (so this stage is skipped if nevtrequested <= 0)
      if (nevtrequested > 0) then
         if (rank.eq.0) then
            if (kpart.ne.klord) then
               write(6,*) 'LHE events not available beyond LO'
               stop
            endif
            if (doipsgen) then
               write(6,*) 'LHE events not available for this process'
               stop
            endif
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
       
