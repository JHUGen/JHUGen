      subroutine mcfm_exit(itmx,xinteg,xinteg_err)
      implicit none
      include 'types.f'
************************************************************************
*                                                                      *
*  This routine should perform the final processing and print-outs     *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'mpif.h'
      include 'efficiency.f'
      include 'kprocess.f'
      include 'PDFerrors.f'
      include 'kpart.f'
      include 'outputflags.f'
      include 'bypart.f'
      include 'histo.f'
      include 'taucut.f'
      include 'mpicommon.f'
      integer ierr
      integer(kind=8) collect,sihis(nplot,maxnbin),siuscore(nplot),
     & sioscore(nplot),sient(nplot)
      double precision slord_bypart(-1:1,-1:1),shist(nplot,maxnbin),
     & shavg(nplot),shint(nplot),shsig(nplot)

c--- APPLgrid - to use grids
c      include 'constants.f'
c      include 'ptilde.f'
c      include 'APPLinclude.f'
c--- APPLgrid -end

      logical:: bin
      integer:: j,k,itmx,iu
      real(dp):: xinteg,xinteg_err,minPDFxsec,maxPDFxsec
      real(dp):: PDFarray(0:1000),PDFcentral,PDFerror,PDFperror,
     & PDFnerror
      real(dp):: lordnorm,rescale
      real(dp):: ggpart,gqpart,qgpart,qqpart,qqbpart,
     & gqbpart,qbgpart,qbqbpart,qbqpart
      common/finalpart/ggpart,gqpart,qgpart,qqpart,qqbpart,
     & gqbpart,qbgpart,qbqbpart,qbqpart
      common/bin/bin
      real(dp):: PDFMCav, PDFMCer, sum1,sum2
      character*15 part,kpartstring

      part=kpartstring(kpart)

      if (rank.eq.0) then
c--- Print-out the value of the integral and its error
      write(6,*)
      if (xinteg < 1d7) then
        write(6,53)'Value of final ',trim(part),' integral is',
     &   xinteg,' +/-',xinteg_err, ' fb'
      else
        write(6,53)'Value of final ',trim(part),' integral is',
     &   xinteg/1d6,' +/-',xinteg_err/1d6, ' nb'
        write(6,*) '(WARNING: result in nanobarns)'
      endif

c--- for gg->H+X processes, also write out the cross section
c---  normalized by sigma(gg->H, finite mt)/sigma(gg->H, mt-> infinity)
      if ( (kcase == kggfus0) .or. (kcase == kggfus1)
     & .or.(kcase == kggfus2) .or. (kcase == kggfus3)
     & .or.(kcase == kHWWjet) .or. (kcase == kHWW2jt)
     & .or.(kcase == kHWW3jt) .or. (kcase == kHWW_4l)
     & .or.(kcase == kHWW2lq) .or. (kcase == kHWWdkW)
     & .or.(kcase == kHZZjet) .or. (kcase == kHZZ2jt)
     & .or.(kcase == kHZZ3jt) .or. (kcase == kHZZpjt)
     & .or.(kcase == kHZZ_jj) .or. (kcase == kHZZ_4l)
     & .or.(kcase == kHZZqgI) ) then
        call finitemtcorr(rescale)
        write(6,*)
        write(6,*) 'Cross section normalized by the ratio'
        write(6,*) 'sigma(gg->H, finite mt)/sigma(gg->H, mt-> infinity)'
        write(6,*) '(i.e. exact for gg->H process, but '//
     &               'approx. for gg->H+n jets, n=1,2,3)'
        write(6,*)
        write(6,53)' Rescaled ',trim(part),' integral is',
     &   xinteg*rescale,' +/-',xinteg_err*rescale, ' fb'
        write(6,'(a25,f9.5,a2)') '   (Rescaling factor is ',rescale,')'
      endif

   53 format(a15,a9,a12,G14.6,a4,G13.5,a3)
      endif
c--- Print-out a summary of the effects of jets and cuts

      call mpi_reduce(ntotshot,collect,1,mpi_integer8,
     .     mpi_sum,0,mpi_comm_world,ierr)
      ntotshot=collect
      call mpi_reduce(ntotzero,collect,1,mpi_integer8,
     .     mpi_sum,0,mpi_comm_world,ierr)
      ntotzero=collect
      call mpi_reduce(njetzero,collect,1,mpi_integer8,
     .     mpi_sum,0,mpi_comm_world,ierr)
      njetzero=collect
      call mpi_reduce(ncutzero,collect,1,mpi_integer8,
     .     mpi_sum,0,mpi_comm_world,ierr)
      ncutzero=collect

      if (rank == 0) write(6,*)
      if ((rank == 0) .and. (usescet .eqv. .false.)) then
      if (ntotshot > 0) then
        write(6,*) 'Total number of shots       : ',ntotshot
        write(6,*) 'Total no. failing cuts      : ',ntotzero
        write(6,*) 'Number failing jet cuts     : ',njetzero
        write(6,*) 'Number failing process cuts : ',ncutzero
        write(6,*)
        call flush(6)

c--- Calculate the actual number of shots that were passed
c--- through the jet and cut routines
        ntotshot=ntotshot-(ntotzero-njetzero-ncutzero)
        write(6,54) 'Jet efficiency : ',
     &    100._dp-100._dp*real(njetzero,dp)/real(ntotshot,dp)
        write(6,54) 'Cut efficiency : ',
     &    100._dp-100._dp*real(ncutzero,dp)/real((ntotshot-njetzero),dp)
        write(6,54) 'Total efficiency : ',
     &    100._dp-100._dp*real((njetzero+ncutzero),dp)/real(ntotshot,dp)
        write(6,*)
        endif
      endif

      call mpi_reduce(lord_bypart,slord_bypart,9,mpi_double_precision,
     .     mpi_sum,0,mpi_comm_world,ierr)
      if ((rank.eq.0) .and. (usescet .eqv. .false.))then
      lord_bypart(:,:)=slord_bypart(:,:)
      lordnorm=0._dp
      do j=-1,1
      do k=-1,1
        lordnorm=lordnorm+lord_bypart(j,k)
      enddo
      enddo
      ggpart=lord_bypart( 0, 0)/lordnorm
      gqpart=lord_bypart( 0,+1)/lordnorm
      gqbpart=lord_bypart( 0,-1)/lordnorm
      qgpart=lord_bypart(+1, 0)/lordnorm
      qbgpart=lord_bypart(-1, 0)/lordnorm
      qqpart=lord_bypart(+1,+1)/lordnorm
      qbqbpart=lord_bypart(-1,-1)/lordnorm
      qqbpart=lord_bypart(+1,-1)/lordnorm
      qbqpart=lord_bypart(-1,+1)/lordnorm
      write(6,*) 'Contribution from parton sub-processes:'
      write(6,*) '---------------------------------------'
      write(6,55) '   GG    ',ggpart*xinteg,ggpart*100._dp
      write(6,55) '   GQ    ',gqpart*xinteg,gqpart*100._dp
      write(6,55) '   GQB   ',gqbpart*xinteg,gqbpart*100._dp
      write(6,55) '   QG    ',qgpart*xinteg,qgpart*100._dp
      write(6,55) '   QBG   ',qbgpart*xinteg,qbgpart*100._dp
      write(6,55) '   QQ    ',qqpart*xinteg,qqpart*100._dp
      write(6,55) '   QBQB  ',qbqbpart*xinteg,qbqbpart*100._dp
      write(6,55) '   QQB   ',qqbpart*xinteg,qqbpart*100._dp
      write(6,55) '   QBQ   ',qbqpart*xinteg,qbqpart*100._dp
      write(6,*) '---------------------------------------'
      call flush(6)

   54 format(a20,f6.2,'%')
   55 format(4x,a9,' |',G18.5,f8.2,'%')
      endif
c--- If we've calculated PDF errors, present results using
c--- new implementation of PDF uncertainty (9/2013)

      if (PDFerrors) then
      call mpi_reduce(PDFxsec,PDFarray,1001,
     .     mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
      if (rank.eq.0) then
!        PDFarray(:)=PDFxsec(:)
        call computepdfuncertainty(PDFarray,PDFcentral,
     &   PDFperror,PDFnerror,PDFerror)
c------ loop over output units
        open(unit=91,status='unknown',file='pdfuncertainty.res')
        do iu=6,91,85
        write(iu,*)
        write(iu,58) '********* PDF uncertainty analysis *********'
        write(iu,58) '*                                          *'
        write(iu,57) 'Central value',PDFcentral
        write(iu,58) '*                                          *'
        write(iu,58) '*        Absolute PDF uncertainties        *'
        write(iu,57) '   Symmetric +/-',PDFerror
        write(iu,57) '   +ve direction',PDFperror
        write(iu,57) '   -ve direction',PDFnerror
        write(iu,58) '*                                          *'
        write(iu,58) '*        Relative PDF uncertainties        *'
        write(iu,60) '   Symmetric +/-',PDFerror/PDFcentral*100._dp
        write(iu,60) '   +ve direction',PDFperror/PDFcentral*100._dp
        write(iu,60) '   -ve direction',PDFnerror/PDFcentral*100._dp
        write(iu,58) '*                                          *'
        write(iu,58) '********************************************'
        enddo
        close(91)
      endif

   56 format('* PDF error set ',i3,' -->',f15.3,' fb  *')
   57 format('*   ',a16,f14.3,' fb      *')
   58 format(a44)
   59 format('*   ',a16,f14.3,'         *')
   60 format('*   ',a16,f14.2,' %       *')
      endif

      if (bin) then
c--- Finalize the histograms, if we're not filling ntuples instead
      call mpi_reduce(hist,shist,nplot*maxnbin,mpi_double_precision,
     .     mpi_sum,0,mpi_comm_world,ierr)
         hist(:,:)=shist(:,:)
      call mpi_reduce(havg,shavg,nplot,mpi_double_precision,
     .     mpi_sum,0,mpi_comm_world,ierr)
         havg(:)=shavg(:)
      call mpi_reduce(hint,shint,nplot,mpi_double_precision,
     .     mpi_sum,0,mpi_comm_world,ierr)
         hint(:)=shint(:)
      call mpi_reduce(hsig,shsig,nplot,mpi_double_precision,
     .     mpi_sum,0,mpi_comm_world,ierr)
         hsig(:)=shsig(:)
      call mpi_reduce(ihis,sihis,nplot*maxnbin,mpi_integer8,
     .     mpi_sum,0,mpi_comm_world,ierr)
         ihis(:,:)=sihis(:,:)
      call mpi_reduce(iuscore,siuscore,nplot,mpi_integer8,
     .     mpi_sum,0,mpi_comm_world,ierr)
         iuscore(:)=siuscore(:)
      call mpi_reduce(ioscore,sioscore,nplot,mpi_integer8,
     .     mpi_sum,0,mpi_comm_world,ierr)
         ioscore(:)=sioscore(:)
      call mpi_reduce(ient,sient,nplot,mpi_integer8,
     .     mpi_sum,0,mpi_comm_world,ierr)
         ient(:)=sient(:)
         if (creatent .eqv. .false.) then
            if (dswhisto .eqv. .false.) then
c---  Traditional MCFM histograms
               if (rank.eq.0) call histofin(xinteg,xinteg_err,0,itmx)
            else
c---  DSW histograms - store the information
               call dswhbook(200,'Sigma   ',one,zip,ten)
               call dswhfill(200,0.5_dp,xinteg)
               call dswhfill(200,1.5_dp,xinteg_err)
c---  DSW histograms - output and close file
               call NTfinalize
            endif
         endif
      endif
c--- APPLgrid - creating grid
c      if (creategrid)then
c       call write_grid(xinteg)
c      endif
c--- APPLgrid - end

c--- Write out references
      if (rank == 0) then
        call writereference()
      endif

      return

      end
