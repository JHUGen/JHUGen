      subroutine mcfm_exit(itmx,xinteg,xinteg_err)
************************************************************************
*                                                                      *
*  This routine should perform the final processing and print-outs     *
*                                                                      *
************************************************************************
      implicit none
      include 'efficiency.f'
      include 'process.f'
      include 'PDFerrors.f'
      include 'part.f'
      include 'outputflags.f'
      include 'bypart.f'

c--- APPLgrid - to use grids
c      include 'constants.f'
c      include 'ptilde.f'
c      include 'APPLinclude.f'
c--- APPLgrid -end
      
      logical bin
      integer j,k,itmx,iu
      double precision xinteg,xinteg_err,minPDFxsec,maxPDFxsec
      double precision PDFarray(0:1000),PDFcentral,PDFerror,PDFperror,
     & PDFnerror
      double precision lordnorm,rescale
      double precision ggpart,gqpart,qgpart,qqpart,qqbpart,
     & gqbpart,qbgpart,qbqbpart,qbqpart
      common/finalpart/ggpart,gqpart,qgpart,qqpart,qqbpart,
     & gqbpart,qbgpart,qbqbpart,qbqpart
      common/bin/bin
      double precision PDFMCav, PDFMCer, sum1,sum2
      
c--- Print-out the value of the integral and its error
      write(6,*) 
      if (xinteg .lt. 5d7) then 
        write(6,53)'Value of final ',part,' integral is',
     .   xinteg,' +/-',xinteg_err, ' fb'
      else 
        write(6,53)'Value of final ',part,' integral is',
     .   xinteg/1d6,' +/-',xinteg_err/1d6, ' nb'
        write(6,*) '(WARNING: result in nanobarns)'
      endif

c--- for gg->H+X processes, also write out the cross section
c---  normalized by sigma(gg->H, finite mt)/sigma(gg->H, mt-> infinity)
      if (((case(1:5) .eq. 'ggfus') .or. (case(1:3) .eq. 'HWW')
     & .or.(case(1:3) .eq. 'HZZ')) .and. (case .ne. 'HWWint')
     &  .and. (case .ne. 'HWW_tb') .and. (case .ne. 'HWWH+i')
     &  .and. (case .ne. 'HZZint') .and. (case .ne. 'HZZ_tb') 
     &  .and. (case .ne. 'HZZH+i') ) then
        call finitemtcorr(rescale)
        write(6,*)
      write(6,*) 'Cross section normalized by the ratio'
      write(6,*) 'sigma(gg->H, finite mt)/sigma(gg->H, mt-> infinity)'
      write(6,*) '(i.e. exact for gg->H process, but '//
     .               'approx. for gg->H+n jets, n=1,2,3)'
        write(6,*)
        write(6,53)' Rescaled ',part,' integral is',
     .   xinteg*rescale,' +/-',xinteg_err*rescale, ' fb'   
        write(6,'(a25,f7.3,a2)') '   (Rescaling factor is ',rescale,')'  
      endif
     
   53 format(a15,a4,a12,G14.6,a4,G14.6,a3)

c--- Print-out a summary of the effects of jets and cuts
      write(6,*) 
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
     .  100d0-100d0*dfloat(njetzero)/dfloat(ntotshot)
      write(6,54) 'Cut efficiency : ',
     .  100d0-100d0*dfloat(ncutzero)/dfloat((ntotshot-njetzero))
      write(6,54) 'Total efficiency : ',
     .  100d0-100d0*dfloat((njetzero+ncutzero))/dfloat(ntotshot)
      write(6,*) 
      
      lordnorm=0d0
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
      write(6,55) '   GG    ',ggpart*xinteg,ggpart*100d0
      write(6,55) '   GQ    ',gqpart*xinteg,gqpart*100d0
      write(6,55) '   GQB   ',gqbpart*xinteg,gqbpart*100d0
      write(6,55) '   QG    ',qgpart*xinteg,qgpart*100d0
      write(6,55) '   QBG   ',qbgpart*xinteg,qbgpart*100d0
      write(6,55) '   QQ    ',qqpart*xinteg,qqpart*100d0
      write(6,55) '   QBQB  ',qbqbpart*xinteg,qbqbpart*100d0
      write(6,55) '   QQB   ',qqbpart*xinteg,qqbpart*100d0
      write(6,55) '   QBQ   ',qbqpart*xinteg,qbqpart*100d0
      write(6,*) '---------------------------------------'
      call flush(6)

   54 format(a20,f6.2,'%')
   55 format(4x,a9,' |',G18.5,f8.2,'%')

c--- If we've calculated PDF errors, present results using   
c--- new implementation of PDF uncertainty (9/2013)
      if (PDFerrors) then
        PDFarray(:)=PDFxsec(:)
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
        write(iu,60) '   Symmetric +/-',PDFerror/PDFcentral*100d0
        write(iu,60) '   +ve direction',PDFperror/PDFcentral*100d0
        write(iu,60) '   -ve direction',PDFnerror/PDFcentral*100d0
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
 
      if (bin) then
c--- Finalize the histograms, if we're not filling ntuples instead
         if (creatent .eqv. .false.) then
            if (dswhisto .eqv. .false.) then
c---  Traditional MCFM histograms
               call histofin(xinteg,xinteg_err,0,itmx)
            else
c---  DSW histograms - store the information
               call dswhbook(200,'Sigma   ',1.0d0,0.0d0,10.0d0)
               call dswhfill(200,0.5d0,xinteg)
               call dswhfill(200,1.5d0,xinteg_err)
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

      return
      
      end
