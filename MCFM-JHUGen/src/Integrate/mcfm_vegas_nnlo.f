      subroutine mcfm_vegas_nnlo(myinit,myitmx,myncall,mybin,xinteg,xerr)
      implicit none
      include 'types.f'
************************************************************************
*                                                                      *
*  Modified version of mcfm_vegas that computes NNLO                   *
*  cross-section using dipole subtraction for NLO contribution         *
*  and SCET for the NNLO coefficient                                   *
*                                                                      *
*  This routine should perform the sweeps of vegasnr                   *
*                                                                      *
*    Input parameters:                                                 *
*       myinit  :  the vegasnr routine entry point                     *
*       myitmx  :  the number of vegasnr sweeps                        *
*      myncall  :  the number of iterations per sweep                  *
*          bin  :  whether or not the results should be histogrammed   *
*                                                                      *
*    Returned variables:                                               *
*       xinteg  :  value of integration                                *
*         xerr  :  integration error                                   *
*                                                                      *
*    coeffonly =  .true.   -->  NNLO coefficient only                  *
*    coeffonly =  .false.  -->  full NNLO result                       *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'gridinfo.f'
      include 'realwt.f'
      include 'scale.f'
      include 'facscale.f'
      include 'vegas_common.f'
      include 'PDFerrors.f'
      include 'frag.f'
      include 'reset.f'
      include 'masses.f'
      include 'kprocess.f'
      include 'ipsgen.f'
      include 'kpart.f'
      include 'taucut.f'
      include 'nproc.f'
      include 'mpicommon.f'
      integer myitmx,myinit,i,j,k,mynproc,nprocabove
      integer(kind=8) myncall,myncall_save
      integer:: ierr
      logical:: mybin,bin,coeffonly_save
      real(dp):: sigr,sdr,sigdk,sddk,chidk,
     & sigfrag,sdfrag,chifrag,sigWdk,sdWdk,chiWdk,
     & xreal,xreal2,xinteg,xerr,adjust,myscale,myfacscale,
     & mymb,sumsig,sumsigr,sumsigf,sumsd,sumsdr,sumsdf,
     & sigips(4),sdips(4),xcallwt,
     & sig(5),sd(5),chi(5)
      integer mykpart
      character*3 getstr,psgen
      common/mykpart/mykpart
      common/bin/bin
      common/xreal/xreal,xreal2
      real(dp):: lowint,virtint,realint,fragint,scetint,qtint
      real(dp):: region(2*mxdim),lord_bypart(-1:1,-1:1)
      logical:: first,myreadin
      common/bypart/lord_bypart
      external lowint,virtint,realint,fragint,scetint,qtint
      data first/.true./
      save first,sigips,sdips
           
c--- Initialize all integration results to zero, so that the
c--- total of all contributions may be combined at the end

c--- NLO: virt
c--- NLO: real
c--- NNLO: below cut
c--- NNLO: above cut virt
c--- NNLO: above cut real

      sig(:)=zip
      sd(:)=zip
      
      lord_bypart(:,:)=0d0

      if (PDFerrors) then
        PDFxsec(:)=zip
      endif

c--- Controls behaviour of gen_njets: need to reset phase-space
c--- boundaries when going from virt to real (using tota)
c--- need to reset scale also, for special scalestart values
      reset=.false.
      scalereset=.false.

c--- Put the vegasnr parameters in the common block
      itmx=myitmx
      ncall=myncall
      bin=mybin
      
c--- Store value of kpart in mykpart, which will be retained;
c--- also store value of scale in myscale, which will be retained;
c--- part and scale can be changed to make sure that the tota option works.
      mykpart=kpart
      myscale=scale
      myfacscale=facscale
      mynproc=nproc
      mymb=mb
      myncall_save=myncall

c--- skip NLO calculation if coeffonly=.true.
      if (coeffonly) goto 77
      
c-------------------------------------------------------------------------------
c------------------------------ NLO calculation --------------------------------
c-------------------------------------------------------------------------------

      usescet=.false.
      abovecut=.false.
      coeffonly=.false.

c--- set up the grid info for virtual
      if (first .and. (myinit == 1)) then
c-- special input name for virtual grid
          ingridfile='dvegas_virt_below_'//ingridfile
          myreadin=readin
      else
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
          outgridfile='dvegas_virt_below.grid'
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            outgridfile='dvegas_virt_below_PS'//psgen(1:1)//'.grid'
          endif
        else
          readin=.true.
          writeout=.false.
          ingridfile='dvegas_virt_below.grid'
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            ingridfile='dvegas_virt_below_PS'//psgen(1:1)//'.grid'
          endif
        endif
      endif

c--- Virtual integration should have one extra dimensions
c--- (added and then taken away)
      kpart=kvirt
      reset=.true.
      scalereset=.true.
      ndim=ndim+1
      call boundregion(ndim,region)
      call vegasnr(region,ndim,virtint,myinit,myncall,myitmx,
     &             nprn,sig(1),sd(1),chi(1))
      ndim=ndim-1
            
c--- set up the grid info for real integration
      if (first .and. (myinit == 1)) then
c-- special input name for real grid
        ingridfile(8:11)='real'
        readin=myreadin
      else
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
          outgridfile='dvegas_real_below.grid'          
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            outgridfile='dvegas_real_below_PS'//psgen(1:1)//'.grid'
        endif
        else
          readin=.true.
          writeout=.false.
          ingridfile='dvegas_real_below.grid'
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            ingridfile='dvegas_real_below_PS'//psgen(1:1)//'.grid'
        endif
        endif
      endif

c--- Real integration should have three extra dimensions
      scale=myscale
      facscale=myfacscale
      kpart=kreal
      reset=.true.
      xreal=0d0
      xreal2=0d0
      adjust=(real(ndim+3,dp))/(real(ndim+1,dp))
      ncall=int(real(myncall,dp)**adjust,kind=8)/2
      if (ncall > myncall*10) ncall=myncall*10
      if (rank == 0) then
        write(6,*) 'Adjusting number of points for real to',ncall
      endif
      ndim=ndim+3
      call boundregion(ndim,region)
      call vegasnr(region,ndim,realint,myinit,ncall,myitmx,
     &            nprn,sig(2),sd(2),chi(2))
      ndim=ndim-3
      if (rank == 0) write(6,*) 
      ncall=myncall

   77 continue
c-------------------------------------------------------------------------------
c-------------------------- SCET NNLO calculation ------------------------------
c-------------------------------------------------------------------------------

      coeffonly_save=coeffonly

      usescet=.true.
      coeffonly=.true.
      
c--- set up SCET variables
      call setupscet(nprocabove)

c--- SCET below-cut contribution
c--- integration should have two extra dimensions (added and then taken away)
      if (first .and. (myinit == 1)) then
c-- special input name for SCET grid
          ingridfile='dvegas_scet_below_'//ingridfile
          myreadin=readin
      else
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
          outgridfile='dvegas_scet_below.grid'
        else
          readin=.true.
          writeout=.false.
          ingridfile='dvegas_scet_below.grid'
        endif
      endif
      kpart=mykpart
      reset=.true.
      scalereset=.true.
      ndim=ndim+2
      abovecut=.false.
      call boundregion(ndim,region)
      call vegasnr(region,ndim,scetint,myinit,myncall,myitmx,
     &             nprn,sig(3),sd(3),chi(3))
! Uncomment these lines (and comment-out the above two) to switch
! from jettiness-slicing to QT-slicing
!      call vegasnr(region,ndim,qtint,myinit,myncall,myitmx,
!     &             nprn,sig(3),sd(3),chi(3))
      ndim=ndim-2
      nproc=nprocabove
c--- scale up number of points more for W/Z/H+jet processes
      if ( (nprocabove == 22) .or. (nprocabove == 27)
     &.or. (nprocabove == 44) .or. (nprocabove == 270)
     &.or. (nprocabove == 271) .or. (nprocabove == 272) ) then
        myncall=myncall*20
      else
        myncall=myncall*10
      endif
      abovecut=.true.
      call chooser
      if (first .eqv. .true.) then
        readin=.false.
        writeout=.true.
        outgridfile='dvegas_scet_above.grid'
      else
        readin=.true.
        writeout=.false.
        ingridfile='dvegas_scet_above.grid'
      endif
            
c--- set up the grid info for virtual
      if (first .and. (myinit == 1)) then
c-- special input name for virtual grid
          ingridfile='dvegas_virt_above_'//ingridfile
          myreadin=readin
      else
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
          outgridfile='dvegas_virt_above.grid'
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            outgridfile='dvegas_virt_above_PS'//psgen(1:1)//'.grid'
          endif
        else
          readin=.true.
          writeout=.false.
          ingridfile='dvegas_virt_above.grid'
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            ingridfile='dvegas_virt_above_PS'//psgen(1:1)//'.grid'
          endif
        endif
      endif

c--- Virtual integration should have one extra dimensions
c--- (added and then taken away)
      kpart=kvirt
      reset=.true.
      scalereset=.true.
      ndim=ndim+1
      call boundregion(ndim,region)
      call vegasnr(region,ndim,virtint,myinit,myncall,myitmx,
     &             nprn,sig(4),sd(4),chi(4))
      ndim=ndim-1
            
c--- set up the grid info for real integration
      if (first .and. (myinit == 1)) then
c-- special input name for real grid
        ingridfile(8:11)='real'
        readin=myreadin
      else
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
          outgridfile='dvegas_real_above.grid'          
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            outgridfile='dvegas_real_above_PS'//psgen(1:1)//'.grid'
        endif
        else
          readin=.true.
          writeout=.false.
          ingridfile='dvegas_real_above.grid'
          if ( (mynproc == 301) .or. (mynproc == 302)
     &     .or.(mynproc == 303) .or. (mynproc == 304)
     &     .or.(mynproc == 370) .or. (mynproc == 371)
     &     .or.(doipsgen)) then
            ingridfile='dvegas_real_above_PS'//psgen(1:1)//'.grid'
        endif
        endif
      endif

c--- Real integration should have three extra dimensions
      scale=myscale
      facscale=myfacscale
      kpart=kreal
      reset=.true.
      xreal=0d0
      xreal2=0d0
      adjust=(real(ndim+3,dp))/(real(ndim+1,dp))
      ncall=int(real(myncall,dp)**adjust,kind=8)/2
      if (ncall > myncall*40) ncall=myncall*40
      if (rank == 0) then
        write(6,*) 'Adjusting number of points for real to',ncall
      endif
      ndim=ndim+3
      call boundregion(ndim,region)
      call vegasnr(region,ndim,realint,myinit,ncall,myitmx,
     &            nprn,sig(5),sd(5),chi(5))
      ndim=ndim-3
      if (rank == 0) write(6,*) 
      ncall=myncall

c--- return nproc to the value from the input file
      nproc=mynproc
      if (first) call chooser

c--- calculate integration variables to be returned
      xinteg=zip
      xerr=zip
      do k=1,5
        xinteg=xinteg+sig(k)
        xerr=xerr+sd(k)**2
      enddo
      xerr=sqrt(xerr)
      
c--- return part, scale, myncall and coeffonly to their real values
      kpart=mykpart
      scale=myscale
      myncall=myncall_save
      coeffonly=coeffonly_save
      first=.false.

      return
      end
