      subroutine mcfm_vegas(myinit,myitmx,myncall,mybin,xinteg,xerr)
************************************************************************
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
*    This version is modified to work with photon processes.           *
*    At NLO one can have:                                              *
*                                                                      *
*    virt -> Born + virtual corrections + integrated QCD counterterms  *
*    real -> Real + subtraction terms                                  *  
*                [+ QED photon-quark subtractions if frag=.true.]      *
*    frag -> Fragmentation contributions, incl. integrated QED dipoles *
*                                                                      *
*    tota -> virt + real [frag=.false.]                                *
*    tota -> virt + real + frag [frag=.true.]                          *
*                                                                      *
************************************************************************

      implicit none
      include 'gridinfo.f'
      include 'realwt.f'
      include 'scale.f'
      include 'facscale.f'
      include 'vegas_common.f'
      include 'PDFerrors.f'
      include 'frag.f'
      include 'reset.f'
      include 'masses.f'
      include 'process.f'
      include 'ipsgen.f'
      include 'part.f'
      integer myitmx,myncall,myinit,i,j,k,nproc,mynproc,myncall_save
      logical mybin,bin
      double precision sig,sd,chi,sigr,sdr,sigdk,sddk,chidk,
     & sigfrag,sdfrag,chifrag,sigWdk,sdWdk,chiWdk,
     & xreal,xreal2,xinteg,xerr,adjust,myscale,myfacscale,
     & mymb,sumsig,sumsigr,sumsigf,sumsd,sumsdr,sumsdf,
     & sigips(4),sdips(4),xcallwt
      character*4 mypart
      character*3 getstr,psgen
      common/nproc/nproc
      common/mypart/mypart
      common/bin/bin
      common/xreal/xreal,xreal2
      double precision lowint,virtint,realint,fragint
      double precision region(2*mxdim),lord_bypart(-1:1,-1:1)
      logical first,myreadin
      common/bypart/lord_bypart
      external lowint,virtint,realint,fragint
      data first/.true./
      save first,sigips,sdips
           
c--- Initialize all integration results to zero, so that the
c--- total of virt and real may be combined at the end for 'tota'
      sig=0d0
      sigr=0d0
      sigdk=0d0
      sigWdk=0d0
      sigfrag=0d0
      sd=0d0
      sdr=0d0
      sddk=0d0
      sdWdk=0d0
      sdfrag=0d0
      xreal=0d0
      xreal2=0d0
      
c--- integer controlling stages of Z+gamma+jet and Z+gamma+gamma processes
      ipsgen=1
      
      do j=-1,1
      do k=-1,1
        lord_bypart(j,k)=0d0
      enddo
      enddo
      if (PDFerrors) then
        do i=0,maxPDFsets
          PDFxsec(i)=0d0
        enddo
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
      
c--- Store value of part in mypart, which will be retained;
c--- also store value of scale in myscale, which will be retained;
c--- part and scale can be changed to make sure that the tota option works.
      mypart=part
      myscale=scale
      myfacscale=facscale
      mynproc=nproc
      mymb=mb
      myncall_save=myncall

c--- special process: 5FNS + 4FNS: initialization
c---   for process 421 (426) we will actually run
c---   process 411 (416) [5FNS] and then add on the result of process 401 (406) [4FNS]
      if ((mynproc .eq. 421) .or. (mynproc .eq. 426)) then
        nproc=mynproc-10
      sumsig=0d0
      sumsigr=0d0
      sumsd=0d0
      sumsdr=0d0
      endif

c--- special handling of multiple PS generation for
c--- Z+gamma+jet and Z+gamma+gamma processes: we will loop
c--- over ipsgen=1,..,n where n=2 for Z+gamma+jet and n=4 for  Z+gamma+gamma
      if ( (mynproc .eq. 301) .or. (mynproc .eq. 302)
     & .or.(mynproc .eq. 303) .or. (mynproc .eq. 304)
     & .or.(mynproc .eq. 370) .or. (mynproc .eq. 371)
     & .or.(doipsgen)) then
      sumsig=0d0
      sumsigr=0d0
      sumsigf=0d0
      sumsd=0d0
      sumsdr=0d0
      sumsdf=0d0
      endif

   77 continue     
      if ((mynproc .eq. 421) .or. (mynproc .eq. 426)) then
        mb=mymb
      call chooser
        reset=.true.
        scalereset=.true.
c--- special write-out/read-in for 5FNS + 4FNS process
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
        if (nproc .ge. 411) then
            outgridfile='dvegas_5FNS_'//part//'.grid'
        else
            outgridfile='dvegas_4FNS_'//part//'.grid'
        endif        
        else
          readin=.true.
          writeout=.false.
        if (nproc .ge. 411) then
            ingridfile='dvegas_5FNS_'//part//'.grid'
        else
            ingridfile='dvegas_4FNS_'//part//'.grid'
        endif        
        endif
      endif    
      
c--- special handling of multiple PS generation for
c--- Z+gamma+jet and Z+gamma+gamma processes 
      if ( (mynproc .eq. 301) .or. (mynproc .eq. 302)
     & .or.(mynproc .eq. 303) .or. (mynproc .eq. 304)
     & .or.(mynproc .eq. 370) .or. (mynproc .eq. 371)
     & .or.(doipsgen)) then
        reset=.true.
        scalereset=.true.
      psgen=getstr(ipsgen)
      write(6,*) '********* Phase space region ',ipsgen,' *********'
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
          outgridfile='dvegas_'//part//'_PS'//psgen(1:1)//'.grid'
        else
          readin=.true.
          writeout=.false.
          ingridfile='dvegas_'//part//'_PS'//psgen(1:1)//'.grid'
        endif
      endif
      
c--- Basic lowest-order integration
      if (part .eq. 'lord') then
       call boundregion(ndim,region)
       call vegasnr(region,ndim,lowint,myinit,myncall,myitmx,
     .               nprn,sig,sd,chi)
      endif


c---- REMOVE THIS PIECE EVENTUALLY
      if (part .eq. 'frit') then
         ndim=ndim+1
         call boundregion(ndim,region) 
         call vegasnr(region,ndim,lowint,myinit,myncall,myitmx,
     &                 nprn,sig,sd,chi)
         ndim=ndim-1
      endif
c---- REMOVE THIS PIECE EVENTUALLY

c--- If we're doing the tota integration, then set up the grid info
      if ((mypart .eq. 'tota') .or. (mypart .eq. 'todk')) then        
        if (first .and. (myinit .eq. 1)) then
c-- special input name for virtual grid
            ingridfile='dvegas_virt_'//ingridfile
            myreadin=readin
        else
          if (first .eqv. .true.) then
            readin=.false.
            writeout=.true.
            outgridfile='dvegas_virt.grid'          
            if ( (mynproc .eq. 301) .or. (mynproc .eq. 302)
     &       .or.(mynproc .eq. 303) .or. (mynproc .eq. 304)
     &       .or.(mynproc .eq. 370) .or. (mynproc .eq. 371)
     &       .or.(doipsgen)) then
              outgridfile='dvegas_virt_PS'//psgen(1:1)//'.grid'
          endif
          else
            readin=.true.
            writeout=.false.
            ingridfile='dvegas_virt.grid'
            if ( (mynproc .eq. 301) .or. (mynproc .eq. 302)
     &       .or.(mynproc .eq. 303) .or. (mynproc .eq. 304)
     &       .or.(mynproc .eq. 370) .or. (mynproc .eq. 371)
     &       .or.(doipsgen)) then
              ingridfile='dvegas_virt_PS'//psgen(1:1)//'.grid'
          endif
          endif
        endif
      endif        

c--- Virtual integration should have one extra dimension
c--- (added and then taken away)
      if (  (mypart .eq. 'virt') .or. (mypart .eq. 'tota')
     . .or. (mypart .eq. 'todk') )  then
        part='virt'        
        reset=.true.
        scalereset=.true.
        ndim=ndim+1
        call boundregion(ndim,region)
        call vegasnr(region,ndim,virtint,myinit,myncall,myitmx,
     .              nprn,sig,sd,chi)
        ndim=ndim-1
      endif
            
c--- If we're doing the tota integration, then set up the grid info
      if ((mypart .eq. 'tota') .or. (mypart .eq. 'todk')) then
        if (first .and. (myinit .eq. 1)) then
c-- special input name for real grid
          ingridfile(8:11)='real'
          readin=myreadin
        else
          if (first .eqv. .true.) then
            readin=.false.
            writeout=.true.
            outgridfile='dvegas_real.grid'          
            if ( (mynproc .eq. 301) .or. (mynproc .eq. 302)
     &       .or.(mynproc .eq. 303) .or. (mynproc .eq. 304)
     &       .or.(mynproc .eq. 370) .or. (mynproc .eq. 371)
     &       .or.(doipsgen)) then
              outgridfile='dvegas_real_PS'//psgen(1:1)//'.grid'
          endif
          else
            readin=.true.
            writeout=.false.
            ingridfile='dvegas_real.grid'
            if ( (mynproc .eq. 301) .or. (mynproc .eq. 302)
     &       .or.(mynproc .eq. 303) .or. (mynproc .eq. 304)
     &       .or.(mynproc .eq. 370) .or. (mynproc .eq. 371)
     &       .or.(doipsgen)) then
              ingridfile='dvegas_real_PS'//psgen(1:1)//'.grid'
          endif
          endif
        endif        
      endif 

c--- Real integration should have three extra dimensions
c--- 'realwt' is a special option that in general should be false
c--- ('realwt' true samples the integral according to the
c---   unsubtracted real emission weight)
      if (mypart .eq. 'real') then
        part='real'
        scalereset=.true.
        if (realwt) then
          nprn=0
        endif
        xreal=0d0
        xreal2=0d0
        ndim=ndim+3
        call boundregion(ndim,region)
        call vegasnr(region,ndim,realint,myinit,myncall,myitmx,
     .              nprn,sigr,sdr,chi)
        ndim=ndim-3
        write(6,*) 
        ncall=myncall
        if (realwt) then
          sigr=xreal
          sdr=dsqrt(abs((xreal2-xreal**2)/dfloat(ncall)))
          write(6,*) itmx,' iterations of ',ncall,' calls'
          write(6,*) 'Value of subtracted integral',sigr
          write(6,*) 'Error on subtracted integral',sdr
        endif
      endif
      if ((mypart .eq. 'tota') .or. (mypart .eq. 'todk')) then
        scale=myscale
        facscale=myfacscale
        part='real'
        reset=.true.
        if (realwt) then
          nprn=0
        endif
        xreal=0d0
        xreal2=0d0
        adjust=(dfloat(ndim+3))/(dfloat(ndim+1))
        ncall=int(dfloat(myncall)**adjust)/2
c--- cap on number of points, may be ndim-dependent, to ensure PS not too large
c        if (ndim .lt. 20) then        
c        if (ncall .gt. myncall*20) ncall=myncall*20
c      else
        if (ncall .gt. myncall*10) ncall=myncall*10
c      endif
        write(6,*) 'Adjusting number of points for real to',ncall
        ndim=ndim+3
        call boundregion(ndim,region)
        call vegasnr(region,ndim,realint,myinit,ncall,myitmx,
     .              nprn,sigr,sdr,chi)
        ndim=ndim-3
        write(6,*) 
        ncall=myncall

        if (realwt) then
          sigr=xreal
          sdr=dsqrt(abs((xreal2-xreal**2)/dfloat(ncall)))
          write(6,*) itmx,' iterations of ',ncall,' calls'
          write(6,*) 'Value of subtracted integral',sigr
          write(6,*) 'Error on subtracted integral',sdr
        endif
      endif      

c--- If we're doing the todk integration, then set up the grid info
      if (mypart .eq. 'todk') then
        if (first .and. (myinit .eq. 1)) then
c-- special input name for real grid
          ingridfile(8:11)='redk'
          readin=myreadin
        else
          if (first .eqv. .true.) then
            readin=.false.
            writeout=.true.
            outgridfile='dvegas_redk.grid'          
          else
            readin=.true.
            writeout=.false.
            ingridfile='dvegas_redk.grid'
          endif
        endif        
      endif 
      
      if (mypart .eq. 'todk')  then
        scale=myscale
        nproc=nproc+1
        call chooser
        part='real'
        reset=.true.
        if (realwt) then
          nprn=0
        endif
        xreal=0d0
        xreal2=0d0
        adjust=(dfloat(ndim+3))/(dfloat(ndim+1))
        ncall=int(dfloat(myncall)**adjust)/2
        write(6,*) 'Adjusting number of points for real to',ncall
        ndim=ndim+3
        call boundregion(ndim,region)
        call vegasnr(region,ndim,realint,myinit,ncall,myitmx,
     .              nprn,sigdk,sddk,chidk)
        ndim=ndim-3
        write(6,*) 
        ncall=myncall
        nproc=nproc-1
        call chooser

        if (realwt) then
          sigdk=xreal
          sddk=dsqrt(abs((xreal2-xreal**2)/dfloat(ncall)))
          write(6,*) itmx,' iterations of ',ncall,' calls'
          write(6,*) 'Value of subtracted integral',sigdk
          write(6,*) 'Error on subtracted integral',sddk
        endif
      endif      

c--- If we're doing the todk integration for ttbar production with
c--- a hadronic W decay, then set up the grid info
      if ((mypart .eq. 'todk') .and. (case .eq. 'tt_bbh')) then
        if (first .and. (myinit .eq. 1)) then
c-- special input name for real grid
          ingridfile(8:11)='rWdk'
          readin=myreadin
        else
          if (first .eqv. .true.) then
            readin=.false.
            writeout=.true.
            outgridfile='dvegas_rWdk.grid'          
          else
            readin=.true.
            writeout=.false.
            ingridfile='dvegas_rWdk.grid'
          endif
        endif        
      endif 
      
      if ((mypart .eq. 'todk') .and. (case .eq. 'tt_bbh')) then
        scale=myscale
        nproc=nproc+2
        call chooser
        part='real'
        reset=.true.
        if (realwt) then
          nprn=0
        endif
        xreal=0d0
        xreal2=0d0
        adjust=(dfloat(ndim+3))/(dfloat(ndim+1))
        ncall=int(dfloat(myncall)**adjust)/2
        write(6,*) 'Adjusting number of points for real to',ncall
        ndim=ndim+3
        call boundregion(ndim,region)
        call vegasnr(region,ndim,realint,myinit,ncall,myitmx,
     .              nprn,sigWdk,sdWdk,chiWdk)
        ndim=ndim-3
        write(6,*) 
        ncall=myncall
        nproc=nproc-2
        call chooser

        if (realwt) then
          sigWdk=xreal
          sdWdk=dsqrt(abs((xreal2-xreal**2)/dfloat(ncall)))
          write(6,*) itmx,' iterations of ',ncall,' calls'
          write(6,*) 'Value of subtracted integral',sigdk
          write(6,*) 'Error on subtracted integral',sddk
        endif
      endif      

c--- Fragmentation contribution: only if frag is set to .true.
c--- If we're doing the tota integration and , then set up the grid info
      if ((mypart .eq. 'tota') .and. (frag)) then
        if (first .and. (myinit .eq. 1)) then
c-- special input name for frag grid
          ingridfile(8:11)='frag'
          readin=myreadin
        else
          if (first .eqv. .true.) then
            readin=.false.
            writeout=.true.
            outgridfile='dvegas_frag.grid'          
            if ( (mynproc .eq. 301) .or. (mynproc .eq. 302)
     &       .or.(mynproc .eq. 303) .or. (mynproc .eq. 304)
     &       .or.(mynproc .eq. 370) .or. (mynproc .eq. 371)
     &       .or.(doipsgen)) then
              outgridfile='dvegas_frag_PS'//psgen(1:1)//'.grid'
          endif
          else
            readin=.true.
            writeout=.false.
            ingridfile='dvegas_frag.grid'
            if ( (mynproc .eq. 301) .or. (mynproc .eq. 302)
     &       .or.(mynproc .eq. 303) .or. (mynproc .eq. 304)
     &       .or.(mynproc .eq. 370) .or. (mynproc .eq. 371)
     &       .or.(doipsgen)) then
              ingridfile='dvegas_frag_PS'//psgen(1:1)//'.grid'
            endif
          endif
        endif        
      endif 
      
c--- special handling of multiple PS generation for
c--- Z+gamma+jet and Z+gamma+gamma processes:
c---  fragmentation pieces only contribute for
c---   ipsgen=1 (Z+gamma+jet)
c---   ipsgen=1,2 (Z+gamma+gamma) 
c---   ipsgen=1,3 (W+gamma+gamma) 
      if (((mynproc .eq. 301) .and. (ipsgen .gt. 2)) .or. 
     &    ((mynproc .eq. 302) .and. (ipsgen .gt. 1)) .or.
     &    ((mynproc .eq. 370) .and. ((ipsgen .eq. 2).or.(ipsgen.eq.4))) 
     & .or.((mynproc .eq. 371) .and. ((ipsgen .eq. 2).or.(ipsgen.eq.4)))
     &     )then
        sigfrag=0d0
        sdfrag=0d0
        goto 33
      endif

      if ( ((mypart .eq. 'tota') .and. (frag))
     & .or. (mypart .eq. 'frag') ) then
         part='frag'
         scale=myscale
         facscale=myfacscale
         if (mypart .eq. 'frag') then
         scalereset=.true.
         else
         reset=.true.
         endif
         ncall=myncall
c         write(6,*) 'Adjusting number of points for frag to',ncall
         ndim=ndim+1
         call boundregion(ndim,region) 
         fragint_mode=.true.            ! for isolation
         call vegasnr(region,ndim,fragint,myinit,myncall,myitmx,
     &                 nprn,sigfrag,sdfrag,chifrag)
         fragint_mode=.false.            ! for isolation
         ndim=ndim-1
       rescale=.false.       ! turn rescaling off again
      endif

   33 continue 

c--- special process: 5FNS + 4FNS: initialization
c---   for process 421 (426) we will actually run
c---   process 411 (416) [5FNS] and then add on the result of process 401 (406) [4FNS]
      if ((mynproc .eq. 421) .or. (mynproc .eq. 426)) then 
        sumsig=sumsig+sig
        sumsigr=sumsigr+sigr
        sumsd=sumsd+sd**2
        sumsdr=sumsdr+sdr**2
        if (nproc .ge. 411) then
          nproc=nproc-10
        goto 77
      endif
c--- return nproc to the value from the input file
        nproc=mynproc       
      sig=sumsig
      sigr=sumsigr
      sd=dsqrt(sumsd)
      sdr=dsqrt(sumsdr)
      endif
      
c--- special handling of multiple PS generation for
c--- Z+gamma+jet and Z+gamma+gamma processes: we will loop
c--- over ipsgen=1,..,n where n=2 for Z+gamma+jet and n=4 for  Z+gamma+gamma
      if ( (mynproc .eq. 301) .or. (mynproc .eq. 302)
     & .or.(mynproc .eq. 303) .or. (mynproc .eq. 304)
     & .or.(mynproc .eq. 370) .or. (mynproc .eq. 371)
     & .or.(doipsgen)) then
c--- store cross section and standard deviation info for this region on first run
        if (myinit .eq. 0) then
          sigips(ipsgen)=sig+sigr+sigfrag
          sdips(ipsgen)=sqrt(sd**2+sdr**2+sdfrag**2)
        endif
        sumsig=sumsig+sig
        sumsigr=sumsigr+sigr
        sumsigf=sumsigf+sigfrag
        sumsd=sumsd+sd**2
        sumsdr=sumsdr+sdr**2
        sumsdf=sumsdf+sdfrag**2
        ipsgen=ipsgen+1
        if ((((nproc .eq. 302).or.(nproc .eq. 304)).and.(ipsgen .le. 2))
     &  .or.(((nproc .eq. 301).or.(nproc .eq. 303)).and.(ipsgen .le. 4))
     &  .or.(((nproc .eq. 370).or.(nproc .eq. 371)).and.(ipsgen .le. 4))
     &  .or.((doipsgen).and.(ipsgen .le. maxipsgen))
     &     ) then
c--- W+2 photons case
          if ((nproc .eq. 370).or.(nproc .eq. 371)) then
            if ((myinit .eq. 1) .and. (ipsgen .gt. 1)) then
c--- dynamic scaling of points on second run according to uncertainties found in first
              xcallwt=(sdips(ipsgen)/sdips(1))
              if (xcallwt .lt. 0.2d0) xcallwt=0.2d0
              if (xcallwt .gt. 5d0) xcallwt=5d0
              myncall=int(float(myncall_save)*xcallwt)
c              write(6,*) '>>>>>> Dynamic scaling factor = ',
c     &                    xcallwt,' for ipsgen=',ipsgen
            endif
            goto 77 
          endif
c--- change number of PS points in different regions
          if ((nproc.eq.301).or.(nproc.eq.303)
     &    .or.(nproc.eq.370).or.(nproc.eq.371)) then
             if (ipsgen.eq.2) then
                myncall=10*myncall
             elseif (ipsgen.eq.3) then
                myncall=myncall/10
                myncall=5*myncall
             elseif (ipsgen.eq.4) then
                myncall=myncall
             endif
          endif
          myncall=myncall
          goto 77
      endif
      myncall=myncall
      sig=sumsig
      sigr=sumsigr
      sigfrag=sumsigf
      sd=dsqrt(sumsd)
      sdr=dsqrt(sumsdr)
      sdfrag=dsqrt(sumsdf)
      endif
      
c--- calculate integration variables to be returned
      xinteg=sig+sigr+sigdk+sigWdk+sigfrag
      xerr=dsqrt(sd**2+sdr**2+sddk**2+sdWdk**2+sdfrag**2)      
      
c--- return part, scale and myncall to their real values
      part=mypart
      scale=myscale
      myncall=myncall_save
      first=.false.
      
      return
      end
      
      
      subroutine boundregion(idim,region)
c--- Initializes integration region [0,1] for each variable
c--- in the idim-dimensional integration range
      implicit none
      include 'mxdim.f'
      integer i,idim
      double precision region(2*mxdim)
      
      do i=1,idim
      region(i)=0d0
      region(i+idim)=1d0
      enddo
      
      return
      end
      
      
