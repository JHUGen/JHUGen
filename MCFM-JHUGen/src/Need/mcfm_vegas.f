      subroutine mcfm_vegas(myinit,myitmx,myncall,mybin,xinteg,xerr)
      implicit none
      include 'types.f'
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
      include 'bypart.f'
      include 'nproc.f'
      integer:: myitmx,myncall,myinit,i,j,k,mynproc,myncall_save
      logical:: mybin,bin
      real(dp):: sig,sd,chi,sigr,sdr,sigdk,sddk,chidk,
     & sigfrag,sdfrag,chifrag,sigWdk,sdWdk,chiWdk,
     & xreal,xreal2,xinteg,xerr,adjust,myscale,myfacscale,
     & mymb,sumsig,sumsigr,sumsigf,sumsd,sumsdr,sumsdf,
     & sigips(4),sdips(4),xcallwt
      integer mykpart
      character*3 getstr,psgen
      common/mykpart/mykpart
      common/bin/bin
      common/xreal/xreal,xreal2
      real(dp):: lowint,virtint,realint,fragint
      real(dp):: region(2*mxdim)
      logical:: first,myreadin
      external lowint,virtint,realint,fragint
      data first/.true./
      save first,sigips,sdips
           
c--- Initialize all integration results to zero, so that the
c--- total of virt and real may be combined at the end for ktota
      sig=0._dp
      sigr=0._dp
      sigdk=0._dp
      sigWdk=0._dp
      sigfrag=0._dp
      sd=0._dp
      sdr=0._dp
      sddk=0._dp
      sdWdk=0._dp
      sdfrag=0._dp
      xreal=0._dp
      xreal2=0._dp
c--- integer:: controlling stages of Z+gamma+jet and Z+gamma+gamma processes
      ipsgen=1
      
      do j=-1,1
      do k=-1,1
        lord_bypart(j,k)=0._dp
      enddo
      enddo
      if (PDFerrors) then
        do i=0,maxPDFsets
          PDFxsec(i)=0._dp
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
      
c--- Store value of part in mykpart, which will be retained;
c--- also store value of scale in myscale, which will be retained;
c--- part and scale can be changed to make sure that the tota option works.
      mykpart=part
      myscale=scale
      myfacscale=facscale
      mynproc=nproc
      mymb=mb
      myncall_save=myncall

c--- special process: 5FNS + 4FNS: initialization
c---   for process 421 (426) we will actually run
c---   process 411 (416) [5FNS] and then add on the result of process 401 (406) [4FNS]
      if ((mynproc == 421) .or. (mynproc == 426)) then
        nproc=mynproc-10
      sumsig=0._dp
      sumsigr=0._dp
      sumsd=0._dp
      sumsdr=0._dp
      endif

c--- special handling of multiple PS generation for
c--- Z+gamma+jet and Z+gamma+gamma processes: we will loop
c--- over ipsgen=1,..,n where n=2 for Z+gamma+jet and n=4 for  Z+gamma+gamma
      if ( (mynproc == 301) .or. (mynproc == 302)
     & .or.(mynproc == 303) .or. (mynproc == 304)
     & .or.(mynproc == 370) .or. (mynproc == 371)
     & .or.(doipsgen)) then
      sumsig=0._dp
      sumsigr=0._dp
      sumsigf=0._dp
      sumsd=0._dp
      sumsdr=0._dp
      sumsdf=0._dp
      endif

   77 continue     
      if ((mynproc == 421) .or. (mynproc == 426)) then
        mb=mymb
      call chooser
        reset=.true.
        scalereset=.true.
c--- special write-out/rea.e-_dpin for 5FNS + 4FNS process
        if (first .eqv. .true.) then
          readin=.false.
          writeout=.true.
        if (nproc >= 411) then
            outgridfile='dvegas_5FNS_'//part//'.grid'
        else
            outgridfile='dvegas_4FNS_'//part//'.grid'
        endif        
        else
          readin=.true.
          writeout=.false.
        if (nproc >= 411) then
            ingridfile='dvegas_5FNS_'//part//'.grid'
        else
            ingridfile='dvegas_4FNS_'//part//'.grid'
        endif        
        endif
      endif    
      
c--- special handling of multiple PS generation for
c--- Z+gamma+jet and Z+gamma+gamma processes 
      if ( (mynproc == 301) .or. (mynproc == 302)
     & .or.(mynproc == 303) .or. (mynproc == 304)
     & .or.(mynproc == 370) .or. (mynproc == 371)
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
      if (kpart==klord) then
!       call boundregion(ndim,region)
!       call vegasnr(region,ndim,lowint,myinit,myncall,myitmx,
!     &               nprn,sig,sd,chi)
       call vegas_Pomp(lowint,sig,0._dp,0._dp,ndim,myncall,myitmx,myinit)
      endif


c--- If we're doing the tota integration, then set up the grid info
      if ((mykpart==ktota) .or. (mykpart==ktodk)) then
        if (first .and. (myinit == 1)) then
c-- special input name for virtual grid
            ingridfile='dvegas_virt_'//ingridfile
            myreadin=readin
        else
          if (first .eqv. .true.) then
            readin=.false.
            writeout=.true.
            outgridfile='dvegas_virt.grid'          
            if ( (mynproc == 301) .or. (mynproc == 302)
     &       .or.(mynproc == 303) .or. (mynproc == 304)
     &       .or.(mynproc == 370) .or. (mynproc == 371)
     &       .or.(doipsgen)) then
              outgridfile='dvegas_virt_PS'//psgen(1:1)//'.grid'
          endif
          else
            readin=.true.
            writeout=.false.
            ingridfile='dvegas_virt.grid'
            if ( (mynproc == 301) .or. (mynproc == 302)
     &       .or.(mynproc == 303) .or. (mynproc == 304)
     &       .or.(mynproc == 370) .or. (mynproc == 371)
     &       .or.(doipsgen)) then
              ingridfile='dvegas_virt_PS'//psgen(1:1)//'.grid'
          endif
          endif
        endif
      endif        

c--- Virtual integration should have one extra dimension
c--- (added and then taken away)
      if (  (mykpart==kvirt) .or. (mykpart==ktota)
     & .or. (mykpart==ktodk) )  then
        part=kvirt        
        reset=.true.
        scalereset=.true.
        ndim=ndim+1
!        call boundregion(ndim,region)
!        call vegasnr(region,ndim,virtint,myinit,myncall,myitmx,
!     &              nprn,sig,sd,chi)
       call vegas_Pomp(virtint,sig,0._dp,0._dp,ndim,myncall,myitmx,myinit)
        ndim=ndim-1
      endif
            
c--- If we're doing the tota integration, then set up the grid info
      if ((mykpart==ktota) .or. (mykpart==ktodk)) then
        if (first .and. (myinit == 1)) then
c-- special input name for real grid
          ingridfile(8:11)=kreal
          readin=myreadin
        else
          if (first .eqv. .true.) then
            readin=.false.
            writeout=.true.
            outgridfile='dvegas_real.grid'          
            if ( (mynproc == 301) .or. (mynproc == 302)
     &       .or.(mynproc == 303) .or. (mynproc == 304)
     &       .or.(mynproc == 370) .or. (mynproc == 371)
     &       .or.(doipsgen)) then
              outgridfile='dvegas_real_PS'//psgen(1:1)//'.grid'
          endif
          else
            readin=.true.
            writeout=.false.
            ingridfile='dvegas_real.grid'
            if ( (mynproc == 301) .or. (mynproc == 302)
     &       .or.(mynproc == 303) .or. (mynproc == 304)
     &       .or.(mynproc == 370) .or. (mynproc == 371)
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
      if (mykpart==kreal) then
        part=kreal
        scalereset=.true.
        if (realwt) then
          nprn=0
        endif
        xreal=0._dp
        xreal2=0._dp
        ndim=ndim+3
!        call boundregion(ndim,region)
!        call vegasnr(region,ndim,realint,myinit,myncall,myitmx,
!     &              nprn,sigr,sdr,chi)
        call vegas_Pomp(realint,sig,0._dp,0._dp,ndim,myncall,myitmx,myinit)
        ndim=ndim-3
        write(6,*) 
        ncall=myncall
        if (realwt) then
          sigr=xreal
          sdr=sqrt(abs((xreal2-xreal**2)/real(ncall,dp)))
          write(6,*) itmx,' iterations of ',ncall,' calls'
          write(6,*) 'Value of subtracted integral',sigr
          write(6,*) 'Error on subtracted integral',sdr
        endif
      endif
      if ((mykpart==ktota) .or. (mykpart==ktodk)) then
        scale=myscale
        facscale=myfacscale
        part=kreal
        reset=.true.
        if (realwt) then
          nprn=0
        endif
        xreal=0._dp
        xreal2=0._dp
        adjust=(real(ndim+3,dp))/(real(ndim+1,dp))
        ncall=int(real(myncall,dp)**adjust)/2
c--- cap on number of points, may be ndim-dependent, to ensure PS not too large
c        if (ndim < 20) then        
c        if (ncall > myncall*20) ncall=myncall*20
c      else
        if (ncall > myncall*10) ncall=myncall*10
c      endif
        write(6,*) 'Adjusting number of points for real to',ncall
        ndim=ndim+3
!        call boundregion(ndim,region)
!        call vegasnr(region,ndim,realint,myinit,ncall,myitmx,
!     &              nprn,sigr,sdr,chi)
        call vegas_Pomp(realint,sig,0._dp,0._dp,ndim,ncall,myitmx,myinit)
        ndim=ndim-3
        write(6,*) 
        ncall=myncall

        if (realwt) then
          sigr=xreal
          sdr=sqrt(abs((xreal2-xreal**2)/real(ncall,dp)))
          write(6,*) itmx,' iterations of ',ncall,' calls'
          write(6,*) 'Value of subtracted integral',sigr
          write(6,*) 'Error on subtracted integral',sdr
        endif
      endif      

c--- If we're doing the todk integration, then set up the grid info
      if (mykpart==ktodk) then
        if (first .and. (myinit == 1)) then
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
      
      if (mykpart==ktodk)  then
        scale=myscale
        nproc=nproc+1
        call chooser
        part=kreal
        reset=.true.
        if (realwt) then
          nprn=0
        endif
        xreal=0._dp
        xreal2=0._dp
        adjust=(real(ndim+3,dp))/(real(ndim+1,dp))
        ncall=int(real(myncall,dp)**adjust)/2
        write(6,*) 'Adjusting number of points for real to',ncall
        ndim=ndim+3
!        call boundregion(ndim,region)
!        call vegasnr(region,ndim,realint,myinit,ncall,myitmx,
!     &              nprn,sigdk,sddk,chidk)
        call vegas_Pomp(realint,sig,0._dp,0._dp,ndim,ncall,myitmx,myinit)
        ndim=ndim-3
        write(6,*) 
        ncall=myncall
        nproc=nproc-1
        call chooser

        if (realwt) then
          sigdk=xreal
          sddk=sqrt(abs((xreal2-xreal**2)/real(ncall,dp)))
          write(6,*) itmx,' iterations of ',ncall,' calls'
          write(6,*) 'Value of subtracted integral',sigdk
          write(6,*) 'Error on subtracted integral',sddk
        endif
      endif      

c--- If we're doing the todk integration for ttbar production with
c--- a hadronic W decay, then set up the grid info
      if ((mykpart==ktodk) .and. (kcase==ktt_bbh)) then
        if (first .and. (myinit == 1)) then
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
      
      if ((mykpart==ktodk) .and. (kcase==ktt_bbh)) then
        scale=myscale
        nproc=nproc+2
        call chooser
        part=kreal
        reset=.true.
        if (realwt) then
          nprn=0
        endif
        xreal=0._dp
        xreal2=0._dp
        adjust=(real(ndim+3,dp))/(real(ndim+1,dp))
        ncall=int(real(myncall,dp)**adjust)/2
        write(6,*) 'Adjusting number of points for real to',ncall
        ndim=ndim+3
!        call boundregion(ndim,region)
!        call vegasnr(region,ndim,realint,myinit,ncall,myitmx,
!     &              nprn,sigWdk,sdWdk,chiWdk)
      call vegas_Pomp(realint,sig,0._dp,0._dp,ndim,ncall,myitmx,myinit)
        ndim=ndim-3
        write(6,*) 
        ncall=myncall
        nproc=nproc-2
        call chooser

        if (realwt) then
          sigWdk=xreal
          sdWdk=sqrt(abs((xreal2-xreal**2)/real(ncall,dp)))
          write(6,*) itmx,' iterations of ',ncall,' calls'
          write(6,*) 'Value of subtracted integral',sigdk
          write(6,*) 'Error on subtracted integral',sddk
        endif
      endif      

c--- Fragmentation contribution: only if frag is set to .true.
c--- If we're doing the tota integration and , then set up the grid info
      if ((mykpart==ktota) .and. (frag)) then
        if (first .and. (myinit == 1)) then
c-- special input name for frag grid
          ingridfile(8:11)=kfrag
          readin=myreadin
        else
          if (first .eqv. .true.) then
            readin=.false.
            writeout=.true.
            outgridfile='dvegas_frag.grid'          
            if ( (mynproc == 301) .or. (mynproc == 302)
     &       .or.(mynproc == 303) .or. (mynproc == 304)
     &       .or.(mynproc == 370) .or. (mynproc == 371)
     &       .or.(doipsgen)) then
              outgridfile='dvegas_frag_PS'//psgen(1:1)//'.grid'
          endif
          else
            readin=.true.
            writeout=.false.
            ingridfile='dvegas_frag.grid'
            if ( (mynproc == 301) .or. (mynproc == 302)
     &       .or.(mynproc == 303) .or. (mynproc == 304)
     &       .or.(mynproc == 370) .or. (mynproc == 371)
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
      if (((mynproc == 301) .and. (ipsgen > 2)) .or. 
     &    ((mynproc == 302) .and. (ipsgen > 1)) .or.
     &    ((mynproc == 370) .and. ((ipsgen == 2).or.(ipsgen==4))) 
     & .or.((mynproc == 371) .and. ((ipsgen == 2).or.(ipsgen==4)))
     &     )then
        sigfrag=0._dp
        sdfrag=0._dp
        goto 33
      endif

      if ( ((mykpart==ktota) .and. (frag))
     & .or. (mykpart==kfrag) ) then
         part=kfrag
         scale=myscale
         facscale=myfacscale
         if (mykpart==kfrag) then
         scalereset=.true.
         else
         reset=.true.
         endif
         ncall=myncall
c         write(6,*) 'Adjusting number of points for frag to',ncall
         ndim=ndim+1
!         call boundregion(ndim,region) 
         fragint_mode=.true.            ! for isolation
!         call vegasnr(region,ndim,fragint,myinit,myncall,myitmx,
!     &                 nprn,sigfrag,sdfrag,chifrag)
         call vegas_Pomp(fragint,sig,0._dp,0._dp,ndim,myncall,myitmx,myinit)
         fragint_mode=.false.            ! for isolation
         ndim=ndim-1
       rescale=.false.       ! turn rescaling off again
      endif

   33 continue 

c--- special process: 5FNS + 4FNS: initialization
c---   for process 421 (426) we will actually run
c---   process 411 (416) [5FNS] and then add on the result of process 401 (406) [4FNS]
      if ((mynproc == 421) .or. (mynproc == 426)) then 
        sumsig=sumsig+sig
        sumsigr=sumsigr+sigr
        sumsd=sums.e+_dpsd**2
        sumsdr=sumsdr+sdr**2
        if (nproc >= 411) then
          nproc=nproc-10
        goto 77
      endif
c--- return nproc to the value from the input file
        nproc=mynproc       
      sig=sumsig
      sigr=sumsigr
      sd=sqrt(sumsd)
      sdr=sqrt(sumsdr)
      endif
      
c--- special handling of multiple PS generation for
c--- Z+gamma+jet and Z+gamma+gamma processes: we will loop
c--- over ipsgen=1,..,n where n=2 for Z+gamma+jet and n=4 for  Z+gamma+gamma
      if ( (mynproc == 301) .or. (mynproc == 302)
     & .or.(mynproc == 303) .or. (mynproc == 304)
     & .or.(mynproc == 370) .or. (mynproc == 371)
     & .or.(doipsgen)) then
c--- store cross section and standard deviation info for this region on first run
        if (myinit == 0) then
          sigips(ipsgen)=sig+sigr+sigfrag
          sdips(ipsgen)=sqrt(sd**2+sdr**2+sdfrag**2)
        endif
        sumsig=sumsig+sig
        sumsigr=sumsigr+sigr
        sumsigf=sumsigf+sigfrag
        sumsd=sums.e+_dpsd**2
        sumsdr=sumsdr+sdr**2
        sumsdf=sumsdf+sdfrag**2
        ipsgen=ipsgen+1
        if ((((nproc == 302).or.(nproc == 304)).and.(ipsgen <= 2))
     &  .or.(((nproc == 301).or.(nproc == 303)).and.(ipsgen <= 4))
     &  .or.(((nproc == 370).or.(nproc == 371)).and.(ipsgen <= 4))
     &  .or.((doipsgen).and.(ipsgen <= maxipsgen))
     &     ) then
c--- W+2 photons case
          if ((nproc == 370).or.(nproc == 371)) then
            if ((myinit == 1) .and. (ipsgen > 1)) then
c--- dynamic scaling of points on second run according to uncertainties found in first
              xcallwt=(sdips(ipsgen)/sdips(1))
              if (xcallwt < 0.2_dp) xcallwt=0.2_dp
              if (xcallwt > 5._dp) xcallwt=5._dp
              myncall=int(float(myncall_save)*xcallwt)
c              write(6,*) '>>>>>> Dynamic scaling factor = ',
c     &                    xcallwt,' for ipsgen=',ipsgen
            endif
            goto 77 
          endif
c--- change number of PS points in different regions
          if ((nproc==301).or.(nproc==303)
     &    .or.(nproc==370).or.(nproc==371)) then
             if (ipsgen==2) then
                myncall=10*myncall
             elseif (ipsgen==3) then
                myncall=myncall/10
                myncall=5*myncall
             elseif (ipsgen==4) then
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
      sd=sqrt(sumsd)
      sdr=sqrt(sumsdr)
      sdfrag=sqrt(sumsdf)
      endif
      
c--- calculate integration variables to be returned
      xinteg=sig+sigr+sigdk+sigWdk+sigfrag
      xerr=sqrt(sd**2+sdr**2+sddk**2+sdWdk**2+sdfrag**2)      
      
c--- return part, scale and myncall to their real values
      part=mykpart
      scale=myscale
      myncall=myncall_save
      first=.false.
      
      return
      end
      
      
      subroutine boundregion(idim,region)
      implicit none
      include 'types.f'
c--- Initializes integration region [0,1] for each variable
c--- in the idim-dimensional integration range
      
      include 'mxdim.f'
      integer:: i,idim
      real(dp):: region(2*mxdim)
      
      do i=1,idim
      region(i)=0._dp
      region(i+idim)=1._dp
      enddo
      
      return
      end
      
      
