      subroutine chooser
      implicit none
      include 'types.f'
c---- Note added 4/21/03
c---- plabel set to 'ig' (for 'ignore') means that this
c---- particle should not be subject to any cuts, so that the
c---- total cross-section comes out correctly when the BR is removed
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'vegas_common.f'
      include 'zerowidth.f'
      include 'removebr.f'
      include 'bbproc.f'
      include 'nwz.f'
      include 'kprocess.f'
      include 'flags.f'
      include 'heavyflav.f'
      include 'nflav.f'
      include 'nodecay.f'
      include 'stopscales.f'
      include 'scale.f'
      include 'facscale.f'
      include 'nlooprun.f'
      include 'b0.f'
      include 'colstruc.f'
      include 'stopbmass.f'
      include 'fourthgen.f'
      include 'anomcoup.f'
      include 'srdiags.f'
      include 'clustering.f'
      include 'frag.f'
      include 'plabel.f'
      include 'interference.f'
      include 'couple.f'
      include 'kpart.f'
      include 'hdecaymode.f'
      include 'breit.f'
      include 'mcfmplotinfo.f'
      include 'lastphot.f'
      include 'new_pspace.f'
      include 'dm_params.f'
      include 'swapxz.f'
      include 'verbose.f'
      include 'runstring.f'
      include 'vdecayid.f'
      include 'ipsgen.f'
      include 'nuflav.f'
      include 'nqcdjets.f'
      include 'bitflags.f'
      include 'nproc.f'
      include 'notag.f'
      include 'VVstrong.f'
      include 'taucut.f'
      include 'mpicommon.f'
      include 'noglue.f'
      include 'toploopgaga.f'
      real(dp):: wwbr,zzbr,tautaubr,gamgambr,zgambr,Rcut,Rbbmin,
     & alphas,cmass,bmass
      real(dp):: br,BrnRat,brwen,brzee,brznn,brtau,brtop,brcharm
      integer:: mproc,j,isub,ilomomenta
      logical :: reset_alphaEM
      character*100 pname
      character*1 order
      character*72 string
      character*82 pwrite
      real(dp):: Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/Rbbmin/Rbbmin
      common/Rcut/Rcut
      common/BrnRat/BrnRat
      common/isub/isub
      common/ilomomenta/ilomomenta
      common/qmass/cmass,bmass
      data hdecaymode/'xxxx'/
     
      do j=1,mxpart
      plabel(j)=''
      enddo

      do j=1,50
      mcfmplotinfo(j)=-1
      enddo

      string='process.DAT' 
      open(unit=21,file=string,status='old',err=43)
      call checkversion(21,string)
      
      if (verbose) write(6,*) 'Chooser:process chosen by nproc=',nproc

      do j=1,600
      read(21,*,err=44) mproc,pname,order
      
      if (rank.eq.0) then
         if (nproc < 0) then 
            write(6,*) mproc,pname 
         endif
      endif
      if (mproc == nproc) go to 42
      if (pname == 'EOF') go to 44
      enddo
      goto 44

 42   continue
      if (rank.eq.0) then
         if (verbose) then
       write(6,*)
       write(6,*) '*************************** f(p1)+f(p2) --> *****'//
     .           '*************************************'
      j=index(pname,'[')
      if (j > 101) j=101
      if (j > 0) then
        pwrite=adjustl(pname(19:j-1))
        write(6,99) pwrite
        pwrite=adjustl(pname(j:len(pname)))
        write(6,99) pwrite
      else
        pwrite=adjustl(pname(19:100))
        write(6,99) pwrite
      endif
      write(6,*) '*************************************************'//
     & '*************************************'
      write(6,*)
      endif
      endif
      close(unit=21)

c--- check no. of momenta appearing in LO process, fill ilomomenta common block
      if (index(pname,'p3') > 0)  ilomomenta=3
      if (index(pname,'p4') > 0)  ilomomenta=4
      if (index(pname,'p5') > 0)  ilomomenta=5
      if (index(pname,'p6') > 0)  ilomomenta=6
      if (index(pname,'p7') > 0)  ilomomenta=7
      if (index(pname,'p8') > 0)  ilomomenta=8
      if (index(pname,'p9') > 0)  ilomomenta=9
      if (index(pname,'p10') > 0) ilomomenta=10
      if (index(pname,'p11') > 0) ilomomenta=11
      if (index(pname,'p12') > 0) ilomomenta=12
      
      plabel(1)='pp'
      plabel(2)='pp'

c--- the default behaviour is to remove no branching ratio
      BrnRat=1._dp

c--- set up most parameters
      call coupling
      
      notag=0
      nqcdjets=0
      isub=0
      bbproc=.false.
      nodecay=.false.
      nfonly=.false.
      caonly=.false.
      fourthgen=.false.
      rescale=.false.
      doipsgen=.false.
      reset_alphaEM=.false.
      lastphot=-1
c--- default is no interference contributions from identical fermions
      interference=.false.
      vsymfact=1._dp
c--- default is weak process
      VVstrong=.false.

c-- Rbbmin is an additional variable, added so that the separation
c-- between two b jets can be controlled separately from the Delta_R
c-- cut between any other types of jet
c-- Default behaviour: the same value as for the other jets
      Rbbmin=Rcut
      
c-----------------------------------------------------------------------

      if ((nproc == 1) .or. (nproc == 6)) then
        kcase=kW_only
        mass3=wmass
        width3=wwidth
        n3=1
        ndim=4
        nqcdjets=0

        mcfmplotinfo= (/ 34, (0,j=1,49) /)

c---W^+
        if     (nproc == 1) then
C-- 1  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))'
C--    '  f(p1)+f(p2) --> W^+ (for total Xsect)' (removebr=.true.)
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='pp'
          nwz=1
c---W^-
        elseif (nproc == 6) then
c-- 6  '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))'
c--    '  f(p1)+f(p2) --> W^- (for total Xsect)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='na'
          plabel(5)='pp'
          nwz=-1
        else
          call nprocinvalid()
        endif

c--- simple hack for W-prime cross-section
        if (runstring(1:6) == 'wprime') then
          read(runstring(7:10),'(f4.1)') mass3
          mass3=mass3*1d3
          width3=wwidth*(mass3/wmass) ! scale W width linearly with W' mass
c--- width fit from Fig. 10 of Altarelli et al, Z. Phys. C45, 109 (1989)
          width3=4._dp+(mass3-100._dp)/30._dp
          wmass=mass3
          wwidth=width3
c--- sum over both W+ and W- 
          nwz=2
          write(6,*) '<<<<<< WARNING <<<<<'
          write(6,*) ' computing W-prime cross section'
          write(6,*) ' summed over both +ve and -ve charges'
          write(6,*) ' for W-prime of mass, width (GeV): ',wmass,wwidth
        endif
        
c--- total cross-section
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif

c-----------------------------------------------------------------------

      elseif ((nproc == 11) .or. (nproc == 16)) then
        kcase=kW_1jet
        nqcdjets=1
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        plabel(5)='pp'
        plabel(6)='pp'

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 11) then
C-- 11 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+f(p5)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) + f(p5)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 16) then
c-- 16 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+f(p5)'
c--    '  f(p1)+f(p2) --> W^- (no BR) + f(p5)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc == 12) .or. (nproc == 17)) then
        kcase=kWbfrmc
        nqcdjets=1
        ndim=7
        n2=0
        n3=1
        nflav=4
        plabel(5)='bq'
        plabel(6)='pp'

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 12) then
c-- 12 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+bbar(p5)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) + bbar(p5)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 17) then
c-- 17 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+b(p5)'
c--    '  f(p1)+f(p2) --> W^- (no BR) + b(p5)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
        mass3=wmass
        width3=wwidth
        mass2=mb
        if (Vcb == zip) Vcb=0.041_dp
        write(6,*) 'Setting Vcb=0.041 for this process'
        if (Vub == zip) Vub=0.00347_dp
        write(6,*) 'Setting Vub=0.00347 for this process'

c-----------------------------------------------------------------------

      elseif ((nproc == 13) .or. (nproc == 18)) then
        kcase=kW_cjet
c--- this process works best using the new PS generation
        new_pspace=.true.
c--- this process works best using the new PS generation
        nqcdjets=1
        ndim=7
        mb=0
        n2=0
        n3=1
        nflav=3
        plabel(5)='bq'
        plabel(6)='pp'

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 13) then
c-- 13 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+cbar(p5)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) + cbar(p5)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 18) then
c-- 18 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+c(p5)'
c--    '  f(p1)+f(p2) --> W^- (no BR) + c(p5)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c--- change W mass (after couplings and BRs already calculated)
c      if     (runstring(4:8) == 'mw_80') then
c        wmass=80.4_dp
c      elseif (runstring(4:8) == 'mw200') then
c        wmass=200._dp
c      elseif (runstring(4:8) == 'mw400') then
c        wmass=400._dp
c      endif
      
c--- change charm mass
c      if     (runstring(9:13) == 'mc1.3') then
c        mc=1.3_dp
c        mcsq=mc**2
c      elseif (runstring(9:13) == 'mc4.5') then
c        mc=4.5_dp
c        mcsq=mc**2
c      elseif (runstring(9:13) == 'mc20.') then
c        mc=20._dp
c        mcsq=mc**2
c      endif
      
        mass3=wmass
        width3=wwidth
        mass2=mc

c-----------------------------------------------------------------------

      elseif ((nproc == 14) .or. (nproc == 19)) then
        kcase=kWcjet0
        nqcdjets=1
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        mass2=zip
        nflav=3
        plabel(5)='bq'
        plabel(6)='pp'

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 14) then
c-- 13 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+cbar(p5) [massless]'
c--    '  f(p1)+f(p2) --> W^+ (no BR) + cbar(p5) [massless]' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 19) then
c-- 18 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+c(p5) [massless]'
c--    '  f(p1)+f(p2) --> W^- (no BR) + c(p5) [massless]' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc == 20) .or. (nproc == 25)) then
        kcase=kWbbmas
        swapxz=.true.
        write(6,*) 'mb=',mb
        nqcdjets=2
        flav=5
        bbproc=.true.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
 
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 20) then
c-- 20 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+b~(p6) [massive]'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +b(p5)+b~(p6) [massive]' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 25) then
c-- 25 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + b(p5)+b~(p6) [massive]'
c--    '  f(p1)+f(p2) --> W^- (no BR) +b(p5)+b~(p6) [massive]' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
 
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc == 21) .or. (nproc == 26)) then
        kcase=kWbbbar
        write(6,*) 'mb=',mb
        nqcdjets=2
        bbproc=.true.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        mb=0
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 21) then
c-- 21 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+b~(p6)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +b(p5)+b~(p6)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 26) then
c-- 26 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + b(p5)+b~(p6)'
c--    '  f(p1)+f(p2) --> W^- (no BR) +b(p5)+b~(p6)' (removebr=.true.)
          nwz=-1 
          plabel(3)='el'
          plabel(4)='na'
        endif
             
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc == 22) .or. (nproc == 27)) then
        kcase=kW_2jet
        nqcdjets=2
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 22) then
c-- 22 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+f(p5)+f(p6)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +f(p5)+f(p6)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 27) then
c-- 27 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + f(p5)+f(p6)'
c--    '  f(p1)+f(p2) --> W^- (no BR) +f(p5)+f(p6)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc == 23) .or. (nproc == 28)) then
        kcase=kW_3jet
        nqcdjets=3
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=13
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 23) then
c-- 23 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+f(p5)+f(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +f(p5)+f(p6)+f(p7)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 28) then
c-- 28 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + f(p5)+f(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> W^- (no BR) +f(p5)+f(p6)+f(p7)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc == 24) .or. (nproc == 29)) then
        kcase=kWbbjet
        write(6,*) 'mb=',mb
        nqcdjets=3
        bbproc=.true.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        mb=zip
        ndim=13
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 24) then
c-- 24 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+b~(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> W^+ (no BR) +b(p5)+b~(p6)+f(p7)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 29) then
c-- 29 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + b(p5)+b~(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> W^- (no BR) +b(p5)+b~(p6)+f(p7)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc > 30) .and. (nproc <= 35)) then
        kcase=kZ_only
        nqcdjets=0
        nwz=0
        mass3=zmass
        width3=zwidth
        n3=1
        ndim=4
        plabel(5)='pp'

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 31) then
c-- 31 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))'
c--    '  f(p1)+f(p2) --> Z^0 (for total Xsect)' (removebr=.true.)
          call checkminzmass(1)
          plabel(3)='el'
          plabel(4)='ea'
          q1=-1._dp
          l1=le
          r1=re
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee
          endif
        elseif (nproc == 32) then
c-- 32 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4)))'
          plabel(3)='nl'
          plabel(4)='na'
          q1=zip
          l1=ln*sqrt(3._dp)
          r1=rn*sqrt(3._dp)
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brznn
          endif
        elseif (nproc == 33) then
c-- 33 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4))'
          call checkminzmass(1)
          plabel(3)='qb'
          plabel(4)='ab'
          q1=Q(5)*sqrt(xn)
          l1=l(5)*sqrt(xn)
          r1=r(5)*sqrt(xn)
        elseif (nproc == 34) then
c-- 34 '  f(p1)+f(p2) --> Z^0(-->3*(d(p3)+d~(p4)))'
          call checkminzmass(1)
          plabel(3)='qb'
          plabel(4)='ab'
          q1=Q(1)*sqrt(3._dp*xn)
          l1=l(1)*sqrt(3._dp*xn)
          r1=r(1)*sqrt(3._dp*xn)
        elseif (nproc == 35) then
c-- 35 '  f(p1)+f(p2) --> Z^0(-->2*(u(p3)+u~(p4)))'
          call checkminzmass(1)
          plabel(3)='qb'
          plabel(4)='ab'
          q1=Q(2)*sqrt(2._dp*xn)
          l1=l(2)*sqrt(2._dp*xn)
          r1=r(2)*sqrt(2._dp*xn)
        else
          call nprocinvalid()
        endif 

c-----------------------------------------------------------------------

        elseif (nproc == 36) then
        kcase=kttZbbl
        nwz=1
        ndim=16
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth
        bbproc=.true.
        
        mcfmplotinfo= (/ 34, 78, 345, 678, (0,j=1,46) /)
        
c--  36 '  f(p1)+f(p2) -> Z -> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+e^-(p7)+nu~(p8)'
          nqcdjets=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=(brwen*brtop)**2
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            plabel(7)='ig'
            plabel(8)='ig'
            nqcdjets=0
            bbproc=.false.
          endif


c-----------------------------------------------------------------------
      
      elseif ((nproc >= 41) .and. (nproc <= 43)) then
        kcase=kZ_1jet
        nqcdjets=1
        nwz=0
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
        plabel(5)='pp'
        plabel(6)='pp'

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 41) then
c-- 41 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+f(p5)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +f(p5)' (removebr=.true.)
          call checkminzmass(1)
          plabel(3)='el'
          plabel(4)='ea'
          q1=-1._dp
          l1=le
          r1=re
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee
          endif
        elseif (nproc == 42) then
c-- 42 '  f(p1)+f(p2) --> Z_0(-->3*(nu(p3)+nu~(p4)))-(sum over 3 nu)+f(p5)'
            plabel(3)='nl'
            plabel(4)='na'
            q1=zip
            l1=ln*sqrt(3._dp)
            r1=rn*sqrt(3._dp)
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
              BrnRat=brznn
            endif
        elseif (nproc == 43) then
c-- 43 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4))+f(p5)'
          call checkminzmass(1)
          plabel(3)='qb'
          plabel(4)='ab'
          q1=Q(5)*sqrt(xn)
          l1=l(5)*sqrt(xn)
          r1=r(5)*sqrt(xn)
        endif
      
c-----------------------------------------------------------------------
      
      elseif (nproc == 44) then
c-- 44 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+f(p5)+f(p6)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +f(p5)+f(p6)' (removebr=.true.)
        kcase=kZ_2jet
        call checkminzmass(1)
        ndim=10
        n2=0
        n3=1
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        q1=-1._dp
        l1=le
        r1=re
        nwz=0   
        mass3=zmass
        width3=zwidth
       
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
        endif
      
      elseif (nproc == 45) then
c-- 45 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+f(p5)+f(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +f(p5)+f(p6)+f(p7)' (removebr=.true.)
        kcase=kZ_3jet
        call checkminzmass(1)
        ndim=13
        n2=0
        n3=1
        nqcdjets=3
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        q1=-1._dp
        l1=le
        r1=re
        nwz=0   
        mass3=zmass
        width3=zwidth
       
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
        endif

c-----------------------------------------------------------------------
      
      elseif (nproc == 46) then
c-- 46 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))+f(p5)+f(p6)'
        kcase=kZ_2jet
        ndim=10
        n2=0
        n3=1
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='na'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        q1=zip
        l1=ln*sqrt(3._dp)
        r1=rn*sqrt(3._dp)
        nwz=0   
        mass3=zmass
        width3=zwidth
       
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brznn
        endif

      elseif (nproc == 47) then
c-- 47 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))+f(p5)+f(p6)+f(p7)'
        kcase=kZ_3jet
        ndim=13
        n2=0
        n3=1
        nqcdjets=3
        plabel(3)='nl'
        plabel(4)='na'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        q1=zip
        l1=ln*sqrt(3._dp)
        r1=rn*sqrt(3._dp)
        nwz=0   
        mass3=zmass
        width3=zwidth
       
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brznn
        endif

c-----------------------------------------------------------------------
          
      elseif (nproc == 50) then
c-- 50 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b~(p5)+b(p6) (massive)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +b~(p5)+b(p6) (massive)' (removebr=.true.)
        kcase=kZbbmas
        call checkminzmass(1)
        write(6,*) 'mb=',mb
        bbproc=.true.
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        q1=-1._dp
        l1=le
        r1=re
        ndim=10
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
      
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        flav=5
        nflav=4

c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
        endif

      elseif ((nproc >= 51) .and. (nproc <= 53)
     &   .or. (nproc >= 56) .and. (nproc <= 58)) then
        kcase=kZbbbar
        call checkminzmass(1)
        bbproc=.true.
        mb=zip
        nqcdjets=2
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        q1=-1._dp
        l1=le
        r1=re
        ndim=10
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
        if     (nproc <= 53) then
          flav=5
          nflav=4
        else
          flav=4
          nflav=3
        endif

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     ((nproc == 51) .or. (nproc == 56)) then
c-- 51 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+b~(p6)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +b(p5)+b~(p6)' (removebr=.true.)
c-- 56 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c(p5)+c~(p6)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +c(p5)+c~(p6)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='ea'
          q1=-1._dp
          l1=le
          r1=re
c--- total cross-section             
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee
          endif
        elseif (nproc == 52) then
c-- 52 '  f(p1)+f(p2) --> Z_0(-->3*(nu(p3)+nu~(p4)))+b(p5)+b~(p6)'
          plabel(3)='nl'
          plabel(4)='na'
          q1=zip
          l1=ln*sqrt(3._dp)
          r1=rn*sqrt(3._dp)
        elseif (nproc == 53) then
c-- 53 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4))+b(p5)+b~(p6)'
          plabel(3)='qb'
          plabel(4)='ab'
          q1=Q(5)*sqrt(xn)
          l1=l(5)*sqrt(xn)
          r1=r(5)*sqrt(xn)
        endif
        
      elseif (nproc == 54) then
c-- 54 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+b~(p6)+f(p7)'
c--    '  f(p1)+f(p2) --> Z^0 (no BR) +b(p5)+b~(p6)+f(p7)' (removebr=.true.)
        kcase=kZbbjet
        ndim=13
        bbproc=.true.
        mb=zip
        nqcdjets=3
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        q1=-1._dp
        l1=le
        r1=re
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
        flav=5
        nflav=4

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
        endif

c-----------------------------------------------------------------------
          
      elseif (nproc/10 == 6) then
        kcase=kWWqqbr
        call readcoup
        nqcdjets=0  
        plabel(7)='pp'
        nwz=1
        ndim=10
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

        mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)
        
c--- include singly resonant diagrams if zerowidth=.false. , but only
c---  as long as anomtgc=.false. too
        srdiags=((zerowidth .eqv. .false.)
     &     .and. ( anomtgc  .eqv. .false.))
      
c--- zero srdiags for CDFdijet calculation 
c        if(runstring(1:10)=='cdf_Wdijet') then 
c           srdiags=.false.
c        endif

        if    ((nproc == 61) .or. (nproc == 69)) then
c--  61 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6))'
c--     '  f(p1)+f(p2) --> W^+ + W^- (for total Xsect)' (removebr=.true.)
          if (nproc == 69) then
c--  69 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6)) [no pol]'
            kcase=kWWnpol
          write(*,*)'Setting zerowidth to true for process 69'
          zerowidth = .true.
          endif
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='el'
          plabel(6)='na'
          l1=1._dp
c--- total cross-section             
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen**2
          endif
        elseif (nproc == 62) then
c--  62 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->q(p5)+q~(p6))'
c--- note: scattering diagrams are NOT included, only couplings change
          kcase=kWWqqbr
          nqcdjets=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='qj'
          plabel(6)='qj'
          plabel(7)='pp'
          l1=sqrt(xn*2._dp)

        elseif (nproc == 63) then 
c--  63 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->q(p5)+q~(p6)) [rad.in.dk]'
          kcase=kWWqqdk
          nqcdjets=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='qj'
          plabel(6)='qj'
          plabel(7)='pp'
          l1=sqrt(xn*2._dp)

        elseif (nproc == 64) then
c--  64 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+W^+(--> q(p5)+ q~(p6))'
c--- note: scattering diagrams are NOT included, only couplings change
          kcase=kWWqqbr
          nqcdjets=2
          plabel(5)='qj'
          plabel(6)='qj'
          plabel(3)='el'
          plabel(4)='na'
          plabel(7)='pp'
          l1=sqrt(xn*2._dp)
        elseif (nproc == 65) then
c--  65 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+W^+(--> q(p5)+ q~(p6)),[rad.in.dk]'
c--- note: scattering diagrams are NOT included, only couplings change
          kcase=kWWqqdk
          nqcdjets=2
          plabel(5)='qj'
          plabel(6)='qj'
          plabel(3)='el'
          plabel(4)='na'
          plabel(7)='pp'
          l1=sqrt(xn*2._dp)
        elseif (nproc == 66) then
c--  66 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+W^-(-->e^-(p5)+nu~(p6))+f(p7)'
c--     '  f(p1)+f(p2) --> W^+ + W^- + f(p7) (for total Xsect)' (removebr=.true.)
          kcase=kWW_jet
          nflav=4
          nqcdjets=1
          ndim=13
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='el'
          plabel(6)='na'
          plabel(7)='pp'
          plabel(8)='pp'
          l1=1._dp
c--- total cross-section             
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen**2
          endif
        endif

c-----------------------------------------------------------------------

      elseif ((nproc >= 70) .and. (nproc <= 80)) then
        kcase=kWZbbar
        call checkminzmass(2)
        call readcoup
        plabel(7)='pp'
        ndim=10
        mb=0
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)
        
c--- include singly resonant diagrams if zerowidth=.false. , but only
c---  as long as anomtgc=.false. too
        srdiags=((zerowidth .eqv. .false.)
     &     .and. ( anomtgc  .eqv. .false.))
      
c--- Zero srdiags for CDFdijet calculation 
c        if(runstring(1:10)=='cdf_Wdijet') then 
c           srdiags=.false.
c        endif

        if (nproc <= 75) then
c-- W^+Z
          nwz=+1

          if     (nproc == 71) then             
c--  71 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->e^-(p5)+e^+(p6))'
c--     '  f(p1)+f(p2) --> W^+ (for total Xsect) + Z^0 ' (removebr=.true.)
            plabel(3)='nl'
            plabel(4)='ea'
            plabel(5)='ml'
            plabel(6)='ma'
            q1=-1._dp
            l1=le
            r1=re
c--- total cross-section             
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              plabel(5)='ig'
              plabel(6)='ig'
              call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
              BrnRat=brwen*brzee
            endif
          elseif (nproc == 72) then
c--  72 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->nu_e(p5)+nu~_e(p6))'
            plabel(3)='nl'
            plabel(4)='ea'
            plabel(5)='nl'
            plabel(6)='na'
            plabel(7)='pp'
            q1=zip
            l1=ln*sqrt(3._dp)
            r1=rn*sqrt(3._dp)
          elseif (nproc == 73) then
c--  73 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->b(p5)+b~(p6))'
            bbproc=.true.
            nqcdjets=2
            plabel(3)='nl'
            plabel(4)='ea'
            plabel(5)='bq'
            plabel(6)='ba'
            plabel(7)='pp'
            q1=Q(5)*sqrt(xn)
            l1=l(5)*sqrt(xn)
            r1=r(5)*sqrt(xn)
          elseif (nproc == 74) then
c--  74 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->3*(d(p5)+d~(p6)))'
            nqcdjets=2
            plabel(3)='nl'
            plabel(4)='ea'
            plabel(5)='qj'
            plabel(6)='qj'
            plabel(7)='pp'
            q1=Q(5)*sqrt(3._dp*xn)
            l1=l(5)*sqrt(3._dp*xn)
            r1=r(5)*sqrt(3._dp*xn)
          elseif (nproc == 75) then
c--  75 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->2*(u(p5)+u~(p6)))'
            nqcdjets=2
            plabel(3)='nl'
            plabel(4)='ea'
            plabel(5)='qj'
            plabel(6)='qj'
            plabel(7)='pp'
            q1=Q(4)*sqrt(2._dp*xn)
            l1=l(4)*sqrt(2._dp*xn)
            r1=r(4)*sqrt(2._dp*xn)
          else
            call nprocinvalid()
          endif 

        elseif (nproc >= 76) then
c-- W^-Z
          nwz=-1

          if     (nproc == 76) then
c--  76 '  f(p1)+f(p2) --> W^-(-->mu^-(p3)+nu~(p4))+Z^0(-->e^-(p5)+e^+(p6))'
c--     '  f(p1)+f(p2) --> W^- + Z^0 (for total Xsect)' (removebr=.true.)
            plabel(3)='el'
            plabel(4)='na'
            plabel(5)='ml'
            plabel(6)='ma'
            plabel(7)='pp'
            q1=-1._dp
            l1=le
            r1=re
c--- total cross-section             
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              plabel(5)='ig'
              plabel(6)='ig'
              call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
              BrnRat=brwen*brzee
            endif
          elseif (nproc == 77) then
c--  77 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->nu(p5)+nu~(p6))'
            plabel(3)='el'
            plabel(4)='na'
            plabel(5)='nl'
            plabel(6)='na'
            plabel(7)='pp'
            q1=zip
            l1=ln*sqrt(3._dp)
            r1=rn*sqrt(3._dp)
          elseif (nproc == 78) then
c--  78 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->b(p5)+b~(p6))'
            bbproc=.true.
            nqcdjets=2
            plabel(3)='el'
            plabel(4)='na'
            plabel(5)='bq'
            plabel(6)='ba'
            plabel(7)='pp'
            q1=Q(5)*sqrt(xn)
            l1=l(5)*sqrt(xn)
            r1=r(5)*sqrt(xn)
          elseif (nproc == 79) then
c--  79 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->3*(d(p5)+d~(p6))'
            nqcdjets=0
            plabel(3)='el'
            plabel(4)='na'
            plabel(5)='qj'
            plabel(6)='qj'
            plabel(7)='pp'
            q1=Q(5)*sqrt(3._dp*xn)
            l1=l(5)*sqrt(3._dp*xn)
            r1=r(5)*sqrt(3._dp*xn)
          elseif (nproc == 80) then
c--  80 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+Z^0(-->2*(u(p5)+u~(p6)))'
            nqcdjets=0
            plabel(3)='el'
            plabel(4)='na'
            plabel(5)='qj'
            plabel(6)='qj' 
            plabel(7)='pp'
            q1=Q(4)*sqrt(2._dp*xn)
            l1=l(4)*sqrt(2._dp*xn)
            r1=r(4)*sqrt(2._dp*xn)
          else
            call nprocinvalid()
          endif 

        endif
            
c-----------------------------------------------------------------------

      elseif ((nproc > 80) .and. (nproc <= 90)) then
        kcase=kZZlept
        call checkminzmass(1)
        if ((nproc == 81) .or. (nproc == 83)
     &  .or.(nproc == 90)) call checkminzmass(2)
        call readcoup
        plabel(7)='pp'
        nqcdjets=0
        nwz=0
        ndim=10
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=zmass
        width3=zwidth
        q1=-1._dp
        l1=le
        r1=re
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)
        
c--- only include singly-resonant diagrams when not in zerowidth approx.        
      if (zerowidth) then
        srdiags=.false.
      else
        srdiags=.true.
      endif
      
        if (nproc == 81 .or. nproc == 86) then
c--  81 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->mu^-(p5)+mu^+(p6))'
c--     '  f(p1)+f(p2) --> Z^0 + Z^0 (for total Xsect)' (removebr=.true.)
c--  86 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6))+Z^0(-->mu^-(p3)+mu^+(p4)) (NO GAMMA*)'
c--     '  f(p1)+f(p2) --> Z^0 + Z^0 (for total Xsect) (NO GAMMA*)' (removebr=.true.)
          plabel(3)='el' 
          plabel(4)='ea'
          plabel(5)='ml'
          plabel(6)='ma'
          q2=-1._dp
          l2=le
          r2=re

c--- check runstring to change from (e,mu) Z decays to (e,e) or (mu,mu)
c          if     (index(runstring,'ELEL') > 0) then
c            plabel(5)='el'
c            plabel(6)='ea'
c            interference=.true.
c            vsymfact=0.25_dp
c          elseif (index(runstring,'MUMU') > 0) then
c            plabel(3)='ml'
c            plabel(4)='ma'
c            interference=.true.
c            vsymfact=0.25_dp
c          endif

c--- if vector boson decays specified, initialize appropriately
          if (vdecayid) then
            call setvdecay(34,0)
            call setvdecay(56,0)
            if (plabel(3) == plabel(5)) then
              if (plabel(3) .ne. 'nl') then
c------ for both Z decays to neutrinos, neglect interference effects for simplicity
                interference=.true.
                vsymfact=0.25_dp
              endif
            endif
          endif
          

          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2._dp*brzee**2  ! factor of 2 for identical particles
          endif
        elseif (nproc == 82 .or. nproc == 87) then
c--  82 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->3*(nu(p5)+nu~(p6)))'
c--  87 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->3*(nu(p5)+nu~(p6))) [no gamma^*]'
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='nl'
          plabel(6)='na'
          q2=zip
          l2=ln*sqrt(3._dp)
          r2=rn*sqrt(3._dp)
        elseif (nproc == 83 .or. nproc == 88) then
c--  83 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->b(p5)+b~(p6))'
c--  88 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->b(p5)+b~(p6)) [no gamma^*]'
          mb=0
          bbproc=.true.
          nqcdjets=2
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          q2=Q(5)*sqrt(xn)
          l2=l(5)*sqrt(xn)
          r2=r(5)*sqrt(xn)
        elseif ((nproc == 84) .or. (nproc == 89))  then
c--  84 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4)) + Z^0(-->3*(nu(p5)+nu~(p6)))'
c--  89 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4)) + Z^0(-->3*(nu(p5)+nu~(p6))) [no gamma^*]'
          mb=0
          bbproc=.true.
          nqcdjets=2
          plabel(3)='bq'
          plabel(4)='ba'
          plabel(5)='nl'
          plabel(6)='na'
          q2=zip
          l2=ln*sqrt(3._dp)
          r2=rn*sqrt(3._dp)
          q1=Q(5)*sqrt(xn)
          l1=l(5)*sqrt(xn)
          r1=r(5)*sqrt(xn)
      elseif (nproc == 85) then
c---  85 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->3*(nu(p5)+nu~(p6)))+f(p7)'
        kcase=kZZ_jet
        nqcdjets=1
        ndim=13
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='nl'
          plabel(6)='na'
          plabel(7)='pp'
          q2=zip
          l2=ln*sqrt(3._dp)
          r2=rn*sqrt(3._dp)
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2._dp*brzee*brznn  ! factor of 2 for identical particles
          endif
        elseif (nproc == 90) then
c--  90 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))'
          interference=.true.
          vsymfact=0.25_dp
          plabel(3)='el' 
          plabel(4)='ea'
          plabel(5)='el'
          plabel(6)='ea'
          q2=-1._dp
          l2=le
          r2=re
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2._dp*brzee**2  ! factor of 2 for identical particles
          endif
        
          mcfmplotinfo= (/ 34, 56, 36, 45, 3456, (0,j=1,45) /)

        else
          call nprocinvalid()
        endif 

c-- remove gamma^* if necessary
        if ((nproc > 85) .and. (nproc < 90)) then
          q1=zip
          q2=zip
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc == 91) .or. (nproc == 96)) then
        kcase=kWHbbar
        hdecaymode='tlta'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(5)='tl'
        plabel(6)='ta'
        plabel(7)='pp'
        
        ndim=10
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)
        
        if     (nproc == 91) then
c--  91 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) H(-->tau-(p5)+tau+(p6)) '
c--     '  f(p1)+f(p2) --> W+ + H (for total Xsect)' (removebr=.true.)
          plabel(3)='nl'
          plabel(4)='ea'
          nwz=1
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen*tautaubr
          endif
        elseif (nproc == 96) then
c--  96 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+ H(-->tau-(p5)+tau+(p6))' 
c--     '  f(p1)+f(p2) --> W- + H (for total Xsect)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='na'
          nwz=-1
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen*tautaubr
            bbproc=.false.
            nqcdjets=0
          endif
        else
            call nprocinvalid()
        endif

c-----------------------------------------------------------------------

      elseif ((nproc == 92) .or. (nproc == 97)
     &   .or. (nproc == 920) .or. (nproc == 970)) then
        kcase=kWHbbar
        if ((nproc == 920) .or. (nproc == 970)) kcase=kWHbbdk
        hdecaymode='bqba'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        nqcdjets=2
        notag=2
        
        ndim=10
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)
        
        if     ((nproc == 92) .or. (nproc == 920)) then
c--  92 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) H(-->b(p5)+b~(p6)) '
c--     '  f(p1)+f(p2) --> W+ + H (for total Xsect)' (removebr=.true.)
          plabel(3)='nl'
          plabel(4)='ea'
          nwz=1
c--- uncomment to sum over both W+ and W- 
c          nwz=2
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen*br
            nqcdjets=0
            notag=0
          endif
        elseif ((nproc == 97) .or. (nproc == 970)) then
c--  97 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+ H(-->b(p5)+b~(p6))' 
c--     '  f(p1)+f(p2) --> W- + H (for total Xsect)' (removebr=.true.)
          plabel(3)='el'
          plabel(4)='na'
          nwz=-1
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen*br
            nqcdjets=0
            notag=0
          endif
        else
            call nprocinvalid()
        endif

c-----------------------------------------------------------------------

      elseif ((nproc == 93) .or. (nproc == 98)) then
        kcase=kWHgaga
        mb=zip
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        hdecaymode='gaga'
        plabel(5)='ga'
        plabel(6)='ga'
        plabel(7)='pp'
        
        ndim=10
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)
        
        if     (nproc == 93) then
c '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + H(-->gamma(p5)+gamma(p6)' 'N'
          plabel(3)='nl'
          plabel(4)='ea'
          nwz=+1
        elseif (nproc == 98) then
c '  f(p1)+f(p2) --> W^-(-->nu(p3)+e^+(p4)) + H(-->gamma(p5)+gamma(p6)' 'N'
          plabel(3)='el'
          plabel(4)='na'
          nwz=-1
        endif

c--- uncomment to sum over both W+ and W- 
c        nwz=2

        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'             
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*gamgambr
        endif

c-----------------------------------------------------------------------

      elseif ((nproc == 94) .or. (nproc == 99)) then
        kcase=kWH__WW
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        hdecaymode='wpwm'
        nqcdjets=0
        plabel(5)='nl'
        plabel(6)='ea'
        plabel(7)='el'
        plabel(8)='na'
        plabel(9)='pp'
        
        mcfmplotinfo= (/ 34, 56, 78, 5678, (0,j=1,46) /)
        
        ndim=16
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=wmass
        width3=wwidth

        if     (nproc == 94) then
c--- 94 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + H(-->W^+(nu(p3),e^+(p4))W^-(e^-(p5),nub(p6)))'
          plabel(3)='nl'
          plabel(4)='ea'
          nwz=+1
        elseif (nproc == 99) then
          plabel(3)='el'
          plabel(4)='na'
          nwz=-1
        endif

        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          plabel(8)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen**3*wwbr
        endif

c--- print warning if we're below threshold
        if (hmass < 2._dp*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
      if (removebr) then
      write(6,*)
      write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
        stop
      endif
      if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif
                 
c-----------------------------------------------------------------------

        elseif     ((nproc == 95) .or. (nproc == 100)) then
C---95  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + H(-->Z(e^-(p5),e^+(p6))+Z(mu^-(p7),mu(p8)))' 
C---100 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + H(-->Z(e^-(p5),e^+(p6))+Z(mu^-(p7),mu(p8)))'
        kcase=kWH__ZZ
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=0
        plabel(5)='el'
        plabel(6)='ea'
        plabel(7)='el'
        plabel(8)='ea'
        plabel(9)='pp'

        ndim=16
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=wmass
        width3=wwidth

        l1=le
        r1=re
        l2=le
        r2=re

        mcfmplotinfo= (/ 34, 56, 78, 5678, (0,j=1,46) /)
        
        if     (nproc == 95) then
C---93  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + H(-->Z(e^-(p5),e^+(p6))+Z(mu^-(p7),mu(p8)))' 
          plabel(3)='nl'
          plabel(4)='ea'
          nwz=+1
        elseif (nproc == 100) then
C--98  '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + H(-->Z(e^-(p5),e^+(p6))+Z(mu^-(p7),mu(p8)))'
          plabel(3)='el'
          plabel(4)='na'
          nwz=-1
        endif

        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          plabel(8)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=2._dp*brwen*brzee**2*zzbr  ! factor of 2 for identical particles
        endif

c--- print warning if we're below threshold
        if (hmass < 2._dp*zmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
      if (removebr) then
      write(6,*)
      write(6,*) 'Cannot remove H->ZZ BR, not defined below threshold'
        stop
      endif
      if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif

c-----------------------------------------------------------------------

      elseif (((nproc >= 101) .and. (nproc <= 105))
     &    .or. (nproc==110) .or. (nproc==1010)) then
         kcase=kZHbbar
         if(nproc.eq.1010) kcase=kZHbbdk
        hdecaymode='bqba'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=2
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        nqcdjets=2
        notag=2

        ndim=10
        nwz=0 
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)
        
        if ((nproc == 101).or.(nproc==1010)) then
c--  101 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + H(-->b(p5)+b~(p6))'
c--      '  f(p1)+f(p2) --> H + Z0 (for total Xsect)' (removebr=.true.)
          hdecaymode='bqba'
          call checkminzmass(1)
          plabel(3)='el'
          plabel(4)='ea'

          q1=-1._dp
          l1=le
          r1=re
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee*br
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            nqcdjets=0
            notag=0
          endif
        elseif (nproc == 102) then
c--  102 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))) + H(-->b(p5)+b~(p6))'
          hdecaymode='bqba'
          plabel(3)='nl'
          plabel(4)='na'
          plabel(5)='bq'
          plabel(6)='ba'
          q1=zip
          l1=ln*sqrt(3._dp)
          r1=rn*sqrt(3._dp)
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brznn*br
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            nqcdjets=0
            notag=0
          endif
        elseif (nproc == 103) then
c--  103 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4)) + H(-->b(p5)+b~(p6))'     
          hdecaymode='bqba'
          call checkminzmass(1)
          nqcdjets=4
          plabel(3)='bq'
          plabel(4)='ba'
          plabel(5)='bq'
          plabel(6)='ba'
          q1=Q(5)*sqrt(xn)
          l1=l(5)*sqrt(xn)
          r1=r(5)*sqrt(xn)
        elseif (nproc == 104) then
c--  104 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + H(-->gamma(p5)+gamma(p6))' 'N'
          call checkminzmass(1)
          kcase=kZHgaga
          nqcdjets=0
          bbproc=.false.
          q1=-1._dp
          l1=le
          r1=re
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='ga'
          plabel(6)='ga'
          hdecaymode='gaga'
           if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee*gamgambr
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'    
          endif
        elseif (nproc == 105) then
c--  105 '  f(p1)+f(p2) --> Z^0(-->-->3*(nu(p3)+nu~(p4))) + H(-->gamma(p5)+gamma(p6))' 'N'
          kcase=kZHgaga
          nqcdjets=0
          bbproc=.false.
          plabel(3)='nl'
          plabel(4)='na'
          plabel(5)='ga'
          plabel(6)='ga'
          hdecaymode='gaga'
          q1=zip
          l1=ln*sqrt(3._dp)
          r1=rn*sqrt(3._dp)
           if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brznn*gamgambr
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'    
           endif
        elseif (nproc == 110) then
c--  104 '  f(p1)+f(p2) --> Z(-->e^-(p3)+e^+(p4)) + H(-->tau(p5)+tau~(p6))' 'N'
          call checkminzmass(1)
          kcase=kZHbbar
          nqcdjets=0
          bbproc=.false.
          q1=-1._dp
          l1=le
          r1=re
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='tl'
          plabel(6)='ta'
          hdecaymode='tlta'
           if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee*tautaubr
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'    
          endif
        else
          call nprocinvalid()
        endif 

c-----------------------------------------------------------------------

      elseif ((nproc >= 106) .and. (nproc <= 108)) then
        kcase=kZH__WW
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=0
        hdecaymode='wpwm'
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='nl'
        plabel(6)='ea'
        plabel(7)='el'
        plabel(8)='na'
        plabel(9)='pp'
        plabel(10)='pp'
        
        ndim=16
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 56, 78, 5678, (0,j=1,46) /)
        
        if (nproc == 106) then
c--  106 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + H(-->W^+(nu(p5),e^+(p6))W^-(e^-(p7),nub(p8)))'
c--      '  f(p1)+f(p2) --> H + Z0 (for total Xsect)' (removebr=.true.)
          call checkminzmass(1)
          plabel(3)='el'
          plabel(4)='ea'
          q1=0._dp
          l1=le
          r1=re
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'    
            plabel(7)='ig'
            plabel(8)='ig'             
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee*brwen**2*wwbr
          endif
        elseif (nproc == 107) then
c--  107 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))) + H(-->W^+(nu(p5),e^+(p6))W^-(e^-(p7),nub(p8)))'
          plabel(3)='nl'
          plabel(4)='na'
          q1=zip
          l1=ln*sqrt(3._dp)
          r1=rn*sqrt(3._dp)
        elseif (nproc == 108) then
c--  108 '  f(p1)+f(p2) --> Z^0(-->b(p3)+b~(p4)) + H(-->W^+(nu(p5),e^+(p6))W^-(e^-(p7),nub(p8)))'     
          call checkminzmass(1)
          nqcdjets=2
          plabel(3)='bq'
          plabel(4)='ba'
          q1=Q(5)*sqrt(xn)
          l1=l(5)*sqrt(xn)
          r1=r(5)*sqrt(xn)
        else
          call nprocinvalid()
        endif 

c--- print warning if we're below threshold
        if (hmass < 2._dp*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
      if (removebr) then
      write(6,*)
      write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
        stop
      endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif
                 
c-----------------------------------------------------------------------

        elseif (nproc == 109) then
c--  109 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4)) + H(-->Z(e^-(p5),e^+(p6))+Z(mu^-(p7),mu(p8)))'
        kcase=kZH__ZZ
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=0
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='ea'
        plabel(7)='el'
        plabel(8)='ea'
        plabel(9)='pp'
        
        ndim=16
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=zmass
        width3=zwidth
        q1=-1._dp
        l1=le
        r1=re
        l2=le
        r2=re

        mcfmplotinfo= (/ 34, 56, 78, 5678, (0,j=1,46) /)
        
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'    
            plabel(7)='ig'
            plabel(8)='ig'             
            plabel(9)='ig'             
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2._dp*brzee**3*zzbr
          endif
c--- print warning if we're below threshold
        if (hmass < 2._dp*zmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
      if (removebr) then
      write(6,*)
      write(6,*) 'Cannot remove H->ZZ BR, not defined below threshold'
        stop
      endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif
                 
c-----------------------------------------------------------------------

       elseif ((nproc == 111) .or. (nproc == 112)
     &    .or. (nproc == 119)) then
        kcase=kggfus0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(5)='pp'
        ndim=4
      
        n2=0
        n3=1
        mass3=hmass
        width3=hwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 111) then
c--  111 '  f(p1)+f(p2) --> H(-->b(p3)+bbar(p4))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'
          nqcdjets=2
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=0
            BrnRat=br
          endif
          
        elseif (nproc == 112) then
c--  112 '  f(p1)+f(p2) --> H(-->tau^-(p3)+tau^+(p4))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)       
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'
          nqcdjets=0
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            BrnRat=tautaubr
          endif

        elseif (nproc == 119) then
c--  119 '  f(p1)+f(p2) --> H(-->gamma^-(p3)+gamma^+(p4))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)       
          hdecaymode='gaga'
          plabel(3)='ga'
          plabel(4)='ga'
          hdecaymode='gaga'
          nqcdjets=0
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            BrnRat=gamgambr
          endif
        endif

      elseif ((nproc == 113) .or. (nproc == 123)
     &   .or. (nproc == 124) .or. (nproc == 125)
     &   .or. (nproc == 126) .or. (nproc == 127)) then
c--  113 '  f(p1)+f(p2) --> H (--> W^+(nu(p3)+e^+(p4)) + W^-(e^-(p5)+nu~(p6)))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)
        if     (nproc == 113) then
          kcase=kHWW_4l  
        elseif (nproc == 123) then
          kcase=kHWW_tb  
        elseif (nproc == 124) then
          kcase=kHWWint
        elseif (nproc == 125) then
          kcase=kHWWHpi
        elseif (nproc == 126) then
          kcase=kggWW4l
        elseif (nproc == 127) then
          kcase=kggWWbx
        endif           
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
c--- widths according to Kauer et al., for comparison with gg2WW
c        if (abs(hmass-140._dp) < 1.e-4_dp) hwidth=0.008235_dp
c        if (abs(hmass-170._dp) < 1.e-4_dp) hwidth=0.3837_dp
c        if (abs(hmass-200._dp) < 1.e-4_dp) hwidth=1.426_dp
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='pp'
        nqcdjets=0
        ndim=10
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

c--- if vector boson decays specified, initialize appropriately
        if ((nproc >= 123) .and. (vdecayid)) then
          call setvdecay(34,+1)
          call setvdecay(56,-1)
        endif

        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)
        
c--- print warning if we're below threshold
        if (hmass < 2._dp*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points'
        if (removebr) then
        write(6,*)
        write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
        stop
        endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen**2*wwbr
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'               
        endif
        
      elseif (nproc == 114) then
      kcase=kHWW2lq  
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
c--- widths according to Kauer et al., for comparison with gg2WW
c        if (abs(hmass-140._dp) < 1.e-4_dp) hwidth=0.008235_dp
c        if (abs(hmass-170._dp) < 1.e-4_dp) hwidth=0.3837_dp
c        if (abs(hmass-200._dp) < 1.e-4_dp) hwidth=1.426_dp
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        nqcdjets=2
        ndim=10
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)
        
c--- print warning if we're below threshold
        if (hmass < 2._dp*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points'
      if (removebr) then
      write(6,*)
      write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
        stop
      endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=2._dp*xn*brwen**2*wwbr
c        if (kpart==ktodk) BrnRat=BrnRat*(1._dp+as/pi)
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          nqcdjets=0
        endif
        
      elseif (nproc == 115) then
      kcase=kHWWdkW  
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        nqcdjets=2
        ndim=10
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=2._dp*xn*brwen**2*wwbr
c        if (kpart==ktodk) BrnRat=BrnRat*(1._dp+as/pi)
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          nqcdjets=0
        endif
        
c--- print warning if we're below threshold
        if (hmass < 2._dp*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points'
        if (removebr) then
        write(6,*)
        write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
        stop
        endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif
        
c-----------------------------------------------------------------------

      elseif (((nproc >= 116) .and. (nproc <= 118))
     &    .or. (nproc == 128)  .or. (nproc == 129)
     &    .or. (nproc == 130)  .or. (nproc == 131)
     &    .or. (nproc == 132)  .or. (nproc == 133) ) then
        if     (nproc <= 118) then
          kcase=kHZZ_4l  
        elseif (nproc == 128) then
          kcase=kHZZ_tb  
        elseif (nproc == 129) then
          kcase=kHZZint
        elseif (nproc == 130) then
          kcase=kHZZHpi
        elseif (nproc == 131) then
          kcase=kggZZ4l
        elseif (nproc == 132) then
          kcase=kggZZbx
        elseif (nproc == 133) then
          kcase=kHZZqgI
        endif
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(7)='pp'
        nqcdjets=0
        nwz=0
        ndim=10
        if(nproc==133) then 
           ndim=13
           nqcdjets=1 
           notag=1
        endif
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)
        
c--- print warning if we're below threshold
        if (hmass < 2._dp*zmass) then
          write(6,*)
          write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
          write(6,*) 'may not yield sensible results - check the number'
          write(6,*) 'of integration points'
          if (removebr) then
        write(6,*)
        write(6,*) 'Cannot remove H->ZZ BR, not defined below threshold'
        stop
          endif
          if (zerowidth) then
          write(6,*) 'zerowidth=.true. and higgs decay below threshold'
          stop
          endif
        endif
        
        if (   (nproc == 116)  .or. (nproc == 128)
     &    .or. (nproc == 129)  .or. (nproc == 130)
     &    .or. (nproc == 131)  .or. (nproc == 132)
     &    .or. (nproc == 133) ) then
c--  116 '  f(p1)+f(p2) --> H(-->Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6))'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)
c--- 128 '  f(p1)+f(p2) --> H(--> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6)) [top, bottom loops, exact]' 'L'
c--- 129 '  f(p1)+f(p2) --> H(--> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6)) [only H, gg->ZZ intf.]' 'L'
c--- 130 '  f(p1)+f(p2) --> H(--> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6)) [H squared and H, gg->ZZ intf.]' 'L'
c--- 131 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6) [gg only, (H + gg->ZZ) squared]' 'L'
c--- 132 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4)) + Z^0(mu^-(p5)+mu^+(p6) [(gg->ZZ) squared]' 'L'
c--- 133 '  f(p1)+f(p2) --> H(--> Z^0(mu^-(p3)+mu^+(p4)) + Z^0(e^-(p5)+e^+(p6) + f(p7)) [intf with SM, no cut on p7]' 'L'
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='ml'
          plabel(6)='ma'
          l1=le
          r1=re
          l2=le
          r2=re
          q1=-1._dp
          q2=-1._dp

c--- check runstring to change from (e,mu) Z decays to (e,e) or (mu,mu)
c          if     (index(runstring,'ELEL') > 0) then
c            plabel(5)='el'
c            plabel(6)='ea'
c            interference=.true.
c            vsymfact=0.25_dp
c          elseif (index(runstring,'MUMU') > 0) then
c            plabel(3)='ml'
c            plabel(4)='ma'
c            interference=.true.
c            vsymfact=0.25_dp
c          endif
           
c--- if vector boson decays specified, initialize appropriately
          if ((nproc >= 128) .and. (vdecayid)) then
            call setvdecay(34,0)
            call setvdecay(56,0)
            if (plabel(3) == plabel(5)) then
              if (plabel(3) .ne. 'nl') then
c------ for both Z decays to neutrinos, neglect interference effects for simplicity
                interference=.true.
                vsymfact=0.25_dp
              endif
            endif
          endif
          
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2._dp*brzee**2*zzbr  ! factor of 2 for identical particles
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'             
          endif

        elseif (nproc == 117) then
c--  117 '  f(p1)+f(p2) --> H(-->Z^0(3*(nu(p3)+nu~(p4)))+ Z^0(mu^-(p5)+mu^+(p6))'
          plabel(3)='nl'
          plabel(4)='na'
          plabel(5)='ml'
          plabel(6)='ma'
          l1=le
          r1=re
          l2=ln*sqrt(3._dp)
          r2=rn*sqrt(3._dp)      
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2._dp*brzee*brznn*zzbr  ! factor of 2 for identical particles
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'             
          endif
        elseif (nproc == 118) then
c--  118 '  f(p1)+f(p2) --> H(-->Z^0(mu^-(p3)+mu^+(p4)) + Z^0(b(p5)+b~(p6))'
          nqcdjets=2
          plabel(3)='ml'
          plabel(4)='ma'
          plabel(5)='bq'
          plabel(6)='ba'
          l1=le
          r1=re
          l2=l(5)*sqrt(xn)
          r2=r(5)*sqrt(xn)
        else
          call nprocinvalid()
        endif

        elseif (nproc == 120) then
        kcase=kHi_Zga
        call checkminzmass(1)
        nqcdjets=0
        plabel(6)='pp'
        ndim=7
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
      
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 345, (0,j=1,48) /)
        
c--  120 '  f(p1)+f(p2) --> H(-->Z^0(mu^-(p3)+mu^+(p4)) + gamma(p5)')'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='ga'
        l1=le
        r1=re
        q1=-1._dp
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee*zgambr 
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
        endif

        elseif (nproc == 121) then
        kcase=kHi_Zga
        nqcdjets=0
        plabel(6)='pp'
        ndim=7
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
      
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 345, (0,j=1,48) /)
        
c--  121 '  f(p1)+f(p2) --> H(-->Z^0(3*(nu(p3)+nu~(p4))) + gamma(p5)')'
c--      '  f(p1)+f(p2) --> H (for total Xsect)' (removebr=.true.)
        plabel(3)='nl'
        plabel(4)='na'
        q1=zip
        l1=ln*sqrt(3._dp)
        r1=rn*sqrt(3._dp)      
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brznn*zgambr
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
        endif

c-----------------------------------------------------------------------

      elseif ( (nproc == 1281)  .or. (nproc == 1291)
     &    .or. (nproc == 1301)  .or. (nproc == 1311)
     &    .or. (nproc == 1321)
     &    .or. (nproc == 1282)  .or. (nproc == 1292)
     &    .or. (nproc == 1302)  .or. (nproc == 1312)
     &    .or. (nproc == 1322) ) then
        if     ((nproc == 1281) .or. (nproc == 1282)) then
          kcase=kHVV_tb  
        elseif ((nproc == 1291) .or. (nproc == 1292)) then
          kcase=kHVVint
        elseif ((nproc == 1301) .or. (nproc == 1302)) then
          kcase=kHVVHpi
        elseif ((nproc == 1311) .or. (nproc == 1312)) then
          kcase=kggVV4l
        elseif ((nproc == 1321) .or. (nproc == 1322)) then
          kcase=kggVVbx
        endif
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='nl'
        plabel(6)='na'
        plabel(7)='pp'
        l1=le
        r1=re
        q1=-1._dp
        l2=ln
        r2=rn
        q2=zip
        nqcdjets=0
        nwz=0
        ndim=10
        n2=1
        n3=1

c-- parameters for phase space
        doipsgen=.true.
        maxipsgen=2

        if   ( (nproc == 1282)  .or. (nproc == 1292)
     &    .or. (nproc == 1302)  .or. (nproc == 1312)
     &    .or. (nproc == 1322) ) then
          nuflav=3
        else
          nuflav=1
        endif

        mcfmplotinfo= (/ 34, 56, 45, 36, 3456, (0,j=1,45) /)
 
c-----------------------------------------------------------------------

      elseif ((nproc >= 136) .and. (nproc <= 138)) then
        kcase=kH_1jet
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        call setmb_msbar
        mb=zip
        ndim=7
        plabel(3)='qb'
        plabel(4)='ab'
        plabel(5)='bq'
        hdecaymode='bqba'
        nqcdjets=1

        n2=0
        n3=1
        mass3=hmass
        width3=hwidth
        
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if ( (nproc == 137) .and.  
     &       ((kpart==kvirt) .or. (kpart==ktota)) ) then
          write(6,*) 'This process number is not suitable for the'
          write(6,*) 'NLO calculation. Please run processes'
          write(6,*) '136 (virtual+real) and 137 (real) separately.'
          stop
        endif
        if ( (nproc == 138) .and. (kpart.ne.kreal) ) then
          write(6,*) 'This process number is not suitable for such a'
          write(6,*) 'calculation. Please run process 138 (real) only.'
          stop
        endif
             
        if     (nproc == 136) then
c--  136 '  f(p1)+f(p2) --> H (no BR) + b(p5) [+g(p6)]'
          isub=1
          plabel(6)='pp'
        elseif (nproc == 137) then
c--  137 '(p1)+f(p2) --> H (no BR) + b~(p5) [+b(p6)]'
          isub=2
          plabel(5)='ba'
          plabel(6)='bq'
        elseif (nproc == 138) then
c--  138 '  f(p1)+f(p2) --> H (no BR) + b(p5) + b~(p6) [both observed]'
          isub=2
          plabel(5)='ba'
          plabel(6)='bq'
          nqcdjets=2
        endif
        
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          BrnRat=br
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc == 141) 
     &   .or. (nproc == 142) 
     &   .or. (nproc == 144) 
     &   .or. (nproc == 145) 
     &   .or. (nproc == 146)
     &   .or. (nproc == 147)
     &   .or. (nproc == 148)
     &   .or. (nproc == 149)
     &   .or. (nproc == 150)
     &   .or. (nproc == 151)) then
        ndim=16
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth
        bbproc=.true.
        
        mcfmplotinfo= (/ 34, 78, 345, 678, (0,j=1,46) /)
        
        if (nproc == 141) then
c--  141 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+e^-(p7)+nu~(p8)'
c--      '  f(p1)+f(p2) --> t t~ (with BR for total Xsect)' (removebr=.true.)
          kcase=ktt_bbl
          nqcdjets=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=(brwen*brtop)**2
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            plabel(7)='ig'
            plabel(8)='ig'
            nqcdjets=0
            bbproc=.false.
          endif
        elseif (nproc == 142) then
c--  142 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+e^-(p7)+nu~(p8) [radiation in top decay]'
          kcase=ktt_ldk
          nqcdjets=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=(brwen*brtop)**2
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            plabel(7)='ig'
            plabel(8)='ig'
            nqcdjets=0
            bbproc=.false.
          endif
        elseif (nproc == 144) then
c--  144 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->b~(p6)+e^-(p7)+nu~(p8))'
c--           (uncorrelated)'
          kcase=ktt_bbu
          nqcdjets=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=(brwen*brtop)**2
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            plabel(7)='ig'
            plabel(8)='ig'
            nqcdjets=0
            bbproc=.false.
          endif
        elseif (nproc == 145) then
c--  145 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->b~(p6)+e^-(p7)+nu~(p8))'
c--           (uncorrelated) [rad.in.top.dk]'
          kcase=ktt_udk
          nqcdjets=2
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=(brwen*brtop)**2
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            plabel(7)='ig'
            plabel(8)='ig'
            nqcdjets=0
            bbproc=.false.
          endif
        elseif (nproc == 146) then
c--  146 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+q(p7)+q~(p8)'
          kcase=ktt_bbh
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='pp'
          plabel(8)='pp'
          plabel(9)='pp'
          nqcdjets=4
        elseif (nproc == 147) then
c--  147 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+q(p7)+q~(p8) [radiation in top decay]'
          kcase=ktt_hdk
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='pp'
          plabel(8)='pp'
          plabel(9)='pp'
          nqcdjets=4
        elseif (nproc == 148) then
c--  148 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6))+q(p7)+q~(p8) [radiation in hadronic W decay]'
          kcase=ktthWdk
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='pp'
          plabel(8)='pp'
          plabel(9)='pp'
          nqcdjets=4
        elseif (nproc == 149) then
c---  149 '  f(p1)+f(p2) --> t(-->q(p3)+q~(p4)+b(p5))+t~(-->b~(p6)+e-(p7)+nu~(p8))'
          kcase=ktt_bbh
          plabel(3)='pp'
          plabel(4)='pp'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          nqcdjets=4
        elseif (nproc == 150) then
c---  150 '  f(p1)+f(p2) --> t(-->q(p3)+q~(p4)+b(p5))+t~(-->b~(p6)+e-(p7)+nu~(p8)) [radiation in top decay]'
          kcase=ktt_hdk
          plabel(3)='pp'
          plabel(4)='pp'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          nqcdjets=4
        elseif (nproc == 151) then
c---  151 '  f(p1)+f(p2) --> t(-->q(p3)+q~(p4)+b(p5))+t~(-->b~(p6)+e-(p7)+nu~(p8)) [radiation in hadronic W decay]'
          kcase=ktthWdk
          plabel(3)='pp'
          plabel(4)='pp'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          nqcdjets=4
        endif 
c-----------------------------------------------------------------------

      elseif (nproc == 143) then
c--  143 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->nu~(p7)+e^-(p8)+b~(p6))+g(p9)'
c--      '  f(p1)+f(p2)-->t(p345)+t~(p678)+g(p9)' (removebr=.true.)
        kcase=kqq_ttg
        nwz=1
        ndim=19
        nqcdjets=3
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='el'
        plabel(8)='na'
        plabel(9)='pp'

        mcfmplotinfo= (/ 34, 78, 345, 678, (0,j=1,46) /)
        
c--- total cross-section             
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=(brwen*brtop)**2
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'               
          plabel(7)='ig'
          plabel(8)='ig'
          nqcdjets=1
        endif

c-----------------------------------------------------------------------

      elseif (nproc == 157) then
c--  157 '  f(p1)+f(p2) --> t t~ (for total Xsect)'
        kcase=ktt_tot
        nqcdjets=0
        ndim=4
        mass2=mt
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
      elseif (nproc == 158) then
c--  158 '  f(p1)+f(p2) --> b b~ (for total Xsect)'
        kcase=kbb_tot
c      nflav=4
        nqcdjets=0
        ndim=4
        mass2=mb
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
      elseif (nproc == 159) then
c--  159 '  f(p1)+f(p2) --> c c~ (for total Xsect)'
        kcase=kcc_tot
c      nflav=3
        nqcdjets=0
        ndim=4
        mass2=mc
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
 
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
      elseif (nproc == 160) then
      if  ((kpart==ktota)
     & .or.(kpart==kvirt)
     & .or.(kpart==kreal)) then
          write(6,*) 'This process number is available only at LO'
          write(6,*) 'Please set part = lord and rerun'
             stop
      endif
c--  160 '  f(p1)+f(p2) --> t t~ +jet (for total Xsect)'
        kcase=ktt_glu
        nqcdjets=1
        ndim=7
        mass2=mt
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'

c-----------------------------------------------------------------------

      elseif ((nproc == 161) .or. (nproc == 163)) then
c--  161 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+q(p6) [t-channel]'
c--      '  f(p1)+f(p2) --> t(no BR) + q(p6)' (removebr=.true.)
        kcase=kbq_tpq
        isub=1
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='qj'
        plabel(7)='pp'
        nwz=+1
      
      if (nproc == 161) then ! usual approach, mb=0 
c--- extra b that can appear at NLO is massless
          masslessb=.true.
      else                     ! proper ACOT, mb>0 (must run 231 LO)
        masslessb=.false.
      endif

c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=zip
        write(6,*) 'Enforcing mb=0 for this process!'
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 345, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1   
        endif

      elseif (nproc == 162) then
c--  162 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+q(p6) [decay]'
        kcase=kttdkay
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='qj'
        plabel(7)='pp'
        nwz=+1

        if (kpart==klord) then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '161 (lord) or process 162 (virt+real).'
          stop
        endif
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=zip
        write(6,*) 'Enforcing mb=0 for this process!'
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
        
        mcfmplotinfo= (/ 34, 345, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
        endif

      elseif ((nproc == 166) .or. (nproc == 168)) then
c--  166 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+q(p6) [t-channel]''
c--      '  f(p1)+f(p2) --> t~(no BR) + q(p6)' (removebr=.true.)
        kcase=kbq_tpq
        isub=1
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ba'
        plabel(6)='qj'
        plabel(7)='pp'
        nwz=-1

      if (nproc == 166) then ! usual approach, mb=0 
c--- extra b that can appear at NLO is massless
          masslessb=.true.
      else                     ! proper ACOT, mb>0 (must run 236 LO)
        masslessb=.false.
      endif

c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=zip
        write(6,*) 'Enforcing mb=0 for this process!'
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 345, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
        endif
        
      elseif (nproc == 167) then
c--  167 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+q(p6) [decay]'
        kcase=kttdkay
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ba'
        plabel(6)='qj'
        plabel(7)='pp'
        nwz=-1

        if (kpart==klord) then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '166 (lord) or process 167 (virt+real).'
          stop
        endif
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        mb=zip
        write(6,*) 'Enforcing mb=0 for this process!'
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
             
        mcfmplotinfo= (/ 34, 345, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
        endif

c-----------------------------------------------------------------------

      elseif (nproc == 171) then
c--  171 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6)) [s-channel]'
c--      '  f(p1)+f(p2) --> t(no BR) + b~(p6)' (removebr=.true.)
        kcase=kt_bbar
        isub=2
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        nwz=1
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
        
        mcfmplotinfo= (/ 34, 345, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
        endif
        
      elseif (nproc == 172) then
c--  172 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6)) [decay]'
        kcase=ktdecay
        nqcdjets=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        nwz=1
        
        if (kpart==klord) then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '171 (lord) or process 172 (virt+real).'
          stop
        endif
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 345, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1           
        endif
        
      elseif (nproc == 176) then
c--  176 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+b(p6)) [s-channel]'
c--      '  f(p1)+f(p2) --> t~(no BR) + b(p6)' (removebr=.true.)
        kcase=kt_bbar
        isub=2
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ba'
        plabel(6)='bq'
        plabel(7)='pp'
        nwz=-1
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
        
        mcfmplotinfo= (/ 34, 345, (0,j=1,48) /)
        
         if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1
        endif
             
      elseif (nproc == 177) then
c--  177 '  f(p1)+f(p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+b(p6)) [decay]'
        kcase=ktdecay
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ba'
        plabel(6)='bq'
        plabel(7)='pp'
        nwz=-1
        
        if (kpart==klord) then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '176 (lord) or process 177 (virt+real).'
          stop
        endif
        
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=9
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 345, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=1
        endif
        
c-----------------------------------------------------------------------

      elseif (nproc == 180) then
c--  180 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(p5)'
        kcase=kW_tndk
        nqcdjets=0
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='ig'
        plabel(6)='pp'
        mass2=mt
        nflav=5
        nwz=-1

        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
             
      elseif (nproc == 181) then
c--  181 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(nu(p5)+e^+(p6)+b(p7))'
        kcase=kW_twdk
        nqcdjets=1
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='nl'
        plabel(6)='ea'
        plabel(7)='bq'
        plabel(8)='pp'
        nflav=5
        nwz=-1

        ndim=13
        mb=zip
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
             
        mcfmplotinfo= (/ 34, 567, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          nqcdjets=0
        endif
             
      elseif (nproc == 182) then
c--  182 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(nu(p5)+e^+(p6)+b(p7)) [decay]'
        
        kcase=kWtdkay
        nqcdjets=1
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='nl'
        plabel(6)='ea'
        plabel(7)='bq'
        plabel(8)='pp'
        nflav=5
        nwz=-1

        if (kpart==klord) then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '181 (lord) or process 182 (virt+real).'
          stop
        endif
        
        ndim=13
        mb=zip
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
             
        mcfmplotinfo= (/ 34, 567, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          nqcdjets=0
        endif
             
      elseif (nproc == 183) then
c--  183 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(nu(p5)+e^+(p6)+b(p7))+b(p8)'
        kcase=kWtbwdk
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='nl'
        plabel(6)='ea'
        plabel(7)='bq'
        plabel(8)='pp'
        nflav=5
        nwz=-1

        ndim=16
c--- (this process can also be used for non-zero mb)
        mb=zip
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
             
        mcfmplotinfo= (/ 34, 567, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          nqcdjets=1
        endif
             
      elseif (nproc == 184) then
c--  184 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(p5)+b(p6) [massive b]'
        kcase=kWtbndk
        nqcdjets=2
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='bq'
        plabel(6)='ba'
        mass2=mt
        nflav=5
        nwz=-1

        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
             
      elseif (nproc == 185) then
c--  185 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+tbar(p5)'
        kcase=kW_tndk
        nqcdjets=0
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='ig'
        plabel(6)='pp'
        mass2=mt
        nflav=5
        nwz=+1
        
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
             
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
      elseif (nproc == 186) then
c--  186 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+t~(e^-(p5)+nu~(p6)+bbar(p7))'
        kcase=kW_twdk
        nqcdjets=1
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='ba'
        plabel(8)='pp'
        nflav=5
        nwz=+1
        
        ndim=13
        mb=0
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 567, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          nqcdjets=0
        endif
             
      elseif (nproc == 187) then
c--  182 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+t~(e^-(p5)+nu~(p6)+bbar(p7)) [decay]'
        
        kcase=kWtdkay
        nqcdjets=1
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='bq'
        plabel(8)='pp'
        nflav=5
        nwz=+1

        if (kpart==klord) then
          write(6,*) 'This process number can not be used for a'
          write(6,*) 'LO calculation. Please run either process'
          write(6,*) '186 (lord) or process 187 (virt+real).'
          stop
        endif
        
        ndim=13
        mb=0
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth
             
        mcfmplotinfo= (/ 34, 567, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          nqcdjets=0
        endif
             
      elseif ((nproc >= 200) .and. (nproc <= 210)) then
        kcase=khttjet
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=1
        plabel(5)='pp'
 
        ndim=7
        n2=0
        n3=1
        mass3=hmass
        width3=hwidth
        
        if     (nproc == 201) then
c--  201 '  f(p1)+f(p2)--> H(-->b(p3)+b~(p4)) + f(p5) [full mt dep.]'
c--      '  f(p1)+f(p2)--> H(p3+p4) + f(p5) (for total Xsect)' (removebr=.true.)
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'  
          nqcdjets=3

          mcfmplotinfo= (/ 34, (0,j=1,49) /)

          if (removebr) then
            BrnRat=br
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=1
          endif

        elseif (nproc == 202) then
c--  202 '  f(p1)+f(p2)--> H (-> tau(p3) tau~(p4)) + f(p5) [full mt dep.]'
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'

          mcfmplotinfo= (/ 34, (0,j=1,49) /)

          if (removebr) then
            BrnRat=tautaubr
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=1
          endif

        elseif ((nproc == 203) .or. (nproc == 204)
     &     .or. (nproc == 210)) then
          kcase=kggfus1
          nqcdjets=1
          plabel(5)='pp'
          plabel(6)='pp'
          ndim=7
      
          n2=0
          n3=1

          mcfmplotinfo= (/ 34, (0,j=1,49) /)

          if     (nproc == 203) then
c--  203 '  f(p1)+f(p2) -->H(-->b(p3)+b~(p4)) + f(p5)'
c--      '  f(p1)+f(p2)--> H(p3+p4) + f(p5) (for total Xsect)' (removebr=.true.)
            hdecaymode='bqba'
            plabel(3)='bq'
            plabel(4)='ba'
            nqcdjets=3
            if (removebr) then
              BrnRat=br
              plabel(3)='ig'
              plabel(4)='ig'
              nqcdjets=1
            endif
          elseif (nproc == 204) then
c--  204 '  f(p1)+f(p2) -->H(-->tau^-(p3)+tau^+(p4)) + f(p5)'
            hdecaymode='tlta'
            plabel(3)='tl'
            plabel(4)='ta'
            nqcdjets=1
            if (removebr) then
              Brnrat=tautaubr
              plabel(3)='ig'
              plabel(4)='ig'
            endif
          elseif (nproc == 210) then
c--  210 '  f(p1)+f(p2) -->H(-->gamma(p3)+gamma(p4)) + f(p5)'
            hdecaymode='gaga'
            plabel(3)='ga'
            plabel(4)='ga'
            nqcdjets=1
            if (removebr) then
              BrnRat=gamgambr
              plabel(3)='ig'
              plabel(4)='ig'
            endif
          endif

        
        elseif (nproc == 206) then
c--  206 '  f(p1)+f(p2)--> A(-->b(p3)+b~(p4)) + f(p5) [full mt dep.]'
c--      '  f(p1)+f(p2)--> A(p3+p4) + f(p5) (for total Xsect)' (removebr=.true.)
          kcase=kattjet
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'  
          nqcdjets=3

          mcfmplotinfo= (/ 34, (0,j=1,49) /)

          if (removebr) then
            BrnRat=br
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=1
          endif

        elseif (nproc == 207) then
c--  207 '  f(p1)+f(p2)--> A (--> tau(p3) tau~(p4)) + f(p5) [full mt dep.]'
         kcase=kattjet
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'

          mcfmplotinfo= (/ 34, (0,j=1,49) /)

          if (removebr) then
            BrnRat=tautaubr
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=1
          endif
        endif

        if     (nproc == 208) then
c-- 208 '  f(p1)+f(p2) --> H(-->W^+(p3,p4)W^-(p5,p6)) + f(p7)'
          kcase=kHWWjet
          ndim=13
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='el'
          plabel(6)='na'
          plabel(7)='pp'
          plabel(8)='pp'
          nqcdjets=1
          n2=1
          n3=1
          mass2=wmass
          width2=wwidth
          mass3=wmass
          width3=wwidth

          mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)        

c--- print warning if we're below threshold
          if (hmass < 2._dp*wmass) then
          write(6,*)
          write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
          write(6,*) 'may not yield sensible results - check the number'
          write(6,*) 'of integration points and the value of zerowidth'
        if (removebr) then
        write(6,*)
      write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
          stop
        endif
          if (zerowidth) then
          write(6,*) 'zerowidth=.true. and higgs decay below threshold'
          stop
          endif
          endif
        
          if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=wwbr*brwen**2
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          endif
        endif

       if ((nproc == 209)) then
          kcase=kHZZjet
          l1=le
          l2=le
          r1=re
          r2=re
          ndim=13
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='ml'
          plabel(6)='ma'
          plabel(7)='pp'
          plabel(8)='pp'
          nqcdjets=1
          n2=1
          n3=1
          mass2=zmass
          width2=zwidth
          mass3=zmass
          width3=zwidth
          call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

          mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)

c--- print warning if we're below threshold
        if (hmass < 2._dp*zmass) then
          write(6,*)
          write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
          write(6,*) 'may not yield sensible results - check the number'
          write(6,*) 'of integration points'
          if (zerowidth) then
          write(6,*) 'zerowidth=.true. and higgs decay below threshold'
          stop
          endif
        endif
        
          l1=le
          r1=re
          l2=le
          r2=re
          if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2._dp*brzee**2*zzbr  ! factor of 2 for identical particles
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
          endif

        endif

c-----------------------------------------------------------------------

      elseif ((nproc == 211) .or. (nproc == 212)) then
        kcase=kqq_Hqq
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nwz=2
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=10
        n2=0
        n3=1
      
        mass3=hmass
        width3=hwidth

        mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)

        if     (nproc == 211) then
c--  211 '  f(p1)+f(p2)--> H(-->b(p3)+b~(p4))+f(p5)+f(p6) [WBF]'
c--      '  f(p1)+f(p2)--> H(p3+p4)+f(p5)+f(p6) [WBF]' (removebr=.true.)
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'
          nqcdjets=4
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=2
            BrnRat=br
          endif
        elseif (nproc == 212) then
c--  212 '  f(p1)+f(p2)--> H(-->tau-(p3)+tau+(p4))+f(p5)+f(p6) [WBF]'
c--      '  f(p1)+f(p2)--> H(p3+p4)+f(p5)+f(p6) [WBF]' (removebr=.true.)
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'
          nqcdjets=2
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            Brnrat=tautaubr
          endif
        endif
          
      elseif (nproc == 213) then
        kcase=kqq_HWW
        mb=zip
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nwz=2
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='pp'
        plabel(8)='pp'
        plabel(9)='pp'
        ndim=16
        nqcdjets=2
c      notag=1 ! If only one jet is required
c      notag=0 ! FOR CHECKING VS 211
      
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)        

c--- print warning if we're below threshold
        if (hmass < 2._dp*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
      if (removebr) then
      write(6,*)
      write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
        stop
      endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif
        
        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=wwbr*brwen**2
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif
          
      elseif (nproc == 214) then
        kcase=kqq_HZZ
        l1=le
        r1=re
        l2=le
        r2=re
        mb=zip
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nwz=2
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='ml'
        plabel(6)='ma'
        plabel(7)='pp'
        plabel(8)='pp'
        plabel(9)='pp'
        ndim=16
        nqcdjets=2
      
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 56, 78, 3456, (0,j=1,46) /)

c--- print warning if we're below threshold
        if (hmass < 2._dp*zmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
      if (removebr) then
      write(6,*)
      write(6,*) 'Cannot remove H->ZZ BR, not defined below threshold'
        stop
      endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif
        
        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=2._dp*zzbr*brzee**2  ! factor of 2 for identical particles
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif
          
      elseif (nproc == 215) then
        kcase=kqq_Hgg
      hdecaymode='gaga'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(3)='ga'
        plabel(4)='ga'
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=10
        nqcdjets=2
      
        n2=0
        n3=1
        mass3=hmass
        width3=hwidth
       
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
      if (removebr) then
      BrnRat=gamgambr
      plabel(3)='ig'
      plabel(4)='ig'
      endif
          
      elseif ((nproc == 216) .or. (nproc == 217)) then
        kcase=kqqHqqg
        mb=zip
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nwz=2
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=13
        n2=0
        n3=1

        mass3=hmass
        width3=hwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
       if     (nproc == 216) then
c-- 216 '  f(p1)+f(p2)--> H(-->b(p3)+b~(p4))+f(p5)+f(p6)+f(p7) [WBF+jet]'
c--     '  f(p1)+f(p2)--> H(p3+p4)+f(p5)+f(p6)+f(p7) [WBF+jet]' (removebr=.true.)
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'
          nqcdjets=5
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=3
            BrnRat=br
          endif
        elseif (nproc == 217) then
c-- 217 '  f(p1)+f(p2)--> H(-->tau-(p3)+tau+(p4))+f(p5)+f(p6)+f(p7) [WBF+jet]'
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'
          nqcdjets=3
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=1
            Brnrat=tautaubr
          endif
        endif

c-----------------------------------------------------------------------
        
      elseif ((nproc == 220) .or. (nproc == 2201)) then
        kcase=kqqZZqq
        if (nproc == 2201) VVstrong=.true.
        call checkminzmass(1)
        call checkminzmass(2)
        mb=zip
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nwz=2
C----charge for 34-line
        q1=-1._dp
        l1=le
        r1=re
C----charge for 56-line
        q2=-1._dp
        l2=le
        r2=re
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='ml'
        plabel(6)='ma'
        plabel(7)='pp'
        plabel(8)='pp'
        ndim=16
        nqcdjets=2
      
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 56, 78, 3456, (0,j=1,46) /)

        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=2._dp*brzee**2  ! factor of 2 for identical particles
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif
          

c-----------------------------------------------------------------------

      elseif (nproc == 221) then
        kcase=ktautau
c--  221 '  f(p1)+f(p2)--> tau^-(-->e^-(p3)+nu~_e(p4)+nu_tau(p5))+tau^+(-->nu~_tau(p6)+nu_e(p7)+e^+(p8))'
c--      '  f(p1)+f(p2)--> tau tau~ [for total Xsect]' (removebr=.true.)
        plabel(3)='el'
        plabel(4)='na'
        plabel(5)='nl'
        plabel(6)='na'
        plabel(7)='nl'
        plabel(8)='ea'
        nqcdjets=0
        nwz=1
        ndim=16
        n2=1
        n3=1
        mass2=mtau
        width2=tauwidth
        mass3=mtau
        width3=tauwidth

        mcfmplotinfo= (/ 34, 78, 345, 678, (0,j=1,46) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brtau**2
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          plabel(8)='ig'
        endif

c-----------------------------------------------------------------------

      elseif ((nproc == 222) .or. (nproc == 2221)) then
        kcase=kqqZZqq
        if (nproc == 2221) VVstrong=.true.
        call checkminzmass(1)
        mb=zip
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nwz=2
C----charge for 34-line
        q1=-1._dp
        l1=le
        r1=re
C----charge for 56-line
        q2=zip
        l2=ln
        r2=rn
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='nl'
        plabel(6)='na'
        plabel(7)='pp'
        plabel(8)='pp'
        ndim=16
        nqcdjets=2
      
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 56, 78, 3456, (0,j=1,46) /)

        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=brzee*brznn
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif
          
c-----------------------------------------------------------------------

      elseif ((nproc == 224) .or. (nproc == 2241)) then
        kcase=kqqWWqq
        if (nproc == 2241) VVstrong=.true.
        mb=zip
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='nl'
        plabel(6)='na'
        plabel(7)='pp'
        plabel(8)='pp'
        ndim=16
        nqcdjets=2
      
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 36, 45, 78, 3456, (0,j=1,46) /)

        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=brwen**2
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif
          
c-----------------------------------------------------------------------

      elseif ((nproc == 226) .or. (nproc == 2261)) then
        kcase=kqqVVqq
        if (nproc == 2261) VVstrong=.true.
        mb=zip
C----charge for 34-line
        q1=-1._dp
        l1=le
        r1=re
C----charge for 56-line
        q2=zip
        l2=ln
        r2=rn
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='nl'
        plabel(6)='na'
        plabel(7)='pp'
        plabel(8)='pp'
        ndim=16
        nqcdjets=2
      
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

c-- parameters for phase space
        doipsgen=.true.
        maxipsgen=2

        mcfmplotinfo= (/ 34, 56, 36, 45, 78, 3456, (0,j=1,44) /)
          
c-----------------------------------------------------------------------

      elseif ((nproc == 228) .or. (nproc == 229)
     &   .or. (nproc == 2281) .or. (nproc == 2291)) then
        if ((nproc == 2281) .or. (nproc == 2291)) VVstrong=.true.
        kcase=kqqWWss
        mb=zip
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(7)='pp'
        plabel(8)='pp'
        ndim=16
        nqcdjets=2
      
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth
c--- nwz will be used to signal charge of W decays
        if ((nproc == 228) .or. (nproc == 2281)) then
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='nl'
          plabel(6)='ma'
        endif
        if ((nproc == 229) .or. (nproc == 2291)) then
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
          plabel(5)='ml'
          plabel(6)='na'
        endif
        mcfmplotinfo= (/ 34, 56, 78, 3456, (0,j=1,46) /)

        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=brwen**2
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif
          
c-----------------------------------------------------------------------

      elseif ((nproc == 223) .or. (nproc == 225)
     &   .or. (nproc == 2231) .or. (nproc == 2251)) then
        if ((nproc == 2231) .or. (nproc == 2251)) VVstrong=.true.
        kcase=kqqWZqq
        mb=zip
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(7)='pp'
        plabel(8)='pp'
        ndim=16
        nqcdjets=2
      
C----fix these charge assignments???
        q1=-1._dp
        l1=le
        r1=re
        q2=-1._dp
        l2=le
        r2=re

        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=zmass
        width3=zwidth
c--- nwz will be used to signal charge of W decays
        if ((nproc == 223) .or. (nproc == 2231)) then
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='ml'
          plabel(6)='ma'
        endif
        if ((nproc == 225) .or. (nproc == 2251)) then
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
          plabel(5)='ml'
          plabel(6)='ma'
        endif
        mcfmplotinfo= (/ 34, 56, 78, 3456, (0,j=1,46) /)

        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=brwen*brzee
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif
          
c-----------------------------------------------------------------------

      elseif ((nproc == 231) .or. (nproc == 236)) then
c--  231 '  f(p1)+f(p2) --> t(p3)+b~(p4)+q(p5) [t-channel]'
c--  236 '  f(p1)+f(p2) --> t~(p3)+b(p4)+q(p5) [t-channel]'
        kcase=kqg_tbq
        nqcdjets=1
        ndim=7
        mass2=mt
        mass3=mb
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
        plabel(6)='pp'
      if (nproc == 231) then
        nwz=+1
      else
        nwz=-1
      endif

c---  in the SM, the logical:: fourthgen should be false
c---  for BSM calculations, it should be true and it indicates that
c---   5 flavours should be used in the PDF and alpha-s
        fourthgen=.false.

      if (fourthgen) then
c--- BSM: full 5 light flavours
          nflav=5
          bmass=4.7_dp  !  set b-mass to its usual value
      else
c--- SM: only 4 light flavours
          nflav=4
          bmass=1001._dp !  enforce 4-flavour running in alfamz.f
        endif
      
c--- set up correct scales and as on heavy and light quark lines
        facscale_H=initfacscale_H
        facscale_L=initfacscale_L
        renscale_H=initrenscale_H
        renscale_L=initrenscale_L
c--- make sure it works even if not specifying separate scales
        if (initrenscale_L == zip) then 
        facscale_H=facscale
        facscale_L=facscale
        renscale_H=scale
        renscale_L=scale
      endif
      
        b0=(xn*11._dp-2._dp*nflav)/6._dp
      as_H=alphas(abs(renscale_H),amz,nlooprun)
      as_L=alphas(abs(renscale_L),amz,nlooprun)

      
c-----------------------------------------------------------------------

      elseif ((nproc == 232) .or. (nproc == 237)) then
c--  232 '  f(p1)+f(p2) --> t(p3)+b~(p4)+q(p5)+q(p6) [t-channel]'
c--  237 '  f(p1)+f(p2) --> t~(p3)+b(p4)+q(p5)+q(p6) [t-channel]'
        kcase=kqgtbqq
        nqcdjets=2
        ndim=10
        mass2=mt
        mass3=mb
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
        plabel(6)='pp'
      if (nproc == 232) then
        nwz=+1
      else
        nwz=-1
      endif

c---  in the SM, the logical:: fourthgen should be false
c---  for BSM calculations, it should be true and it indicates that
c---   5 flavours should be used in the PDF and alpha-s
        fourthgen=.false.

      if (fourthgen) then
c--- BSM: full 5 light flavours
          nflav=5
        bmass=4.7_dp  !  set b-mass to its usual value
      else
c--- SM: only 4 light flavours
          nflav=4
          bmass=1001._dp !  enforce 4-flavour running in alfamz.f
        endif
      
c--- set up correct scales and as on heavy and light quark lines
        facscale_H=initfacscale_H
        facscale_L=initfacscale_L
        renscale_H=initrenscale_H
        renscale_L=initrenscale_L
c--- make sure it works even if not specifying separate scales
        if (initrenscale_L == zip) then 
        facscale_H=facscale
        facscale_L=facscale
        renscale_H=scale
        renscale_L=scale
      endif
      
        b0=(xn*11._dp-2._dp*nflav)/6._dp
      as_H=alphas(abs(renscale_H),amz,nlooprun)
      as_L=alphas(abs(renscale_L),amz,nlooprun)

      
c-----------------------------------------------------------------------

      elseif ((nproc == 233) .or. (nproc == 238)
     &   .or. (nproc == 234) .or. (nproc == 239)) then
c--  233 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6)+q(p7) [t-channel]'
c--  234 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6)+q(p7) [t-channel, rad. in decay]'
c--  238 '  f(p1)+f(p2) --> t~(-->e-(p3)+nu~(p4)+b~(p5))+b(p6)+q(p7) [t-channel]'
c--  239 '  f(p1)+f(p2) --> t~(-->e-(p3)+nu~(p4)+b~(p5))+b(p6)+q(p7) [t-channel, rad. in decay]'
c--      '  f(p1)+f(p2) --> t(no BR) + b~(p6) + q(p7)' (removebr=.true.)
      if ((nproc == 233) .or. (nproc == 238)) then
        kcase=k4ftwdk
      else
        kcase=kdk_4ft
      endif
      nqcdjets=3
      if ((nproc == 233) .or. (nproc == 234)) then
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
      else
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
          plabel(5)='ba'
          plabel(6)='bq'
      endif
      plabel(7)='pp'
      plabel(8)='pp'
            
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=12
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 345, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=2   
        endif

c---  in the SM, the logical:: fourthgen should be false
c---  for BSM calculations, it should be true and it indicates that
c---   5 flavours should be used in the PDF and alpha-s
        fourthgen=.false.

      if (fourthgen) then
c--- BSM: full 5 light flavours
          nflav=5
          bmass=4.7_dp  !  set b-mass to its usual value
      else
c--- SM: only 4 light flavours
          nflav=4
          bmass=1001._dp !  enforce 4-flavour running in alfamz.f
        endif
      
c--- set up correct scales and as on heavy and light quark lines
        facscale_H=initfacscale_H
        facscale_L=initfacscale_L
        renscale_H=initrenscale_H
        renscale_L=initrenscale_L
c--- make sure it works even if not specifying separate scales
        if (initrenscale_L == zip) then 
        facscale_H=facscale
        facscale_L=facscale
        renscale_H=scale
        renscale_L=scale
        endif
      
        b0=(xn*11._dp-2._dp*nflav)/6._dp
        as_H=alphas(abs(renscale_H),amz,nlooprun)
        as_L=alphas(abs(renscale_L),amz,nlooprun)
      
c-----------------------------------------------------------------------

      elseif ((nproc == 235) .or. (nproc == 240)) then
c--  235 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+b~(p6)+q(p7)+f(p8) [t-channel]'
c--  240 '  f(p1)+f(p2) --> t~(-->e-(p3)+nu~(p4)+b~(p5))+b(p6)+q(p7)+f(p8) [t-channel]'
c--      '  f(p1)+f(p2) --> t(no BR) + b~(p6) + q(p7) + f(p8)' (removebr=.true.)
        kcase=k4ftjet
        nqcdjets=4
      if (nproc == 235) then
        nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
      else
        nwz=-1
          plabel(3)='el'
          plabel(4)='na'
          plabel(5)='ba'
          plabel(6)='bq'
      endif
        plabel(7)='pp'
        plabel(8)='pp'
            
c--- ndim is one less than usual, since the top is always on-shell 
        ndim=15
        n3=1
        mass2=mt
        width2=twidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 345, (0,j=1,48) /)
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brtop
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          nqcdjets=3
        endif

c---  in the SM, the logical:: fourthgen should be false
c---  for BSM calculations, it should be true and it indicates that
c---   5 flavours should be used in the PDF and alpha-s
        fourthgen=.false.

      if (fourthgen) then
c--- BSM: full 5 light flavours
          nflav=5
        bmass=4.7_dp  !  set b-mass to its usual value
      else
c--- SM: only 4 light flavours
          nflav=4
          bmass=1001._dp !  enforce 4-flavour running in alfamz.f
        endif
      
c--- set up correct scales and as on heavy and light quark lines
        facscale_H=initfacscale_H
        facscale_L=initfacscale_L
        renscale_H=initrenscale_H
        renscale_L=initrenscale_L
c--- make sure it works even if not specifying separate scales
        if (initrenscale_L == zip) then
        facscale_H=facscale
        facscale_L=facscale
        renscale_H=scale
        renscale_L=scale
      endif
      
        b0=(xn*11._dp-2._dp*nflav)/6._dp
      as_H=alphas(abs(renscale_H),amz,nlooprun)
      as_L=alphas(abs(renscale_L),amz,nlooprun)
      
c-----------------------------------------------------------------------

      elseif ((nproc == 241) .or. (nproc == 246)
     &   .or. (nproc == 242) .or. (nproc == 247)) then
c--  241 '  f(p1)+f(p2) --> t(p3)+b~(p4)+f(p5) [s-channel]'
c--  246 '  f(p1)+f(p2) --> t~(p3)+b(p4)+f(p5) [s-channel]'

cc--- for comparison with C. Oleari's e+e- --> QQbg calculation
c      if (runstring(1:5) == 'carlo') then
cc---     heavy quark mass passed via chars 6 and 7 of runstring
c        read(runstring(6:7),67) imhq
c        mt=real(imhq,dp)
c        mb=mt
c          wmass=zip
c        write(6,*)
c        write(6,*) ' >>> HEAVY QUARK MASS = ',mt,' GeV <<<'
c      endif
c   67   format(i2)
      
      if     ((nproc == 241) .or. (nproc == 246)) then
          kcase=kqq_tbg
          nqcdjets=1
          ndim=7
      elseif ((nproc == 242) .or. (nproc == 247)) then
          kcase=kqqtbgg
          nqcdjets=2
          ndim=10
      else
        write(6,*) 'Unexpected value of nproc in chooser.f!'
        stop
      endif
        mass2=mt
        mass3=mb
      nflav=4
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
        plabel(6)='pp'
      if ((nproc == 241) .or. (nproc == 242)) then
        nwz=+1
      else
        nwz=-1
      endif

c--- set up correct scales and as on heavy and light quark lines
        facscale_H=initfacscale_H
        facscale_L=initfacscale_L
        renscale_H=initrenscale_H
        renscale_L=initrenscale_L
c--- make sure it works even if not specifying separate scales
        if (initrenscale_L == zip) then 
        facscale_H=facscale
        facscale_L=facscale
        renscale_H=scale
        renscale_L=scale
      endif
      
      bmass=1001._dp ! since nflav=4
        b0=(xn*11._dp-2._dp*nflav)/6._dp
      as_H=alphas(abs(renscale_H),amz,nlooprun)
      as_L=alphas(abs(renscale_L),amz,nlooprun)

      
c-----------------------------------------------------------------------

      elseif (nproc == 249) then
c--  e+ e^- -> 3 jets as check of massless limit of process 241

cc--- for comparison with C. Oleari's e+e- --> QQbg calculation
c      if (runstring(1:5) == 'carlo') then
cc---     heavy quark mass passed via chars 6 and 7 of runstring
c        read(runstring(6:7),67) imhq
c        mt=zip
c        mb=zip
c          wmass=zip
c        write(6,*)
c        write(6,*) ' >>> HEAVY QUARK MASS = ',mt,' GeV <<<'
c      endif
      
        kcase=kepem3j
        nqcdjets=1
        ndim=7

        mass2=zip
        mass3=zip
      nflav=4
        n2=0
        n3=0
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='pp'
        plabel(6)='pp'
      nwz=+1
      
      bmass=1001._dp ! since nflav=4
        b0=(xn*11._dp-2._dp*nflav)/6._dp
      
c-----------------------------------------------------------------------

      elseif (nproc == 251) then
c-- 251 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + W^+(-->nu(p5)+e^+(p6))+f(p7)+f(p8)'

        kcase=kWpWp2j
      nqcdjets=2
      ndim=16
      mb=zip
      plabel(3)='nl'
      plabel(4)='ea'
      plabel(5)='nl'
      plabel(6)='ea'
      plabel(7)='pp'
      plabel(8)='pp'
      plabel(9)='pp'
      l1=1._dp

      n2=1
      n3=1
      mass2=wmass
      width2=wwidth
      mass3=wmass
      width3=wwidth

      mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)
        
      write(*,*)'Setting zerowidth to true for process 131'
      zerowidth = .true.
      write(*,*)'Setting removebr to false for process 131'
      removebr = .false.


c-----------------------------------------------------------------------

      elseif (nproc == 252) then
c-- 252 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + W^+(-->nu(p5)+e^+(p6))+f(p7)+f(p8)+f(p9)'

        kcase=kWpWp3j
      nqcdjets=3
      ndim=19
      mb=zip
      plabel(3)='nl'
      plabel(4)='ea'
      plabel(5)='nl'
      plabel(6)='ea'
      plabel(7)='pp'
      plabel(8)='pp'
      plabel(9)='pp'
      plabel(10)='pp'
      l1=1._dp

      n2=1
      n3=1
      mass2=wmass
      width2=wwidth
      mass3=wmass
      width3=wwidth
        
      mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)
        
      write(*,*)'Setting zerowidth to true for process 132'
      zerowidth = .true.
      write(*,*)'Setting removebr to false for process 132'
      removebr = .false.


c-----------------------------------------------------------------------
      elseif ((nproc == 253) .or. (nproc == 254)) then
c-- 253 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + Z(-->e^-(p5)+e^+(p6))+f(p7)+f(p8)'

        kcase=kWpmZjj
        call checkminzmass(2)
      nqcdjets=2
      ndim=16
      mb=zip

        if (nproc == 253) then
        nwz=1
      plabel(3)='nl'
      plabel(4)='ea'
      plabel(5)='el'
      plabel(6)='ea'
      else
      nwz=-1
      plabel(3)='el'
      plabel(4)='na'
      plabel(5)='el'
      plabel(6)='ea'
      endif
      
      plabel(7)='pp'
      plabel(8)='pp'

      n2=1
      n3=1
      mass2=zmass
      width2=zwidth
      mass3=wmass
      width3=wwidth

      mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)
        
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brzee
        endif


c-----------------------------------------------------------------------
      elseif (nproc == 255) then
c-- 255 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + Z(-->e^-(p5)+e^+(p6))+b(p7)+f(p8)'

        kcase=kWpmZbj
        call checkminzmass(2)
      nqcdjets=2
        nwz=1
      ndim=16
      mb=zip
      plabel(3)='nl'
      plabel(4)='ea'
      plabel(5)='el'
      plabel(6)='ea'
      plabel(7)='bq'
      plabel(8)='pp'

      n2=1
      n3=1
      mass2=zmass
      width2=zwidth
      mass3=wmass
      width3=wwidth

      mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)
        
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen*brzee
        endif


c-----------------------------------------------------------------------

      elseif (nproc == 256) then
c-- 256 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + Z(-->e^-(p5)+e^+(p6))+b(p7)+f(p8)'

        kcase=kWpmZbj
        call checkminzmass(2)
      nqcdjets=2
      ndim=16
        nwz=-1
      mb=zip
      plabel(3)='el'
      plabel(4)='na'
      plabel(5)='el'
      plabel(6)='ea'
      plabel(7)='bq'
      plabel(8)='pp'

      n2=1
      n3=1
      mass2=zmass
      width2=zwidth
      mass3=wmass
      width3=wwidth

      mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)
        
c--- total cross-section             
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              plabel(5)='ig'
              plabel(6)='ig'
              call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
              BrnRat=brwen*brzee
            endif

c-----------------------------------------------------------------------
      elseif (nproc == 259) then
c-- 259 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) + Z(-->e^-(p5)+e^+(p6))+b(p7)+b~(p8)'

        kcase=kWpmZbb
        call checkminzmass(2)
      nqcdjets=2
        nwz=1
      ndim=16
      mb=zip
      plabel(3)='nl'
      plabel(4)='ea'
      plabel(5)='el'
      plabel(6)='ea'
      plabel(7)='bq'
      plabel(8)='ba'

      n2=1
      n3=1
      mass2=zmass
      width2=zwidth
      mass3=wmass
      width3=wwidth

      mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)
        
c--- total cross-section             
       if (removebr) then
         plabel(3)='ig'
         plabel(4)='ig'
         plabel(5)='ig'
         plabel(6)='ig'
         call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
         BrnRat=brwen*brzee
       endif

c-----------------------------------------------------------------------

      elseif (nproc == 260) then
c--260 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+ Z(-->e^-(p5)+e^+(p6))+b(p7)+b~(p8)'

        kcase=kWpmZbb
        call checkminzmass(2)
      nqcdjets=2
      ndim=16
        nwz=-1
      mb=zip
      plabel(3)='el'
      plabel(4)='na'
      plabel(5)='el'
      plabel(6)='ea'
      plabel(7)='bq'
      plabel(8)='ba'

      n2=1
      n3=1
      mass2=zmass
      width2=zwidth
      mass3=wmass
      width3=wwidth

      mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)
        
c--- total cross-section             
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              plabel(5)='ig'
              plabel(6)='ig'
              call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
              BrnRat=brwen*brzee
            endif
c-----------------------------------------------------------------------


      elseif ((nproc == 261) .or. (nproc == 266)) then
c--  261 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)'
c--  266 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)[+b~(p6)]'
        kcase=kgQ__ZQ
        nqcdjets=1
        flav=5
        nwz=0
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        isub=1+(nproc-261)/5
        if (nproc == 261) then
          plabel(6)='pp'
        else
          plabel(6)='ba'
        endif
        q1=-1._dp
        l1=le
        r1=re

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
      elseif ((nproc == 262) .or. (nproc == 267)) then
c--  262 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c(p5)'
c--  267 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c(p5)[+c~(p6)]'
        kcase=kgQ__ZQ
        nqcdjets=1
        flav=4
        nwz=0
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
        
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        isub=1+(nproc-262)/5
        if (nproc == 262) then
          plabel(6)='pp'
        else
          plabel(6)='ba'
        endif
        q1=-1._dp
        l1=le
        r1=re
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
      elseif (nproc == 263) then
c--  263 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b~(p5)+b(p6) (1 b-tag)'
        kcase=kZbbmas
        nqcdjets=2
        notag=1

        ndim=10
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        write(6,*) 'mb=',mb
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        q1=-1._dp
        l1=le
        r1=re
        
      flav=5
      nflav=4

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
      elseif (nproc == 264) then
c--  264 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c~(p5)+c(p6) (1 c-tag)'
        kcase=kZccmas
        nqcdjets=2
        notag=1
        
        ndim=10
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth
        
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        mb=mc
        write(6,*) 'mc=',mb
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        q1=-1._dp
        l1=le
        r1=re

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif      
           
c-----------------------------------------------------------------------
          
      elseif ((nproc == 270) .or. (nproc == 271) 
     &   .or. (nproc == 272)) then
      
c--- turn off Higgs decay, for speed
c        nodecay=.true.      
c--- parameters to turn off various pieces, for checking
        f0q=one
        f2q=one
        f4q=one
      
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=10
      
        n2=0
        n3=1

        mass3=hmass
        width3=hwidth
        
        mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)
        
        if     (nproc == 270) then
c-- 270 '  f(p1)+f(p2) --> H(gamma(p3)+gamma(p4))+f(p5)+f(p6)[in heavy top limit]'
c--     '  f(p1)+f(p2) --> H(no BR)+f(p5)+f(p6)[in heavy top limit]' (removebr=.true.)
          hdecaymode='gaga'
          plabel(3)='ga'
          plabel(4)='ga'
          kcase=kgagajj
          nqcdjets=2
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            BrnRat=gamgambr
          endif
          
        elseif     (nproc == 271) then
c-- 271 '  f(p1)+f(p2) --> H(b(p3)+b~(p4))+f(p5)+f(p6)[in heavy top limit]'
c--     '  f(p1)+f(p2) --> H(no BR)+f(p5)+f(p6)[in heavy top limit]' (removebr=.true.)
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'
          kcase=kggfus2
          nqcdjets=4
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=2
            BrnRat=br
          endif
          
        elseif (nproc == 272) then
c-- 272 '  f(p1)+f(p2) --> H(tau-(p3)+tau+(p4))+f(p5)+f(p6)[in heavy top limit]'
c--     '  f(p1)+f(p2) --> H(no BR)+f(p5)+f(p6)[in heavy top limit]' (removebr=.true.)
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'
          kcase=kggfus2
          nqcdjets=2
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            Brnrat=tautaubr
          endif
        endif
                
c-----------------------------------------------------------------------

      elseif     (nproc == 273) then
c-- 273 '  f(p1)+f(p2) -->` H(-->W^+(p3,p4)W^-(p5,p6)) + f(p7) + f(p8)'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

c--- parameters to turn off various pieces, for checking
        f0q=one
        f2q=one
        f4q=one
      
        kcase=kHWW2jt
        ndim=16
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='pp'
        plabel(8)='pp'
        plabel(9)='pp'
        nqcdjets=2
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)
        
c--- print warning if we're below threshold
        if (hmass < 2._dp*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
      if (removebr) then
      write(6,*)
      write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
        stop
      endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif
        
        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=wwbr*brwen**2
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif

c-----------------------------------------------------------------------

      elseif     (nproc == 274) then
c-- 274 f(p1)+f(p2)->H(Z^+(e^-(p3),e^+(p4))Z(mu^-(p5),mu^+(p6)))+f(p7)+f(p8)
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

c--- parameters to turn off various pieces, for checking
        f0q=one
        f2q=one
        f4q=one
      
        kcase=kHZZ2jt
        l1=le
        r1=re
        l2=le
        r2=re
        ndim=16
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='ml'
        plabel(6)='ma'
        plabel(7)='pp'
        plabel(8)='pp'
        plabel(9)='pp'
        nqcdjets=2
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)
        
c--- print warning if we're below threshold
        if (hmass < 2._dp*zmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
      if (removebr) then
      write(6,*)
      write(6,*) 'Cannot remove H->ZZ BR, not defined below threshold'
        stop
      endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif
        
        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=2._dp*brzee**2*zzbr  ! factor of 2 for identical particles
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif

c-----------------------------------------------------------------------

      elseif ((nproc == 275) .or. (nproc == 276)) then

c--- parameters to turn off various pieces, for checking
        f0q=one
        f2q=one
        f4q=one

        kcase=kggfus3
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(5)='pp'
        plabel(6)='pp'
        plabel(7)='pp'
        ndim=13
      
        n2=0
        n3=1

        mass3=hmass
        width3=hwidth
        
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 275) then
c-- 275 '  f(p1)+f(p2) --> H(b(p3)+b~(p4))+f(p5)+f(p6)+f(p7)[in heavy top limit]'
c--     '  f(p1)+f(p2) --> H(no BR)+f(p5)+f(p6)+f(p7)[in heavy top limit]' (removebr=.true.)
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'
          nqcdjets=5
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            nqcdjets=3
            BrnRat=br
          endif
          
        elseif (nproc == 276) then
c-- 276 '  f(p1)+f(p2) --> H(tau-(p3)+tau+(p4))+f(p5)+f(p6)+f(p7)[in heavy top limit]'
c--     '  f(p1)+f(p2) --> H(no BR)+f(p5)+f(p6)+f(p7)[in heavy top limit]' (removebr=.true.)
          hdecaymode='tlta'
          plabel(3)='tl'
          plabel(4)='ta'
          nqcdjets=3
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            BrnRat=tautaubr
          endif
        endif

c-----------------------------------------------------------------------

      elseif     (nproc == 278) then
c-- 278 '  f(p1)+f(p2) --> H(-->W^+(p3,p4)W^-(p5,p6)) + f(p7) + f(p8) + f(p9)'
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

c--- parameters to turn off various pieces, for checking
        f0q=one
        f2q=one
        f4q=one
      
        kcase=kHWW3jt
        ndim=19
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='pp'
        plabel(8)='pp'
        plabel(9)='pp'
        nqcdjets=3
        n2=1
        n3=1
        mass2=wmass
        width2=wwidth
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)
        
c--- print warning if we're below threshold
        if (hmass < 2._dp*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
      if (removebr) then
      write(6,*)
      write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
        stop
      endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif
        
        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=wwbr*brwen**2
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif

c-----------------------------------------------------------------------

      elseif     (nproc == 279) then
c-- 279 f(p1)+f(p2)->H(Z^+(e^-(p3),e^+(p4))Z(mu^-(p5),mu^+(p6)))+f(p7)+f(p8)+f(p9)
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)

c--- parameters to turn off various pieces, for checking
        f0q=one
        f2q=one
        f4q=one
      
        kcase=kHZZ3jt
        l1=le
        r1=re
        l2=le
        r2=re
        ndim=19
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='ml'
        plabel(6)='ma'
        plabel(7)='pp'
        plabel(8)='pp'
        plabel(9)='pp'
        nqcdjets=3
        n2=1
        n3=1
        mass2=zmass
        width2=zwidth
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)
        
c--- print warning if we're below threshold
        if (hmass < 2._dp*zmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->ZZ is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
      if (removebr) then
      write(6,*)
      write(6,*) 'Cannot remove H->ZZ BR, not defined below threshold'
        stop
      endif
        if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif
        
        if (removebr) then
        call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
        BrnRat=2._dp*brzee**2*zzbr  ! factor of 2 for identical particles
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='ig'
        endif

c-----------------------------------------------------------------------

      elseif (nproc == 280) then
             ndim=4 
             kcase=kdirgam
             plabel(3)='ga'
             plabel(4)='pp'
             plabel(5)='pp'
             lastphot=3
             nqcdjets=1
             n3=0
             inclusive=.true.
           write(6,*)
           write(6,*) 'Setting inclusive = .true. '//
     &                  'for direct photon production.'

      elseif (nproc == 282) then
             ndim=7 
             kcase=kgamjet
             plabel(3)='ga'
             plabel(4)='pp'
             plabel(5)='pp'
             plabel(6)='pp'
             lastphot=3
             nqcdjets=2
             n3=0

      elseif (nproc == 283) then
             ndim=4 
             kcase=khflgam
             flav=5
             plabel(3)='ga'
             plabel(4)='bq'
             plabel(5)='pp'
             lastphot=3
             nqcdjets=1
             n3=0
             inclusive=.true.
           
      elseif (nproc == 284) then
             ndim=4 
             kcase=khflgam
             flav=4
             plabel(3)='ga'
             plabel(4)='bq'
             plabel(5)='pp'
             lastphot=3
             nqcdjets=1
             n3=0
             inclusive=.true.
                
      elseif (nproc == 285) then
!==========reset esq^2 so alpha_EM  = 1/137
             reset_alphaEM=.true.
             ndim=4 
             kcase=kgamgam
             plabel(3)='ga'
             plabel(4)='ga'
             plabel(5)='pp'
             lastphot=4
             nqcdjets=0
             n3=0

      elseif (nproc == 2851) then
!==========reset esq^2 so alpha_EM  = 1/137
             reset_alphaEM=.true.
             ndim=4 
             toploopgaga=.false.
             kcase=kgg2gam
c--- this process works best using the new PS generation
             new_pspace=.true.
             plabel(3)='ga'
             plabel(4)='ga'
             plabel(5)='pp'
             lastphot=4
             nqcdjets=0
             n3=0
          elseif (nproc == 2852) then
!===== incl top loops by default
!==========reset esq^2 so alpha_EM  = 1/137
             reset_alphaEM=.true.
             ndim=4
             toploopgaga=.true.
             
             kcase=kgg2gam
c--- this process works best using the new PS generation
             new_pspace=.true.
             plabel(3)='ga'
             plabel(4)='ga'
             plabel(5)='pp'
             lastphot=4
             nqcdjets=0
             n3=0

      elseif (nproc == 286) then
!==========reset esq^2 so alpha_EM  = 1/137
             reset_alphaEM=.true.
             ndim=7 
             kcase=kgmgmjt
c--- this process works best using the new PS generation
c             new_pspace=.true.
c--- this process works best using the new PS generation
             plabel(3)='ga'
             plabel(4)='ga'
             plabel(5)='pp'
             plabel(6)='pp'
             lastphot=4
             nqcdjets=1
             n3=0
     
      elseif (nproc == 287) then
             ndim=7 
             kcase=ktrigam
c--- this process works best using the new PS generation
             new_pspace=.true.
c--- this process works best using the new PS generation
             plabel(3)='ga'
             plabel(4)='ga'
             plabel(5)='ga'
             plabel(6)='pp'
             lastphot=5
             nqcdjets=0
             n3=0

      elseif (nproc == 288) then
             ndim=10 
             kcase=kgmgmjj
             plabel(3)='ga'
             plabel(4)='ga'
             plabel(5)='pp'
             plabel(6)='pp'
             plabel(7)='pp' 
             lastphot=4
             nqcdjets=2
             n3=0
                
      elseif (nproc == 289) then
             ndim=10
             kcase=kfourga
c--- this process works best using the new PS generation
             new_pspace=.true.
             plabel(3)='ga'
             plabel(4)='ga'
             plabel(5)='ga'
             plabel(6)='ga'
             plabel(7)='pp'
             lastphot=6
             nqcdjets=0
             n3=0

c-----------------------------------------------------------------------

c--- These two processes need to be moved to other numbers
      elseif (nproc == 9280) then
c--  280      '  f(p1)+f(p2)--> f(p3)+f(p4)'
             ndim=4 
             kcase=ktwojet
             plabel(3)='pp'
             plabel(4)='pp'
             plabel(5)='pp'
             nqcdjets=2
             n3=0
                
             mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
      elseif (nproc == 9281) then
c--  281      '  f(p1)+f(p2)--> f(p3)+f(p4)+f(p5)'
             ndim=7 
             kcase=kthrjet
             plabel(3)='pp'
             plabel(4)='pp'
             plabel(5)='pp'
             plabel(6)='pp'
             nqcdjets=3
             n3=0
                
c-----------------------------------------------------------------------

      elseif ((nproc == 290) .or. (nproc == 295)) then
        kcase=kWgamma
        nqcdjets=0
        ndim=7
        mb=0
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
        plabel(5)='ga'
        plabel(6)='pp'
        lastphot=5

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 290) then
c-- 290 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+gamma(p5)'
c--     '  f(p1)+f(p2) --> W^+ (no BR) + gamma(p5)' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 295) then
c-- 295 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~+(p4))+gamma(p5)'
c--     '  f(p1)+f(p2) --> W^- (no BR) + gamma(p5)' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
      
        if (zerowidth .eqv. .false.) then
        write(6,*)
          write(6,*) 'Setting removebr to .false. in order to ensure'
        write(6,*) 'lepton-photon singularity can be removed'
        removebr=.false.
      endif
      
c--- total cross-section
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc == 292) .or. (nproc == 297)) then
       
         kcase=kWgajet
         nqcdjets=1
         ndim=10
         mb=0
         rescale=.false.
         n2=0
         n3=1
         mass3=wmass
         width3=wwidth
         plabel(5)='ga'
         plabel(6)='pp'
         lastphot=5
        
         mcfmplotinfo= (/ 34, (0,j=1,49) /)
                 
        if     (nproc == 292) then
c 292 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+(f(p5) --> gamma(p5))'
c     '  f(p1)+f(p2) --> W^+(No BR)+(f(p5) --> gamma(p5)) (removebr =.true.)'
           nwz=1
           plabel(3)='nl'
           plabel(4)='ea'
         
        elseif (nproc == 297) then
c 297 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+(f(p5) -->gamma(p5))'
c     '  f(p1)+f(p2) --> W^-(no BR)+(f(p5) -->gamma(p5)) (removebr=.true.)'
           nwz=-1
           plabel(3)='el'
           plabel(4)='na'
       endif
        
       if (zerowidth .eqv. .false.) then
       write(6,*)
         write(6,*) 'Setting removebr to .false. in order to ensure'
       write(6,*) 'lepton-photon singularity can be removed'
       removebr=.false.
       endif
      
c---  total cross-section             
       if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
       endif

      
c--------------------------------------------------------------------------------------------------

        elseif ((nproc == 300) .or. (nproc == 305)) then
          kcase=kZgamma
          nqcdjets=0
          ndim=7
          n2=0
          n3=1
          mass3=zmass
          width3=zwidth
          nwz=0
          plabel(5)='ga'
          plabel(6)='pp'
          lastphot=5
          
          mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
          if     (nproc == 300) then
c-- 300 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+gamma(p5)'
c--     '  f(p1)+f(p2) --> Z^0 (no BR) +gamma(p5)' (removebr=.true.)
            call checkminzmass(1)
            plabel(3)='el'
            plabel(4)='ea'
            q1=-1._dp
            l1=le
            r1=re
            if (zerowidth .eqv. .false.) then
            write(6,*)
              write(6,*)'Setting removebr to .false. in order to ensure'
            write(6,*)'lepton-photon singularity can be removed'
            removebr=.false.
          endif
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
              BrnRat=brzee
            endif
          elseif (nproc == 305) then
c-- 305 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4)))-(sum over 3 nu)+gamma(p5)'
            plabel(3)='nl'
            plabel(4)='na'
            q1=zip
            l1=ln*sqrt(3._dp)
            r1=rn*sqrt(3._dp)
          endif

c-----------------------------------------------------------------------

        elseif ((nproc == 302) .or. (nproc == 307)) then
          kcase=kZgajet
          nqcdjets=1
          ndim=10
          n2=0
          n3=0
          mass3=zmass
          width3=zwidth
          nwz=0
          plabel(5)='ga'
          plabel(6)='pp'
          plabel(7)='pp'
          lastphot=5
          
          mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
          if     (nproc == 302) then
c-- 302 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))+gamma(p5)+f(p6)'
c--     '  f(p1)+f(p2) --> Z^0 (no BR) +gamma(p5)+jet(p6)' (removebr=.true.)
            call checkminzmass(1)
            plabel(3)='el'
            plabel(4)='ea'
            q1=-1._dp
            l1=le
            r1=re
            if (zerowidth .eqv. .false.) then
            write(6,*)
              write(6,*)'Setting removebr to .false. in order to ensure'
            write(6,*)'lepton-photon singularity can be removed'
            removebr=.false.
          endif
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
              BrnRat=brzee
            endif
          elseif  (nproc == 307) then
c-- 307 '  f(p1)+f(p2) --> Z^0(-->nu(p3)+nu~(p4))+gamma(p5)+f(p6)'
            plabel(3)='nl'
            plabel(4)='na'
            q1=zip
            l1=ln*sqrt(3._dp)
            r1=rn*sqrt(3._dp)
            if (removebr) then
              plabel(3)='ig'
              plabel(4)='ig'
              call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
              BrnRat=brznn
            endif
          endif
c-----------------------------------------------------------------------

      elseif ((nproc == 301) .or. (nproc == 306)) then
c-- 301 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4))+gamma(p5)+gamma(p6)'
c-- 306 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))+gamma(p5)+gamma(p6)'
        kcase=kZ_2gam
        ndim=10
        n2=0
        n3=0

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if (nproc == 301) then
           call checkminzmass(1)
           plabel(3)='el'
           plabel(4)='ea'
           q1=-1._dp
           l1=le
           r1=re
           if (zerowidth .eqv. .false.) then
           write(6,*)
             write(6,*)'Setting removebr to .false. in order to ensure'
           write(6,*)'lepton-photon singularity can be removed'
           removebr=.false.
         endif
c--- total cross-section             
           if (removebr) then
             plabel(3)='ig'
             plabel(4)='ig'
             call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
             BrnRat=brzee
           endif
        elseif (nproc == 306) then
           plabel(3)='nl'
           plabel(4)='na'
           q1=zip
           l1=ln*sqrt(3._dp)
           r1=rn*sqrt(3._dp)
        endif
        plabel(5)='ga'
        plabel(6)='ga'
        plabel(7)='pp'
        lastphot=6
        nwz=0   
        mass3=zmass
        width3=zwidth
        
c-----------------------------------------------------------------------

      elseif ((nproc == 303) .or. (nproc == 308)) then
c-- 303 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4))+gamma(p5)+gamma(p6)+f(p7)'
c-- 308 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))+gamma(p5)+gamma(p6)+f(p7)'
        kcase=kZ2gajt
        nqcdjets=1
        ndim=13
        n2=0
        n3=0

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if (nproc == 303) then
           call checkminzmass(1)
           plabel(3)='el'
           plabel(4)='ea'
           q1=-1._dp
           l1=le
           r1=re
           if (zerowidth .eqv. .false.) then
           write(6,*)
             write(6,*)'Setting removebr to .false. in order to ensure'
           write(6,*)'lepton-photon singularity can be removed'
           removebr=.false.
         endif
c--- total cross-section             
           if (removebr) then
             plabel(3)='ig'
             plabel(4)='ig'
             call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
             BrnRat=brzee
           endif
        elseif (nproc == 308) then
           plabel(3)='nl'
           plabel(4)='na'
           q1=zip
           l1=ln*sqrt(3._dp)
           r1=rn*sqrt(3._dp)
        endif
        plabel(5)='ga'
        plabel(6)='ga'
        plabel(7)='pp'
        lastphot=6
        nwz=0   
        mass3=zmass
        width3=zwidth
                 
c-----------------------------------------------------------------------

      elseif ((nproc == 304) .or. (nproc == 309)) then
c-- 304 '  f(p1)+f(p2) --> Z^0(e^-(p3)+e^+(p4))+gamma(p5)+f(p6)+f(p7)'
c-- 309 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))+gamma(p5)+f(p6)+f(p7)'
        kcase=kZga2jt
        nqcdjets=2
        ndim=13
        n2=0
        n3=0

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if (nproc == 304) then
           call checkminzmass(1)
           plabel(3)='el'
           plabel(4)='ea'
           q1=-1._dp
           l1=le
           r1=re
           if (zerowidth .eqv. .false.) then
           write(6,*)
             write(6,*)'Setting removebr to .false. in order to ensure'
           write(6,*)'lepton-photon singularity can be removed'
           removebr=.false.
         endif
c--- total cross-section             
           if (removebr) then
             plabel(3)='ig'
             plabel(4)='ig'
             call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
             BrnRat=brzee
           endif
        elseif (nproc == 309) then
           plabel(3)='nl'
           plabel(4)='na'
           q1=zip
           l1=ln*sqrt(3._dp)
           r1=rn*sqrt(3._dp)
        endif
        plabel(5)='ga'
        plabel(6)='pp'
        plabel(7)='pp'
        lastphot=5
        nwz=0   
        mass3=zmass
        width3=zwidth
                         
c-----------------------------------------------------------------------

      elseif ((nproc == 311) .or. (nproc == 316)) then
        kcase=kW_bjet
        nqcdjets=2
        flav=5
        isub=1
        
        nflav=5
        mb=zip
        plabel(5)='bq'
        plabel(6)='pp'
        plabel(7)='pp'
        
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 311) then
c--  311 '  f(p1)+b(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+f(p6)'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 316) then
c--  316 '  f(p1)+b(p2) --> W^-(-->e^-(p3)+nu~(p4))+b(p5)+f(p6)'
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
        
        if (removebr) then
c--      '  f(p1)+b(p2) --> W(no BR)+b(p5)+f(p6)' (removebr=.true.)
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc == 321) .or. (nproc == 326)) then
        kcase=kW_bjet
        nqcdjets=2
        flav=4
        isub=1
        
        nflav=4
        mb=zip
        plabel(5)='bq'
        plabel(6)='pp'
        plabel(7)='pp'
        
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 321) then
c--  321 '  f(p1)+b(p2) --> W^+(-->nu(p3)+e^+(p4))+c(p5)+f(p6)'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 326) then
c--  326 '  f(p1)+b(p2) --> W^-(-->e^-(p3)+nu~(p4))+c(p5)+f(p6)'
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
        
        if (removebr) then
c--      '  f(p1)+b(p2) --> W(no BR)+c(p5)+f(p6)' (removebr=.true.)
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc == 331) .or. (nproc == 336)) then
        kcase=kWcjetg
        nqcdjets=2
        nflav=3
        
        plabel(5)='bq'
        plabel(6)='pp'

        ndim=10
        mb=0
        n2=0
        n3=1
        mass2=zip
        mass3=wmass
        width3=wwidth
        
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 331) then
c--  331 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+c(p5)+f(p6) [c-s interaction]'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 336) then
c--  336 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+c(p5)+f(p6) [c-s interaction]'
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
        
        if (removebr) then
c--      '  f(p1)+f(p2) --> W(no BR)+c(p5)+f(p6) [c-s interaction]' (removebr=.true.)
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc == 341) .or. (nproc == 351) 
     &   .or. (nproc == 342) .or. (nproc == 352)) then
        kcase=kZ_bjet
        call checkminzmass(1)
        ndim=10
        n2=0
        n3=1
        nqcdjets=2
        
        mb=zip
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     ((nproc == 341) .or. (nproc == 351)) then
          isub=1        
          plabel(6)='pp'
          plabel(7)='pp'
        elseif ((nproc == 342) .or. (nproc == 352)) then
          isub=2        
          plabel(6)='ba'
          plabel(7)='pp'
        endif
        
        q1=-1._dp
        l1=le
        r1=re
        nwz=0   
        mass3=zmass
        width3=zwidth

        if     ((nproc == 341) .or. (nproc == 342)) then
c--  341 '  f(p1)+b(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+f(p6)'
          flav=5
          nflav=5
        elseif ((nproc == 351) .or. (nproc == 352)) then
c--  351 '  f(p1)+c(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c(p5)+f(p6)'
          flav=4
          nflav=4
        endif
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc == 346) .or. (nproc == 356) 
     &   .or. (nproc == 347) .or. (nproc == 357)) then
        kcase=kZbjetg
        call checkminzmass(1)
        ndim=13
        n2=0
        n3=1
        nqcdjets=3
        
        mb=zip
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='bq'
        
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     ((nproc == 346) .or. (nproc == 356)) then
c--  346 '  f(p1)+b(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+f(p6)+f(p7)'
          isub=1        
          plabel(6)='pp'
          plabel(7)='pp'
        elseif ((nproc == 347) .or. (nproc == 357)) then
c--  347 '  f(p1)+b(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+f(p6)+b~(p7)'
          isub=2        
          plabel(6)='ba'
          plabel(7)='pp'
        endif
        
        q1=-1._dp
        l1=le
        r1=re
        nwz=0   
        mass3=zmass
        width3=zwidth

        if     ((nproc == 346) .or. (nproc == 347)) then
c--  346 '  f(p1)+b(p2) --> Z^0(-->e^-(p3)+e^+(p4))+b(p5)+f(p6)+f(p7)'
          flav=5
          nflav=5
        elseif ((nproc == 356) .or. (nproc == 357)) then
c--  356 '  f(p1)+c(p2) --> Z^0(-->e^-(p3)+e^+(p4))+c(p5)+f(p6)+f(p7)'
          flav=4
          nflav=4
        endif
        
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif (nproc/10 == 36) then
        kcase=kW_only
        n3=1
        ndim=4
        nqcdjets=0
C-- 361 '  c(p1)+sbar(p2) --> W^+(-->nu(p3)+e^+(p4))'
c--- This can be used to calculate cs->W at NLO, in the ACOT-like
c--- scheme where the charm quark appears in the initial state but
c--- the real corrections use a massless charm quark in the final state
c--- (c.f. processes 362 and 363 below)
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='pp'
        nwz=1

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
c--- total cross-section
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
      
c--- change W mass (after couplings and BRs already calculated)
c      if     (runstring(4:8) == 'mw_80') then
c        wmass=80.4_dp
c      elseif (runstring(4:8) == 'mw200') then
c        wmass=200._dp
c      elseif (runstring(4:8) == 'mw400') then
c        wmass=400._dp
c      endif
      
c--- change charm mass
c      if     (runstring(9:13) == 'mc1.3') then
c        mc=1.3_dp
c        mcsq=mc**2
c      elseif (runstring(9:13) == 'mc4.5') then
c        mc=4.5_dp
c        mcsq=mc**2
c      elseif (runstring(9:13) == 'mc20.') then
c        mc=20._dp
c        mcsq=mc**2
c      endif
      
        mass3=wmass
        width3=wwidth
c--- set CKM matrix to remove all elements except for Vcs
        Vud=0._dp
        Vus=zip
        Vub=zip
        Vcd=zip
        Vcs=1._dp
        Vcb=zip

c--- To obtain a complete prediction for cs->W at NLO, including the
c---  the effect of the charm quark mass (as one should, according to
c---  ACOT) processes 362 and 363 should be summed (using ktota)

c--- real matrix elements W+s and W+g, with corresponding (massless)
c---  integrated counterterms and virtual matrix elements
      if (nproc == 362) kcase=kWcsbar

c--- real matrix elements W+c including the charm quark mass,
c---  with the virtual contribution representing the counterterm
c---  consisting of the logarithm convolution      
      if (nproc == 363) kcase=kWcs_ms

        if ( (kpart==klord) .and.  
     &       ((nproc == 362) .or. (nproc == 363)) ) then
          write(6,*) 'This process number is not suitable for the'
          write(6,*) 'LO calculation. Please run process 361'
          write(6,*) 'for the corresponding Born contribution.'
          stop
        endif

c-----------------------------------------------------------------------
      elseif ((nproc == 370) .or. (nproc == 371)) then
c-- 370 '  f(p1)+f(p2) --> W^+(nu(p3)+e^+(p4))+gamma(p5)+gamma(p6)'
c-- 371 '  f(p1)+f(p2) --> W^-(e^-(p3)+nu~(p4))+gamma(p5)+gamma(p6)'
        kcase=kW_2gam
        ndim=10
        n2=0
        n3=0
        lastphot=6
        mcfmplotinfo= (/ 34, 345, 346, 3456, 56, (0,j=1,45) /)
        
        if (nproc == 370) then
           nwz=1
           plabel(3)='nl'
           plabel(4)='ea'
           if (zerowidth .eqv. .false.) then
           write(6,*)
           write(6,*)'Setting removebr to .false. in order to ensure'
           write(6,*)'lepton-photon singularity can be removed'
           removebr=.false.
           endif
c--- total cross-section             
           if (removebr) then
             plabel(3)='ig'
             plabel(4)='ig'
             call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
             BrnRat=brwen
           endif
        elseif (nproc == 371) then
           plabel(3)='el'
           plabel(4)='na'
           nwz=-1
           if (zerowidth .eqv. .false.) then
           write(6,*)
           write(6,*)'Setting removebr to .false. in order to ensure'
           write(6,*)'lepton-photon singularity can be removed'
           removebr=.false.
           endif  
           if (removebr) then
             plabel(3)='ig'
             plabel(4)='ig'
             call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
             BrnRat=brwen
           endif
        endif
        plabel(5)='ga'
        plabel(6)='ga'
        plabel(7)='pp'
        mass3=wmass
        width3=wwidth
        
c-----------------------------------------------------------------------


      elseif ((nproc >= 401) .and. (nproc <= 408)) then
        kcase=kWbbmas
        write(6,*) 'mb=',mb
        flav=5

        bbproc=.false.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        nqcdjets=2
        notag=1
      
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth           

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc <= 403) then
c--- 401  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5) [massive]'
c---      '  f(p1)+f(p2) --> W^+ (no BR) +b(p5) [massive]' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc >= 406) then
c--- 406  '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + b(p5) [massive]'
c---      '  f(p1)+f(p2) --> W^- (no BR) +b(p5) [massive]' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
 
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif

c-----------------------------------------------------------------------

      elseif ((nproc == 411) .or. (nproc == 416)) then
        kcase=kW_bjet
        nqcdjets=2
        flav=5
        isub=1
        
c--- check for Wb+X flag and allow one jet to be untagged in that case
        notag=1
        write(6,*)
        write(6,*)'****************************************************'
        write(6,*)'* WARNING: cuts allow final state of Wb+X          *'
        write(6,*)'****************************************************'
      
        nflav=5
        mb=zip
        plabel(5)='bq'
        plabel(6)='pp'
        plabel(7)='pp'
        
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 411) then
c--  411 '  f(p1)+b(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5)+f(p6)'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 416) then
c--  416 '  f(p1)+b(p2) --> W^-(-->e^-(p3)+nu~(p4))+b(p5)+f(p6)'
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
        
        if (removebr) then
c--      '  f(p1)+b(p2) --> W(no BR)+b(p5)+f(p6)' (removebr=.true.)
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
          plabel(3)='ig'
          plabel(4)='ig'
        endif
        
c-----------------------------------------------------------------------

      elseif ((nproc == 421) .or. (nproc == 426)) then
        kcase=kWbbmas
        write(6,*) 'mb=',mb
        flav=5

        bbproc=.false.
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        nqcdjets=2
        notag=1
      
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth           

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 421) then
c--- 421  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+b(p5) [massive]'
c---      '  f(p1)+f(p2) --> W^+ (no BR) +b(p5) [massive]' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc >= 426) then
c--- 426  '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) + b(p5) [massive]'
c---      '  f(p1)+f(p2) --> W^- (no BR) +b(p5) [massive]' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
 
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif

c-----------------------------------------------------------------------

      elseif ((nproc == 431) .or. (nproc == 436)) then
        kcase=kWbbjem
        write(6,*) 'mb=',mb
        nqcdjets=3
        flav=5
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='pp'
        ndim=13
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
 
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 431) then
c-- 431 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +b(p5)+b~(p6)+f(p7) [massive]'
c--     '  f(p1)+f(p2) --> W^+ (no BR) +b(p5)+b~(p6)+f(p7) [massive]' (removebr=.true.)
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 436) then
c-- 436 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4)) +b(p5)+b~(p6)+f(p7) [massive]'
c--     '  f(p1)+f(p2) --> W^- (no BR) +b(p5)+b~(p6)+f(p7) [massive]' (removebr=.true.)
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
 
c--- total cross-section           
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
           
c-----------------------------------------------------------------------

      elseif ((nproc == 500) .or. (nproc == 510)) then
        kcase=kWttmas
        write(6,*) 'mt=',mt
        flav=6
        plabel(5)='ig'
        plabel(6)='ig'
        plabel(7)='pp'
        ndim=10
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth
 
        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 500) then
C-- 500 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +t(p5)+t~(p6) [massive]' 'N'
          nwz=1
          plabel(3)='nl'
          plabel(4)='ea'
        elseif (nproc == 510) then
C-- 510 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+t(p5)+t~(p6) [massive]' 'N'
          nwz=-1
          plabel(3)='el'
          plabel(4)='na'
        endif
 
c--- total cross-section             
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brwen
        endif
             
c-----------------------------------------------------------------------

      elseif ((nproc == 501) .or. (nproc == 502)
     &   .or. (nproc == 511) .or. (nproc == 512)
     &   .or. (nproc == 503) .or. (nproc == 513)
     &   .or. (nproc == 506) .or. (nproc == 516)) then
        kcase=kqq_ttw
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(11)='pp'
        ndim=20
        n2=0
        n3=1
        mass3=wmass
        width3=wwidth

        mcfmplotinfo= (/ 34, 78, 90, 345, 678, (0,j=1,45) /)
        
        if     (nproc == 501) then 
c-- 501 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(->b~(p6)+e^-(p7)+nu~(p8))+W^+(nu(p9),mu^+(p10))'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='nl'
          plabel(10)='ea'
          nqcdjets=2
        elseif   (nproc == 502) then 
c-- 502 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(->b~(p6)+e^-(p7)+nu~(p8))+W^+(nu(p9),mu^+(p10))[rid]'
          kcase=kttwldk
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='nl'
          plabel(10)='ea'
          nqcdjets=2
        elseif (nproc == 503) then 
c-- 503 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(->b~(p6)+q(p7)+q~(p8))+W^+(nu(p9),mu^+(p10))'
          nwz=+1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(7)='pp'
          plabel(8)='pp'
          plabel(9)='nl'
          plabel(10)='ea'
          nqcdjets=4
        elseif (nproc == 506) then 
c-- 506 '  f(p1)+f(p2) --> t(-->q(p3)+q~(p4)+b(p5))+t~(->b~(p6)+e^-(p7)+nu~(p8))+W^+(nu(p9),mu^+(p10))'
          nwz=+1
          plabel(3)='pp'
          plabel(4)='pp'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='nl'
          plabel(10)='ea'
          nqcdjets=4
        elseif (nproc == 511) then 
c-- 511 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(->b~(p6)+e^-(p7)+nu~(p8))+W^-(mu^-(p9),nu~(p10))'
          nwz=-1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='el'
          plabel(10)='na'
          nqcdjets=2
        elseif (nproc == 512) then 
c-- 512 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(->b~(p6)+e^-(p7)+nu~(p8))+W^-(mu^-(p9),nu~(p10))[rid]'
          kcase=kttwldk
          nwz=-1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='el'
          plabel(10)='na'
          nqcdjets=2
        elseif (nproc == 513) then 
c-- 513 '  f(p1)+f(p2) --> t(-->nu(p3)+e^+(p4)+b(p5))+t~(->b~(p6)+q(p7)+q~(p8))+W^-(mu^-(p9),nu~(p10))'
          nwz=-1
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(7)='pp'
          plabel(8)='pp'
          plabel(9)='el'
          plabel(10)='na'
          nqcdjets=4
        elseif (nproc == 516) then 
c-- 516 '  f(p1)+f(p2) --> t(-->q(p3)+q~(p4)+b(p5))+t~(->b~(p6)+e^-(p7)+nu~(p8))+W^-(mu^-(p9),nu~(p10))'
          nwz=-1
          plabel(3)='pp'
          plabel(4)='pp'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='el'
          plabel(10)='na'
          nqcdjets=4
        endif

        if ((nproc == 501) .or. (nproc == 502)
     & .or. (nproc == 511) .or. (nproc == 512)) then
        if (removebr) then
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=(brtop*brwen)**2*brwen
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            plabel(7)='ig'
            plabel(8)='ig'
            plabel(9)='ig'
            plabel(10)='ig'
            nqcdjets=0
        endif
        endif
        
c-----------------------------------------------------------------------

      elseif (nproc == 529) then
        kcase=kZbbmas
        swapxz=.true.
        call checkminzmass(1)
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='ig'
        plabel(6)='ig'
        q1=-1._dp
        l1=le
        r1=re
        ndim=10
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

        flav=6
        nflav=5

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
c--- total cross-section             
        if (removebr) then
          q1=zip
          plabel(3)='ig'
          plabel(4)='ig'
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=brzee
        endif

c-----------------------------------------------------------------------

      elseif ((nproc == 530) .or. (nproc == 531)) then
        kcase=kqq_ttz
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='el'
        plabel(8)='na'
        nwz=1

        ndim=20
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 78, 90, 345, 678, (0,j=1,45) /)
        
        if     (nproc == 530) then 
c--  530 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->nu~(p7)+e^-(p8)+b~(p6))+Z(e(p9),e~(p10))'
          plabel(9)='el'
          plabel(10)='ea'
          q1=-1._dp
          l1=le
          r1=re
        mb=zip
          if (removebr) then
            q1=zip
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=(brtop*brwen)**2*brzee
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            plabel(7)='ig'
            plabel(8)='ig'
            plabel(9)='ig'
            plabel(10)='ig'
          endif

        elseif (nproc == 531) then 
c--  531 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->nu~(p7)+e^-(p8)+b~(p6))+Z(b(p9),b~(p10))'
          plabel(9)='bq'
          plabel(10)='ba'
          q1=Q(5)*sqrt(xn)
          l1=l(5)*sqrt(xn)
          r1=r(5)*sqrt(xn)
        mb=zip
        endif

      elseif ((nproc == 532) .or. (nproc == 533)) then
        kcase=kqqtthz
        nwz=1

        ndim=20
        n2=0
        n3=1
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 78, 90, 345, 678, (0,j=1,45) /)
        
        nqcdjets=4

        if     (nproc == 532) then 
c--  532 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))+t~(-->q(p7)+q~(p8)+b~(p6))+Z(e(p9),e~(p10))'
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          plabel(7)='pp'
          plabel(8)='pp'
          plabel(9)='el'
          plabel(10)='ea'
          q1=-1._dp
          l1=le
          r1=re
          mb=zip
          if (removebr) then
            q1=zip
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2._dp*xn*(brtop*brwen)**2*brzee
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            plabel(7)='ig'
            plabel(8)='ig'
            plabel(9)='ig'
            plabel(10)='ig'
            nqcdjets=0
          endif
        elseif (nproc == 533) then 
c--  533 '  f(p1)+f(p2)-->t(-->q(p3)+q~(p4)+b(p5))+t~(-->nu~(p7)+e^-(p8)+b~(p6))+Z(e^-(p9),e^+(p10))'
           plabel(3)='pp'
           plabel(4)='pp'
           plabel(5)='bq'
           plabel(6)='ba'
           plabel(7)='na'
           plabel(8)='el'
           plabel(9)='el'
           plabel(10)='ea'
           q1=-1._dp
           l1=le
           r1=re
           mb=zip
           if (removebr) then
            q1=zip
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=2._dp*xn*(brtop*brwen)**2*brzee
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            plabel(7)='ig'
            plabel(8)='ig'
            plabel(9)='ig'
            plabel(10)='ig'
            nqcdjets=0
           endif
        endif

c-----------------------------------------------------------------------

        elseif ((nproc == 540) .or. (nproc == 541)) then
        kcase=kH_tjet
        mb=0
        swapxz=.true.
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=3
        ndim=10
        n2=0
        n3=0
        nflav=5
        plabel(5)='ig'
        plabel(6)='pp'
        hdecaymode='bqba'
        plabel(3)='bq'
        plabel(4)='ba'
        mass3=hmass
        width3=hwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 540) then
c-- 540        '  f(p1)+f(p2) --> H(b(p3)+b(p4))+t(p5)+q(p6)' 
          nwz=+1          
        elseif     (nproc == 541) then
c-- 541        '  f(p1)+f(p2) --> H(b(p3)+b(p4))+t~(p5)+q(p6)' 
          nwz=-1          
        endif

c--- total cross-section
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
        nqcdjets=1
          BrnRat=br
        endif
             
        elseif ((nproc == 544) .or. (nproc == 547)) then
        kcase=kH_tdkj
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        swapxz=.true.
        nqcdjets=3
        ndim=16
        hdecaymode='bqba'
        isub=1
        mb=0
        n2=0
        n3=0
        nflav=5
        mass3=hmass
        width3=hwidth

        mcfmplotinfo= (/ 34, 56, 567, (0,j=1,47) /)
        
        if     (nproc == 544) then
c--- 544   f(p1)+f(p2) --> H(b(p3)+b~(p4))+t(nu(p5)+e^+(p6)+b(p7))+q(p9)' 
        nwz=+1
        swapxz=.true.
        plabel(3)='bq'
        plabel(4)='ba'
        plabel(5)='nu'
        plabel(6)='ea'
        plabel(7)='bq'
        plabel(8)='pp'
        plabel(9)='pp'
        elseif (nproc == 547) then
c-- 547  f(p1)+f(p2) --> H(b(p3)+b~(p4))+t~(e-(p5)+nu(p6)+b~(p7))+q(p8)'
        nwz=-1
        swapxz=.true.
        plabel(3)='bq'
        plabel(4)='ba'
        plabel(5)='el'
        plabel(6)='nb'
        plabel(7)='ba'
        plabel(8)='pp'
        plabel(9)='pp'
        endif
c--- total cross-section             
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          nqcdjets=0
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          BrnRat=br*brtop*brwen
        endif


        elseif ((nproc == 550) .or. (nproc == 551)) then
        kcase=kH_tjet
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=1
        ndim=10
        swapxz=.true.
        mb=0
        n2=0
        n3=0
        nflav=5
        plabel(5)='ig'
        plabel(6)='pp'
        hdecaymode='gaga'
        plabel(3)='ga'
        plabel(4)='ga'
        mass3=hmass
        width3=hwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 550) then
c-- 551        '  f(p1)+f(p2) --> H(ga(p3)+ga(p4))+t(p5)+q(p6)' 
          nwz=+1          
        elseif     (nproc == 551) then
c-- 552        '  f(p1)+f(p2) --> H(ga(p3)+ga(p4))+t~(p5)+q(p6)' 
          nwz=-1          
        endif

c--- total cross-section
        if (removebr) then
          plabel(3)='ig'
          plabel(4)='ig'
          BrnRat=gamgambr
        endif
             
        elseif ((nproc == 554) .or. (nproc == 557)) then
        kcase=kH_tdkj
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=2
        ndim=16
        swapxz=.true.
        hdecaymode='gaga'
        isub=1
        mb=0
        n2=0
        n3=0
        nflav=5
        mass3=hmass
        width3=hwidth

        mcfmplotinfo= (/ 34, 56, 567, (0,j=1,47) /)
        
        if     (nproc == 554) then
c--- 554   f(p1)+f(p2) --> H(ga(p3)+ga(p4))+t(nu(p5)+e^+(p6)+b(p7))+q(p9)' 
        nwz=+1
        plabel(3)='ga'
        plabel(4)='ga'
        plabel(5)='nu'
        plabel(6)='ea'
        plabel(7)='bq'
        plabel(8)='pp'
        plabel(9)='pp'
        elseif (nproc == 557) then
c-- 557  f(p1)+f(p2) --> H(ga(p3)+ga(p4))+t~(e-(p5)+nu(p6)+b~(p7))+q(p8)'
        nwz=-1          
        plabel(3)='ga'
        plabel(4)='ga'
        plabel(5)='el'
        plabel(6)='nb'
        plabel(7)='ba'
        plabel(8)='pp'
        plabel(9)='pp'
        endif

c--- total cross-section
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          nqcdjets=0
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          BrnRat=gamgambr*brtop*brwen
        endif

c-----------------------------------------------------------------------

        elseif ((nproc == 560) .or. (nproc == 561)) then
        kcase=kZ_tjet
        call checkminzmass(1)
        swapxz=.true.
        nqcdjets=1
        ndim=10
        mb=0
        n2=0
        n3=0
        nflav=5
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='ig'
        plabel(6)='pp'
        plabel(7)='pp'
        q1=-1._dp
        l1=le
        r1=re
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 560) then
c-- 560        '  f(p1)+f(p2) --> Z(e-(p3)+e+(p4))+t(p5)+q(p6)' 
          nwz=+1
        elseif (nproc == 561) then
c-- 561        '  f(p1)+f(p2) --> Z(e-(p3)+ep(p4))+t~(p5)+q(p6)' 
          nwz=-1
        endif

c--- total cross-section
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          plabel(3)='ig'
          plabel(4)='ig'
          BrnRat=brzee
        endif

c--- 562 and 563
        elseif ((nproc == 562) .or. (nproc == 563)) then
        kcase=kZt2jet
        call checkminzmass(1)
        swapxz=.true.
        nqcdjets=2
        ndim=13
        mb=0
        n2=0
        n3=0
        nflav=5
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='ig'
        plabel(6)='pp'
        plabel(7)='pp'
        q1=-1._dp
        l1=le
        r1=re
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, (0,j=1,49) /)
        
        if     (nproc == 562) then
c-- 562        '  f(p1)+f(p2) --> Z(e-(p3)+e+(p4))+t(p5)+q(p6)+f(p7)' 
          nwz=+1          
        elseif (nproc == 563) then
c-- 563        '  f(p1)+f(p2) --> Z(e-(p3)+ep(p4))+t~(p5)+q(p6)+f(p7)' 
          nwz=-1          
        endif

c--- total cross-section
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          plabel(3)='ig'
          plabel(4)='ig'
          BrnRat=brzee
        endif

        elseif ((nproc == 564) .or. (nproc == 567)) then
        kcase=kZ_tdkj
        call checkminzmass(1)
        nqcdjets=2
        ndim=15
        swapxz=.true.
        mb=zip
        n2=0
        n3=0
        nflav=5
        q1=-1._dp
        l1=le
        r1=re
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 56, 567, (0,j=1,47) /)
        
        if     (nproc == 564) then
c--- 564 '  f(p1)+f(p2) --> Z(e-(p3)+e+(p4))+t(-->nu(p5)+e^+(p6)+b(p7))+q(p8)''N'
        nwz=+1          
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='nl'
        plabel(6)='ea'
        plabel(7)='bq'
        plabel(8)='pp'
        plabel(9)='pp'
        elseif (nproc == 567) then
c-- 567 '  f(p1)+f(p2) --> Z(e-(p3)+e+(p4))+t~(-->e-(p5)+nu(p6)+b~(p7))+q(p8)''N'
        nwz=-1          
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='ba'
        plabel(8)='pp'
        plabel(9)='pp'
        endif

c--- total cross-section             
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          nqcdjets=1
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          BrnRat=brzee*brtop*brwen
        endif

        elseif ((nproc == 566) .or. (nproc == 569)) then
        kcase=kZtdk2j
        call checkminzmass(1)
        nqcdjets=3
        ndim=18
        swapxz=.true.
        mb=zip
        n2=0
        n3=0
        nflav=5
        q1=-1._dp
        l1=le
        r1=re
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 56, 567, (0,j=1,47) /)
        
        if     (nproc == 566) then
c--- 566 '  f(p1)+f(p2) --> Z(e-(p3)+e+(p4))+t(-->nu(p5)+e^+(p6)+b(p7))+q(p8)+f(p9)'  'L'
        nwz=+1
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='nl'
        plabel(6)='ea'
        plabel(7)='bq'
        plabel(8)='pp'
        plabel(9)='pp'
        elseif (nproc == 569) then
c-- 569 '  f(p1)+f(p2) --> Z(e-(p3)+e+(p4))+t~(-->e-(p5)+nu(p6)+b~(p7))+q(p8)+f(p9)'  'L'
        nwz=-1          
        plabel(3)='el'
        plabel(4)='ea'
        plabel(5)='el'
        plabel(6)='na'
        plabel(7)='ba'
        plabel(8)='pp'
        plabel(9)='pp'
        endif

c--- total cross-section
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          nqcdjets=2
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          BrnRat=brzee*brtop*brwen
        endif

c-----------------------------------------------------------------------

        elseif ((nproc >= 601) .and. (nproc <= 602)) then
        kcase=kHHpair
        mb=0
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=0
        ndim=10
        n2=1
        n3=1
        nflav=5
        mass2=hmass
        width2=hwidth
        mass3=hmass
        width3=hwidth

        write(6,*) 'Higgs pair process only for zerowidth=T'
        zerowidth=.true.

        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)
        
        if     (nproc == 601) then
c-- 601 '  f(p1)+f(p2) --> H(b(p3)+b~(p4))+H(tau^-(p5)+tau^+(p6))' 'L'
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'
          hdecaymode2='tlta'
          plabel(5)='tl'
          plabel(6)='ta'
c--- total cross-section
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            BrnRat=2._dp*br*tautaubr      ! factor of 2 for identical particles

          endif
        elseif     (nproc == 602) then
c-- 602 '  f(p1)+f(p2) --> H(b(p3)+b~(p4))+H(gamma(p5)+gamma(p6))' 'L'
          hdecaymode='bqba'
          plabel(3)='bq'
          plabel(4)='ba'
          hdecaymode2='gaga'
          plabel(5)='ga'
          plabel(6)='ga'
c--- total cross-section
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            BrnRat=2._dp*br*gamgambr      ! factor of 2 for identical particles
          endif
        endif

c-----------------------------------------------------------------------

      elseif ((nproc >= 610) .and. (nproc <= 619)) then
        kcase=kWH1jet
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=1
        plabel(7)='pp'
        plabel(8)='pp'
        ndim=13
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=wmass
        width3=wwidth
        
        mcfmplotinfo= (/ 34, 56, 3456, (0,j=1,47) /)
        
        if (nproc == 610) then
          hdecaymode='tlta'
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='tl'
          plabel(6)='ta'
          nwz=+1
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen*tautaubr
          endif
        elseif (nproc == 611) then
          hdecaymode='bqba'
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          nqcdjets=3
          notag=2
          
          nwz=+1
c--- uncomment to sum over both W+ and W- 
c          nwz=2
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen*br
            nqcdjets=1
            notag=0
          endif

        elseif (nproc == 612) then
          hdecaymode='gaga'
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='ga'
          plabel(6)='ga'
          mb=zip
          nwz=+1
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen*gamgambr
          endif

c--- uncomment to sum over both W+ and W- 
c        nwz=2

        elseif (nproc == 613) then
          hdecaymode='wpwm'
          plabel(3)='nl'
          plabel(4)='ea'
          plabel(5)='nl'
          plabel(6)='ea'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          plabel(10)='pp'
          mcfmplotinfo= (/ 34, 56, 78, 5678, (0,j=1,46) /)
          ndim=19
          nwz=+1
          if (usescet .eqv. .false.) then
            new_pspace=.true.
          endif
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            plabel(7)='ig'
            plabel(8)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen**3*wwbr
          endif

c--- print warning if we're below threshold
          if (hmass < 2._dp*wmass) then
            write(6,*)
            write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
            write(6,*) 'may not yield sensible results - check the number'
            write(6,*) 'of integration points and the value of zerowidth'
            if (removebr) then
              write(6,*)
              write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
              stop
            endif
            if (zerowidth) then
              write(6,*) 'zerowidth=.true. and higgs decay below threshold'
              stop
            endif
          endif
        
        elseif (nproc == 615) then
          hdecaymode='tlta'
          plabel(3)='el'
          plabel(4)='na'
          plabel(5)='tl'
          plabel(6)='ta'
          nwz=-1
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen*tautaubr
          endif

        elseif (nproc == 616) then
          hdecaymode='bqba'
          plabel(3)='el'
          plabel(4)='na'
          plabel(5)='bq'
          plabel(6)='ba'
          nqcdjets=3
          notag=2
          nwz=-1
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen*br
            nqcdjets=1
            notag=0
          endif

        elseif (nproc == 617) then
          hdecaymode='gaga'
          plabel(3)='el'
          plabel(4)='na'
          plabel(5)='ga'
          plabel(6)='ga'
          nwz=-1
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen*gamgambr
          endif

        elseif (nproc == 618) then
          hdecaymode='wpwm'
          plabel(3)='el'
          plabel(4)='na'
          plabel(5)='nl'
          plabel(6)='ea'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          plabel(10)='pp'
          mcfmplotinfo= (/ 34, 56, 78, 5678, (0,j=1,46) /)
          ndim=19
          nwz=-1
          if (usescet .eqv. .false.) then
            new_pspace=.true.
          endif
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            plabel(7)='ig'
            plabel(8)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brwen**3*wwbr
          endif

c--- print warning if we're below threshold
          if (hmass < 2._dp*wmass) then
            write(6,*)
            write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
            write(6,*) 'may not yield sensible results - check the number'
            write(6,*) 'of integration points and the value of zerowidth'
            if (removebr) then
              write(6,*)
              write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
              stop
            endif
            if (zerowidth) then
              write(6,*) 'zerowidth=.true. and higgs decay below threshold'
              stop
            endif
         endif

        endif




c-----------------------------------------------------------------------

      elseif ((nproc >= 620) .and. (nproc <= 629)) then
        kcase=kZH1jet
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        nqcdjets=1
        plabel(7)='pp'
        plabel(8)='pp'
        ndim=13
        n2=1
        n3=1
        mass2=hmass
        width2=hwidth
        mass3=zmass
        width3=zwidth

        mcfmplotinfo= (/ 34, 56, (0,j=1,48) /)
        
        if (nproc == 620) then
          hdecaymode='tlta'
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='tl'
          plabel(6)='ta'
          q1=zip
          l1=le
          r1=re
          nwz=0
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee*tautaubr
          endif
        elseif (nproc == 621) then
          hdecaymode='bqba'
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='bq'
          plabel(6)='ba'
          nqcdjets=3
          notag=2
          
          q1=zip
          l1=le
          r1=re
          nwz=0
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee*br
            nqcdjets=1
            notag=0
          endif
        elseif (nproc == 622) then
          hdecaymode='gaga'
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='ga'
          plabel(6)='ga'
          q1=zip
          l1=le
          r1=re
          nwz=0
          if (removebr) then
            plabel(3)='ig'
            plabel(4)='ig'
            plabel(5)='ig'
            plabel(6)='ig'               
            call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
            BrnRat=brzee*gamgambr
          endif
       elseif (nproc == 623) then
          if (usescet .eqv. .false.) then
             new_pspace=.true.
          endif
          ndim=19
          hdecaymode='wpwm'
          plabel(3)='el'
          plabel(4)='ea'
          plabel(5)='nl'
          plabel(6)='ea'
          plabel(7)='el'
          plabel(8)='na'
          plabel(9)='pp'
          plabel(10)='pp'
       
          q1=zip
          l1=le
          r1=re
          nwz=0
       endif
c-----------------------------------------------------------------------

      elseif (nproc == 640) then
c--  640 '  f(p1)+f(p2)-->t(p3)+t~(p4)+H(p5)'
        kcase=ktottth
        swapxz=.true.
        plabel(3)='ig'
        plabel(4)='ig'
        plabel(5)='ig'
        plabel(6)='pp'
        nwz=1
        n2=0
        n3=0
        ndim=7
        mass2=mt
        
      elseif ((nproc >= 641) .and. (nproc <= 649)) then
c--  641 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))
c--         +t~(-->nu~(p7)+e^-(p8)+b~(p6))+H(b(p9)+b~(p10))'
c--      '  f(p1)+f(p2)-->t(p3+p4+p5)+t~(p6+p7+p8)+H(p9+p10)' (removebr=.true.)
        kcase=kqq_tth
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='el'
        plabel(8)='na'
        plabel(9)='bq'
        plabel(10)='ba'
        hdecaymode='bqba'
        nqcdjets=4

        nwz=1
        swapxz=.true.
        ndim=22
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth

        mcfmplotinfo= (/ 34, 78, 90, 345, 678, (0,j=1,45) /)
        
        if (nproc == 644) then
          plabel(7)='pp'
          plabel(8)='pp'
          nqcdjets=6
        elseif (nproc == 647) then
          plabel(3)='pp'
          plabel(4)='pp'
          nqcdjets=6
        endif         

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=(brtop*brwen)**2*br
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          plabel(8)='ig'
          plabel(9)='ig'
          plabel(10)='ig'
        nqcdjets=0
        endif

      elseif ((nproc >= 651) .and. (nproc <= 659)) then
c--  651 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))
c--         +t~(-->nu~(p7)+e^-(p8)+b~(p6))+H(ga(p9)+ga(p10))'
c--      '  f(p1)+f(p2)-->t(p3+p4+p5)+t~(p6+p7+p8)+H(p9+p10)' (removebr=.true.)
        kcase=kqq_tth
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='el'
        plabel(8)='na'
        plabel(9)='ga'
        plabel(10)='ga'
        hdecaymode='gaga'
        nqcdjets=2

        nwz=1
        ndim=22
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth

        mcfmplotinfo= (/ 34, 78, 90, 345, 678, (0,j=1,45) /)
        
      if (nproc == 654) then
        plabel(7)='pp'
        plabel(8)='pp'
      nqcdjets=4
      elseif (nproc == 657) then
        plabel(3)='pp'
        plabel(4)='pp'
      nqcdjets=4
      endif         

        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=(brtop*brwen)**2*gamgambr
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          plabel(8)='ig'
          plabel(9)='ig'
          plabel(10)='ig'
        nqcdjets=0
        endif

      elseif ((nproc >= 661) .and. (nproc <= 669)) then
c--  661 '  f(p1)+f(p2)-->t(-->nu(p3)+e^+(p4)+b(p5))
c--         +t~(-->nu~(p7)+e^-(p8)+b~(p6))+H)'
c--      '  f(p1)+f(p2)-->t(p3+p4+p5)+t~(p6+p7+p8)+H(WW)' (removebr=.true.)
        kcase=ktth_ww
        call sethparams(br,wwbr,zzbr,tautaubr,gamgambr,zgambr)
        plabel(3)='nl'
        plabel(4)='ea'
        plabel(5)='bq'
        plabel(6)='ba'
        plabel(7)='el'
        plabel(8)='na'
        plabel(9)='nl'
        plabel(10)='ea'
        plabel(11)='el'
        plabel(12)='na'
        hdecaymode='wpwm'
        nqcdjets=2

        nwz=1
        ndim=26
        n2=1
        n3=1
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth

        mcfmplotinfo= (/ 34, 78, 90, 345, 678, (0,j=1,45) /)
        
      if (nproc == 664) then
        plabel(7)='pp'
        plabel(8)='pp'
        nqcdjets=4
      elseif (nproc == 667) then
        plabel(3)='pp'
        plabel(4)='pp'
        nqcdjets=4
      endif         

c--- print warning if we're below threshold
        if (hmass < 2._dp*wmass) then
        write(6,*)
        write(6,*) 'WARNING: Higgs decay H->WW is below threshold and'
        write(6,*) 'may not yield sensible results - check the number'
        write(6,*) 'of integration points and the value of zerowidth'
      if (removebr) then
      write(6,*)
      write(6,*) 'Cannot remove H->WW BR, not defined below threshold'
        stop
      endif
      if (zerowidth) then
        write(6,*) 'zerowidth=.true. and higgs decay below threshold'
        stop
        endif
        endif
                 
        if (removebr) then
          call branch(brwen,brzee,brznn,brtau,brtop,brcharm)
          BrnRat=(brtop*brwen)**2*wwbr*brwen**2
          plabel(3)='ig'
          plabel(4)='ig'
          plabel(5)='ig'
          plabel(6)='ig'
          plabel(7)='ig'
          plabel(8)='ig'
          plabel(9)='ig'
          plabel(10)='ig'
          plabel(11)='ig'
          plabel(12)='ig'
        nqcdjets=0
        endif


      elseif((nproc>=800).and.(nproc<=805)) then
!     800 '  f(p1)+f(p2) --> V-->(X(p3)+X~(p4)) +f(p5) [Vector Mediator] ' 'N '
!     801 '  f(p1)+f(p2) --> A-->(X(p3)+X~(p4)) +f(p5) [Axial Vector Mediator] ' 'N '
!     802 '  f(p1)+f(p2) --> S-->(X(p3)+X~(p4)) +f(p5) [Scalar Mediator] ' 'N '
!     803 '  f(p1)+f(p2) --> PS-->(X(p3)+X~(p4)) +f(p5) [Pseudo Scalar Mediator] ' 'N '
!     804 '  f(p1)+f(p2) --> GG-->(X(p3)+X~(p4)) +f(p5) [Gluonic DM operator] ' 'N '
!     805 '  f(p1)+f(p2) --> S--(X(p3)+X~(p4)) +f(p5) [Scalar Mediator, mt loops] ' 'L'  
         kcase=kdm_jet
         plabel(3)='xm'
         plabel(4)='xa' 
         plabel(5)='pp'
         plabel(6)='pp'
         ndim=7
         nqcdjets=1
         if(nproc==800) then 
            dm_mediator='vector'
         elseif(nproc==801) then 
            dm_mediator='axvect'
         elseif(nproc==802) then 
            dm_mediator='scalar'
         elseif(nproc==803) then 
            dm_mediator='pseudo'
         elseif(nproc==804) then 
            dm_mediator='gluonO'
         elseif(nproc==805) then
            dm_mediator='scalmt'
         endif
         call read_dm_params()
!------- phase space setup
         n3=0
         mass3=medmass
         mass2=xmass
         width3=medwidth 

         
      elseif((nproc>=820).and.(nproc<=823)) then
!     820 '  f(p1)+f(p2) --> V-->(X(p3)+X~(p4)) +gamma(p5) [Vector Mediator] ' 'F '
!     821 '  f(p1)+f(p2) --> A-->(X(p3)+X~(p4)) +gamma(p5) [Axial Vector Mediator] ' 'F '
!     822 '  f(p1)+f(p2) --> S-->(X(p3)+X~(p4)) +gamma(p5) [Scalar Mediator] ' 'F '
!     823 '  f(p1)+f(p2) --> PS-->(X(p3)+X~(p4)) +gamma(p5) [Pseudo Scalar Mediator] ' 'F '
         kcase=kdm_gam
         plabel(3)='xm'
         plabel(4)='xa' 
         plabel(5)='ga'
         plabel(6)='pp'
         lastphot=5
         ndim=7
         nqcdjets=0
         if(nproc==820) then 
            dm_mediator='vector'
         elseif(nproc==821) then 
            dm_mediator='axvect'
         elseif(nproc==822) then
            dm_mediator='scalar'
         elseif(nproc==823) then
            dm_mediator='pseudo'
         endif
         call read_dm_params()
!------- phase space setup
         n3=0
         zerowidth=.true.
         mass3=medmass
         mass2=xmass

      elseif((nproc>=840).and.(nproc<=844)) then
!     840 '  f(p1)+f(p2) --> V-->(X(p3)+X~(p4)) +f(p5)+f(p6) [Vector Mediator] ' 'L '
!     841 '  f(p1)+f(p2) --> A-->(X(p3)+X~(p4)) +f(p5)+f(p6) [Axial Vector Mediator] ' 'L '
!     842 '  f(p1)+f(p2) --> S-->(X(p3)+X~(p4)) +f(p5)+f(p6) [Scalar Mediator] ' 'L '
!     843 '  f(p1)+f(p2) --> PS-->(X(p3)+X~(p4)) +f(p5)+f(p6) [Pseudo Scalar Mediator] ' 'L '
!     844 '  f(p1)+f(p2) --> GG-->(X(p3)+X~(p4)) +f(p5)+f(p6) [Gluonic DM operator] ' 'L '  
         kcase=kdm2jet
         plabel(3)='xm'
         plabel(4)='xa' 
         plabel(5)='pp'
         plabel(6)='pp'
         plabel(7)='pp'
         ndim=10
         nqcdjets=2
         if(nproc==840) then 
            dm_mediator='vector'
         elseif(nproc==841) then 
            dm_mediator='axvect'
         elseif(nproc==842) then 
            dm_mediator='scalar'
         elseif(nproc==843) then 
            dm_mediator='pseudo'
         elseif(nproc==844) then 
            dm_mediator='gluonO'
         endif
         call read_dm_params()
!------- phase space setup
         n2=0
         n3=0
         mass3=medmass
         mass2=xmass

      elseif((nproc>=845).and.(nproc<=848)) then
!     845 '  f(p1)+f(p2) --> V-->(X(p3)+X~(p4)) +gamma(p5)+f(p6) [Vector Mediator] ' 'L '
!     846 '  f(p1)+f(p2) --> A-->(X(p3)+X~(p4)) +gamma(p5)+f(p6) [Axial Vector Mediator] ' 'L '
!     847 '  f(p1)+f(p2) --> S-->(X(p3)+X~(p4)) +gamma(p5)+f(p6) [Scalar Mediator] ' 'L '
!     848 '  f(p1)+f(p2) --> PS-->(X(p3)+X~(p4)) +gamma(p5)+f(p6) [Pseudo Scalar Mediator] ' 'L '   
         kcase=kdm_gaj
         plabel(3)='xm'
         plabel(4)='xa' 
         plabel(5)='ga'
         plabel(6)='pp'
         plabel(7)='pp'
         lastphot=5
         ndim=10
         nqcdjets=1
         if(nproc==845) then 
            dm_mediator='vector'
         elseif(nproc==846) then 
            dm_mediator='axvect'
         elseif(nproc==847) then 
            dm_mediator='scalar'
         elseif(nproc==848) then 
            dm_mediator='pseudo'
         endif
         call read_dm_params()
!------- phase space setup
         n3=0
         mass3=medmass
         mass2=xmass
         

c-----------------------------------------------------------------------


      elseif (nproc/10 >= 90) then
        write(6,*) 'Setting part to lord and zerowidth to false'
        zerowidth=.false.
        kpart=klord
        if     (nproc == 902) then
          kcase=kvlchk2
          nwz=1
          ndim=4
          n3=1
          mass3=wmass
          width3=wwidth
        elseif (nproc == 903) then
          kcase=kvlchk3
          nwz=1
          ndim=7
          n3=1
          mass3=wmass
          width3=wwidth
        elseif (nproc == 904) then
          kcase=kvlchk4
          nwz=1
          ndim=10
          n2=1
          n3=1
          mass2=zmass
          width2=zwidth
          mass3=wmass
          width3=wwidth
        elseif (nproc == 905) then
          kcase=kvlchk5
          nwz=1
          ndim=13
          n2=1
          n3=1
          mass2=hmass
          width2=hwidth
          mass3=wmass
          width3=wwidth
        elseif (nproc == 906) then
          kcase=kvlchk6
          nwz=1
          ndim=16
          n2=1
          n3=1
          mass2=mt
          width2=twidth
          mass3=mt
          width3=twidth
        elseif (nproc == 908) then
          kcase=kvlchk8
          nwz=1
          ndim=22
          n2=1
          n3=1
          mass2=mt
          width2=twidth
          mass3=mt
          width3=twidth
        elseif (nproc == 909) then
          kcase=kvlchkm
          write(6,*) 'mb=',mb
          nwz=1
          ndim=10
          n2=1
          n3=1
          mass2=hmass
          width2=hwidth
          mass3=wmass
          width3=wwidth
        elseif (nproc == 910) then
          kcase=kvlchm3
          write(6,*) 'mt=',mt
          nwz=1
          ndim=7
          n2=0
          n3=0
          mass2=mt
          width2=twidth
          mass3=mt
          width3=twidth
        elseif (nproc == 911) then
          kcase=kvlchwt
          write(6,*) 'mt=',mt

          write(6,*) 'Setting zerowidth = .true.'
          zerowidth=.true.
             
          ndim=13
          mb=0
          n2=1
          n3=1
          mass2=mt
          width2=twidth
          mass3=wmass
          width3=wwidth
        elseif (nproc == 912) then
          kcase=kvlchwn
          write(6,*) 'mt=',mt

          write(6,*) 'Setting zerowidth = .true.'
          zerowidth=.true.
             
          ndim=7
          mb=0
          n2=0
          n3=1
          mass3=wmass
          width3=wwidth
        elseif (nproc == 913) then
          kcase=kvlchwg
          write(6,*) 'mt=',mt

          write(6,*) 'Setting zerowidth = .true.'
          zerowidth=.true.
             
          ndim=16
          mb=0
          n2=1
          n3=1
          mass2=mt
          width2=twidth
          mass3=wmass
          width3=wwidth             
          
        elseif (nproc == 914) then
          kcase=kvlchwh
          write(6,*) 'mt=',mt

          write(6,*) 'Setting zerowidth = .true.'
          zerowidth=.true.
             
          ndim=16
          mb=0
          n2=1
          n3=1
          mass2=mt
          width2=twidth
          mass3=wmass
          width3=wwidth
          
        endif
      else 
        call nprocinvalid()
      endif

c--- set notag (may be modified by user, with care!)
      call setnotag()

c--- set up alpha-s again (in case nflav was changed)
      call coupling2

c--- remove 2 dimensions from integration if decay is not included
      if (nodecay) ndim=ndim-2

c--- report on the removed BR, if necessary
      if (removebr) then
       if (rank == 0) then
        write(6,*)'****************************************************'
       endif
       if (BrnRat .ne. 1._dp) then
        if (rank == 0) then
        write(6,*)'*             Setting zerowidth to .true.          *'
        endif
        zerowidth=.true.
       endif
       if (rank == 0) then
        write(6,98) BrnRat
        write(6,*)'****************************************************'
       endif
      endif
      
c--- if needed for photon proceses reset to 1/137.0
!      if(reset_alphaEM) call reset_aem(aem)
      if(reset_alphaEM) call reset_aem(1._dp/137.0_dp)

c--- initialize arrays that are used in is_functions
      call init_is_functions()

c--- fill up CKM matrix
      call ckmfill(nwz)

c--- set flags to true unless we're doing W+2 jet or Z+2 jet
      if ( ((kcase.ne.kW_2jet) .and. (kcase.ne.kZ_2jet))
     & .or. (kpart==klord) ) then
        Qflag=.true.
        Gflag=.true.
      endif

      return

 43   write(6,*) 'problems opening process.DAT'
      stop

 44   write(6,*) 'Unimplemented process number, nproc = ',nproc, 
     & ' mcfm halted'
      stop
 
 98   format(' *             Brn.Rat. removed = ',  f11.7, '       *')
 99   format(' * ',a82,' *')
     
      end

      subroutine nprocinvalid()
      implicit none
      include 'types.f'
      
      integer:: nproc
      common/nproc/nproc

      write(6,*) 'chooser: Unimplemented case'
      write(6,*) 'nproc=',nproc      
      stop
      
      return 
      end
      
      subroutine checkminzmass(i)
      implicit none
      include 'types.f'
      include 'constants.f'
c--- Checks that the minimum invariant mass specified in the options
c--- file is not zero for boson 34 (i=1) or boson 56 (i=2)
      
      include 'limits.f'
      include 'zerowidth.f'
      integer:: i

c--- if generating exactly on-shell, there's nothing to worry about
      if (zerowidth) return
      
      if ((i == 1) .and. (wsqmin == zip)) then
        write(6,*)
        write(6,*) 'Please set m34min not equal to zero to'
        write(6,*) 'prevent the virtual photon from becoming real.'
        stop
      endif

      if ((i == 2) .and. (bbsqmin == zip)) then
        write(6,*)
        write(6,*) 'Please set m56min not equal to zero to'
        write(6,*) 'prevent the virtual photon from becoming real.'
        stop
      endif
      
      return
      end
      
      
