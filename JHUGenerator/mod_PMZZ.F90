MODULE ModPMZZ
implicit none

contains


FUNCTION CalcMZZProbability(EHat)
use ModCrossSection
use ModKinematics
use ModMisc
use ModParameters
implicit none
real(8) :: DecayWeight,yRnd(1:22),Res,CalcMZZProbability
real(8) :: HiggsDK_Mom(1:4,1:13),Ehat
integer :: HiggsDK_IDUP(1:13),HiggsDK_ICOLUP(1:2,1:13)
integer :: evals

     CalcMZZProbability = 0d0
     do evals=1,PMZZEvals
         call random_number(yRnd)
         if( TauDecays.lt.0 ) then
             DecayWeight = EvalUnWeighted_DecayToVV(yRnd,.false.,EHat,Res,HiggsDK_Mom(1:4,6:9),HiggsDK_IDUP(1:9),HiggsDK_ICOLUP)
         else
             DecayWeight = EvalUnWeighted_DecayToTauTau(yRnd,.false.,EHat,Res,HiggsDK_Mom(1:4,4:13),HiggsDK_IDUP(1:13),HiggsDK_ICOLUP(1:2,1:13))
         endif
         CalcMZZProbability = CalcMZZProbability + DecayWeight
     enddo
     CalcMZZProbability = CalcMZZProbability/dble(PMZZEvals)
     !print *, EHat*100d0, CalcMZZProbability

RETURN
END FUNCTION




SUBROUTINE GetMZZDistribution(ScanMin, ScanMax, nsteps, MZZdistribution)
use ModCrossSection
use ModParameters
implicit none
real(8) :: DecayWeight,DecayWidth,DecayWidth0
real(8) :: Ehat
integer :: nscan, nsteps
real(8) :: ScanMin, ScanMax
real(8) :: MZZdistribution(0:nsteps,1:2)


  do nscan=0, nsteps

     if( nsteps.eq.0 ) then
         EHat = (ScanMin+ScanMax)/2
     else
         EHat = ScanMin + (ScanMax-ScanMin)*nscan/dble(nsteps)
     endif
     DecayWidth = CalcMZZProbability(EHat)
     MZZdistribution(nscan,1:2) = (/EHat, DecayWidth/)

  enddo

RETURN
END SUBROUTINE



SUBROUTINE InitMZZdistribution()
!ScanRange: how far above and below M_Reso to go
!nmax:      number of points above and below M_Reso
use ModMisc
use ModParameters
implicit none
integer :: minindex, maxindex, n, mResoindex
real(8) :: minm4l, maxm4l

    if( ReadPMZZ ) then
        call ReadMZZdistribution(PMZZfile)
    else
        if( PMZZminindex.ne.-1 .or. PMZZmaxindex.ne.-1 ) then
            call Error("InitMZZdistribution called twice!")
        endif

        mResoindex = PMZZsize/2
        PMZZ_mReso = CalcMZZProbability(M_Reso)
        PMZZdistribution(mResoindex,1:2) = (/M_Reso, PMZZ_mReso/)
        PMZZminindex = mResoindex
        PMZZmaxindex = mResoindex
    endif

    !extend out to Gamma in steps of 1 GeV
    call ExtendMZZdistribution(M_Reso+Ga_Reso,   1*GeV)
    call ExtendMZZdistribution(M_Reso-Ga_Reso,   1*GeV)
    !then out to the edge of the distribution in intervals of 5 GeV
    !  with 3 extra points to avoid boundary effects
    call ExtendMZZdistribution(maxinputmHstar+15*GeV, 5*GeV)
    call ExtendMZZdistribution(mininputmHstar-15*GeV, 5*GeV)
    !make sure there are a few data points even for a tiny width
    !so that we don't get NaN
    call ExtendMZZdistribution(M_Reso+2*GeV,     1*GeV)
    call ExtendMZZdistribution(M_Reso-2*GeV,     1*GeV)
    !leave it at that for now

    call WriteMZZdistribution(PMZZfile)

RETURN
END SUBROUTINE


SUBROUTINE ExtendMZZdistribution(extendto, intervalsize)
use ModMisc
use ModParameters
implicit none
integer :: npoints, nsteps, minindex, maxindex
real(8) :: extendto, intervalsize, minm4l, maxm4l, ScanMin, ScanMax

    if( PMZZminindex.eq.-1 .or. PMZZmaxindex.eq.-1 ) call Error("Calling ExtendMZZdistribution without calling InitMZZdistribution first!")
    minm4l = PMZZdistribution(PMZZminindex,1)
    maxm4l = PMZZdistribution(PMZZmaxindex,1)
    extendto = max(extendto, mininputmHstar-3*intervalsize)
    extendto = min(extendto, maxinputmHstar+3*intervalsize)
    if( extendto.gt.maxm4l ) then
        ScanMax = extendto
        npoints = ceiling((extendto - maxm4l)/intervalsize)
        nsteps = npoints-1
        ScanMin = maxm4l + (extendto - maxm4l)/npoints
        minindex = PMZZmaxindex+1
        maxindex = PMZZmaxindex+npoints
        PMZZmaxindex = maxindex
    elseif( extendto.lt.minm4l ) then
        ScanMin = extendto
        npoints = ceiling((minm4l - extendto)/intervalsize)
        nsteps = npoints-1
        ScanMax = minm4l - (minm4l - extendto)/npoints
        maxindex = PMZZminindex-1
        minindex = PMZZminindex-npoints
        PMZZminindex = minindex
    else
        return
    endif

    if( maxindex.gt.PMZZSize.or.minindex.lt.1 ) then
        call Error("Not enough room for the P_decay(m4f) distribution!  Try increasing PMZZsize in mod_Parameters")
    endif
    call GetMZZDistribution(ScanMin, ScanMax, nsteps, PMZZdistribution(minindex:maxindex,1:2))

RETURN
END SUBROUTINE



FUNCTION GetMZZProbability(EHat,intervalsize,usespline)
use ModCrossSection
use ModKinematics
use ModMisc
use ModParameters
implicit none
real(8) :: EHat, BigGamma, GetMZZProbability, PMZZ, intervalsize
logical :: usespline

     GetMZZProbability = 1d0

     if( ReweightDecay ) then

         if( usespline ) then
             if( intervalsize.gt.0d0 ) then
                 !shouldn't need this anymore, with maxinputmHstar and mininputmHstar above
                 call ExtendMZZdistribution(EHat+3*intervalsize, intervalsize)
                 call ExtendMZZdistribution(EHat-3*intervalsize, intervalsize)
             else
                 !shouldn't need this anymore, with maxinputmHstar and mininputmHstar above
                 call ExtendMZZdistribution(EHat+15*GeV,         5*GeV)
                 call ExtendMZZdistribution(EHat-15*GeV,         5*GeV)
             endif
             call EvaluateSpline(EHat, PMZZdistribution(PMZZminindex:PMZZmaxindex,1:2), PMZZmaxindex-PMZZminindex+1, PMZZ)
             GetMZZProbability = GetMZZProbability * PMZZ
         else
             GetMZZProbability = CalcMZZProbability(EHat)
         endif

         if( WidthSchemeIn.eq.3 ) then
             !need to divide by overcompensated factor of m4l*BigGamma from POWHEG
             call HTO_gridHt(EHat/GeV,BigGamma)
             GetMZZProbability = GetMZZProbability / (EHat*BigGamma)
         elseif( WidthSchemeIn.eq.1 ) then
             !divide by overcompensated m4l*gamma_running in the numerator
             GetMZZProbability = GetMZZProbability / (EHat**2*Ga_Reso/M_Reso)
         elseif( WidthSchemeIn.eq.2 ) then
             continue ! numerator is a constant M_Reso*Ga_Reso, do nothing
         else
             call Error("Invalid WidthSchemeIn!")
         endif
     endif

     if( WidthScheme.ne.WidthSchemeIn ) then
         !reweight the propagator
         GetMZZProbability = GetMZZProbability * GetBWPropagator(EHat**2, WidthScheme) / GetBWPropagator(EHat**2, WidthSchemeIn)
     endif


RETURN
END FUNCTION


SUBROUTINE WriteMZZdistribution(outfile)
use ModParameters
implicit none
integer :: i
character(len=*) :: outfile

    open(unit=io_TmpFile,file=outfile,form='unformatted',status='replace')
    do i=PMZZminindex, PMZZmaxindex
        write(io_TmpFile) i, PMZZdistribution(i,1), PMZZdistribution(i,2)
    end do
    close(io_TmpFile)

RETURN
END SUBROUTINE

SUBROUTINE ReadMZZdistribution(infile)
use ModParameters
implicit none
integer :: i, stat
real(8) :: EHat, PMZZ
character(len=*) :: infile

    open(unit=io_TmpFile,file=trim(infile),form='unformatted',access= 'sequential',status='old')
    do while (.true.)
        read(io_TmpFile,iostat=stat,END=99) i, EHat, PMZZ
        PMZZdistribution(i,1:2) = (/EHat, PMZZ/)
        if (PMZZminindex.le.0 .or. PMZZminindex.gt.i) PMZZminindex = i
        if (PMZZmaxindex.lt.i) PMZZmaxindex = i
    end do
99  CONTINUE
    close(io_TmpFile)

RETURN
END SUBROUTINE

SUBROUTINE PrintMZZdistribution()
use ModCrossSection
use ModParameters
implicit none
real(8) :: DecayWidth,DecayWidth0
real(8) :: Ehat,minEhat, maxEhat
integer :: nscan
real(8),parameter :: ScanRange=120d0*GeV


  print *, " finding P_decay(m4l) distribution with ", PMZZEvals, " calls per point"
  DecayWidth0 = GetMZZProbability(M_Reso,-1d0,.false.)
  minEhat = dreal(PrintPMZZ)
  maxEhat = dimag(PrintPMZZ)

  do nscan=0,PrintPMZZIntervals

     EHat = minEhat+nscan*(maxEhat-minEhat)/PrintPMZZIntervals
     if( PrintPMZZIntervals.eq.0 ) EHat = minEhat
     DecayWidth = GetMZZProbability(EHat,-1d0,.false.)
     write(*,"(1F10.5,1PE16.9)") EHat*100d0,DecayWidth/DecayWidth0

  enddo

RETURN

END SUBROUTINE

END MODULE
