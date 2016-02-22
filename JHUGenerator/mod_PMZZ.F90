MODULE ModPMZZ
implicit none

contains


FUNCTION CalcMZZProbability(EHat,Ncalls)
use ModCrossSection
use ModKinematics
use ModMisc
use ModParameters
implicit none
real(8) :: DecayWeight,yRnd(1:22),Res,CalcMZZProbability
real(8) :: HiggsDK_Mom(1:4,1:13),Ehat
integer :: HiggsDK_IDUP(1:13),HiggsDK_ICOLUP(1:2,1:13)
integer :: evals,Ncalls

     CalcMZZProbability = 0d0
     do evals=1,Ncalls
         call random_number(yRnd)
         if( TauDecays.lt.0 ) then
             DecayWeight = EvalUnWeighted_DecayToVV(yRnd,.false.,EHat,Res,HiggsDK_Mom(1:4,6:9),HiggsDK_IDUP(1:9),HiggsDK_ICOLUP)
         else
             DecayWeight = EvalUnWeighted_DecayToTauTau(yRnd,.false.,EHat,Res,HiggsDK_Mom(1:4,4:13),HiggsDK_IDUP(1:13),HiggsDK_ICOLUP(1:2,1:13))
         endif
         CalcMZZProbability = CalcMZZProbability + DecayWeight
     enddo
     CalcMZZProbability = CalcMZZProbability/dble(evals)

RETURN
END FUNCTION




SUBROUTINE GetMZZdistribution(ScanMin, ScanMax, nsteps, doprint, MZZdistribution)
use ModCrossSection
use ModParameters
implicit none
real(8) :: DecayWeight,DecayWidth,DecayWidth0
real(8) :: Ehat
integer :: nscan, nsteps
integer,parameter :: Ncalls=200000
real(8) :: ScanMin, ScanMax
logical :: doprint
real(8) :: MZZdistribution(0:nsteps,1:2)


  if( PMZZ_mReso.lt.0 ) PMZZ_mReso = GetMZZProbability(M_Reso,Ncalls,-1d0,.false.)

  do nscan=0, nsteps

     EHat = ScanMin + (ScanMax-ScanMin)*nscan/dble(nsteps)
     DecayWidth = CalcMZZProbability(EHat,Ncalls)
     MZZdistribution(nscan,1:2) = (/EHat, DecayWidth/PMZZ_mReso/)

  enddo

RETURN
END SUBROUTINE



SUBROUTINE InitMZZdistribution(Ncalls)
!ScanRange: how far above and below M_Reso to go
!nmax:      number of points above and below M_Reso
use ModMisc
use ModParameters
implicit none
integer :: minindex, maxindex, n, mResoindex, Ncalls
real(8) :: minm4l, maxm4l

    if( PMZZminindex.ne.-1 .or. PMZZmaxindex.ne.-1 ) then
        call Error("InitMZZdistribution called twice!")
    endif

    mResoindex = PMZZsize/2
    PMZZ_mReso = GetMZZProbability(M_Reso,Ncalls,-1d0,.true.)
    PMZZdistribution(mResoindex,1:2) = (/M_Reso, PMZZ_mReso/)
    PMZZminindex = mResoindex
    PMZZmaxindex = mResoindex

    !extend out to Gamma in steps of 5 GeV
    call ExtendMZZdistribution(M_Reso+Ga_Reso,   5*GeV,  Ncalls)
    call ExtendMZZdistribution(M_Reso-Ga_Reso,   5*GeV,  Ncalls)
    !out to 3Gamma in steps of 10 GeV
    call ExtendMZZdistribution(M_Reso+3*Ga_Reso, 10*GeV, Ncalls)
    call ExtendMZZdistribution(M_Reso-3*Ga_Reso, 10*GeV, Ncalls)
    !then out to 5Gamma in steps of 20 GeV
    call ExtendMZZdistribution(M_Reso+5*Ga_Reso, 20*GeV, Ncalls)
    call ExtendMZZdistribution(M_Reso-5*Ga_Reso, 20*GeV, Ncalls)
    !leave it at that for now

RETURN
END SUBROUTINE


SUBROUTINE ExtendMZZdistribution(extendto, intervalsize, Ncalls)
use ModMisc
use ModParameters
implicit none
integer :: npoints, nsteps, Ncalls, minindex, maxindex
real(8) :: extendto, intervalsize, minm4l, maxm4l, ScanMin, ScanMax

    if( PMZZminindex.eq.-1 .or. PMZZmaxindex.eq.-1 ) call Error("Calling ExtendMZZdistribution without calling InitMZZdistribution first!")
    minm4l = PMZZdistribution(PMZZminindex,1)
    maxm4l = PMZZdistribution(PMZZmaxindex,1)
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

    call GetMZZDistribution(ScanMin, ScanMax, nsteps, .false., PMZZdistribution(minindex:maxindex,1:2))

RETURN
END SUBROUTINE



FUNCTION GetMZZProbability(EHat,Ncalls,intervalsize,usespline)
use ModCrossSection
use ModKinematics
use ModMisc
use ModParameters
implicit none
real(8) :: EHat, BigGamma, GetMZZProbability, PMZZ, intervalsize
integer :: Ncalls
logical :: usespline

     GetMZZProbability = 1d0

     if( ReweightDecay ) then

         if( usespline ) then
             if( intervalsize.gt.0d0 ) then
                 call ExtendMZZdistribution(EHat, intervalsize, Ncalls)
             else
                 call ExtendMZZdistribution(EHat, 40*GeV, Ncalls)
             endif
             call EvaluateSpline(EHat, PMZZdistribution(PMZZminindex:PMZZmaxindex,1:2), PMZZmaxindex-PMZZminindex+1, PMZZ)
             GetMZZProbability = GetMZZProbability * PMZZ
         else
             GetMZZProbability = CalcMZZProbability(EHat,Ncalls)
         endif

         if( WidthSchemeIn.eq.3 ) then
             !need to divide by overcompensated factor of m4l*BigGamma from POWHEG
             call HTO_gridHt(EHat/GeV,BigGamma)
             GetMZZProbability = GetMZZProbability / (EHat*BigGamma)
         elseif( WidthSchemeIn.eq.1 ) then
             !divide by overcompensated m4l*gamma_running in the numerator
             GetMZZProbability = GetMZZProbability / (EHat**2*Ga_Reso/M_Reso)
         elseif( WidthSchemeIn.eq.2 ) then
             !numerator is a constant M_Reso*Ga_Reso, do nothing
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



SUBROUTINE PrintMZZdistribution()
use ModCrossSection
use ModParameters
implicit none
real(8) :: DecayWeight,DecayWidth,DecayWidth0
real(8) :: Ehat
integer :: nscan
integer,parameter :: nmax=20, Ncalls=200000
real(8),parameter :: ScanRange=120d0*GeV


  DecayWidth0 = GetMZZProbability(M_Reso,Ncalls,-1d0,.false.)

  do nscan=-nmax,+nmax,1

     EHat = M_Reso+ScanRange*nscan/dble(nmax)
     DecayWidth = GetMZZProbability(EHat,Ncalls,-1d0,.false.)
     write(*,"(1F10.5,1PE16.9)") EHat*100d0,DecayWidth/DecayWidth0

  enddo

RETURN

END SUBROUTINE

END MODULE
