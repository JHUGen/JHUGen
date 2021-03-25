MODULE ModTTBHiggs
use ModParameters
implicit none


public :: EvalXSec_PP_TTBH, EvalXSec_PP_BBBH, EvalAmp_GG_TTBH, EvalAmp_QQB_TTBH
public :: InitProcess_TTBH
public :: ExitProcess_TTBH
private

integer,parameter :: ColorlessTag = 1


type :: Particle
   integer  :: PartType
   integer  :: ExtRef
   integer  :: Helicity
   real(8)  :: Mass
   real(8)  :: Mass2
   complex(8) :: Mom(1:4)
   complex(8) :: Pol(1:4)
end type

type :: PtrToParticle
   integer,pointer  :: PartType
   integer,pointer  :: ExtRef
   real(8),pointer  :: Mass
   real(8),pointer  :: Mass2
   integer ,pointer :: Helicity
   complex(8),pointer :: Mom(:)
   complex(8),pointer :: Pol(:)
end type



type :: TreeProcess
   integer :: NumPart
   integer :: NumQua
   integer :: NumSca
   integer :: NumW
   integer :: NumV
   integer :: BosonVertex
   integer,allocatable :: NumGlu(:)
   integer,allocatable :: PartRef(:)
   integer,allocatable :: PartType(:)
   type(PtrToParticle),allocatable :: Quarks(:)
   type(PtrToParticle),allocatable :: Gluons(:)
   type(PtrToParticle) :: Boson
   type(PtrToParticle),allocatable :: Scalars(:)
end type


INTERFACE OPERATOR (.Ndot.)
   module procedure FourVecDot
END INTERFACE OPERATOR (.Ndot.)


real(8), parameter :: PropCut = 1.0d-8
integer, parameter :: Dv=4,Ds=4

type(Particle),save    :: ExtParticles(1:7)
type(TreeProcess),save :: TheTreeAmps_GG_TTBH(1:2)
type(TreeProcess),save :: TheTreeAmps_QQB_TTBH(1:1)


complex(8) :: couplHTT_right_dyn,couplHTT_left_dyn




   CONTAINS


SUBROUTINE EvalXSec_PP_TTBH(Mom,SelectProcess,Res)
implicit none
real(8), intent(in) :: Mom(1:4,1:13)
integer,intent(in) :: SelectProcess! 0=gg, 1=qqb, 2=all
real(8), intent(out) :: Res(-5:5,-5:5)
real(8) :: MatElSq_GG,MatElSq_QQB,MatElSq_QBQ
integer :: iq

   call InitProcess_TTBH()

   MatElSq_QQB = 0d0
   MatElSq_QBQ = 0d0
   MatElSq_GG  = 0d0
   if( SelectProcess.eq.0 ) then
      call EvalAmp_GG_TTBH(Mom(1:4,1:13),MatElSq_GG)
   elseif( SelectProcess.eq.1 ) then
      call EvalAmp_QQB_TTBH(Mom(1:4,1:13),MatElSq_QQB)
      MatElSq_QBQ = MatElSq_QQB
   else
      call EvalAmp_GG_TTBH(Mom(1:4,1:13),MatElSq_GG)
      call EvalAmp_QQB_TTBH(Mom(1:4,1:13),MatElSq_QQB)
      MatElSq_QBQ = MatElSq_QQB
   endif
   do iq=0,5
      if(iq.eq.pdfGlu_) then
         Res(iq,iq) = MatElSq_GG
      else
         Res(iq,-iq) = MatElSq_QQB
         Res(-iq,iq) = MatElSq_QBQ
      endif
   enddo

   call ExitProcess_TTBH()
   RETURN
END SUBROUTINE

SUBROUTINE EvalXSec_PP_BBBH(Mom,SelectProcess,Res)
implicit none
real(8), intent(in) :: Mom(1:4,1:13)
integer,intent(in) :: SelectProcess! 0=gg, 1=qqb, 2=all
real(8), intent(out) :: Res(-5:5,-5:5)
real(8) :: tmptopmass
integer :: tmptopdec
   tmptopmass = M_Top
   tmptopdec = TopDecays

   M_Top = M_Bot
   TopDecays=0
   call EvalXSec_PP_TTBH(Mom,SelectProcess,Res)

   M_Top = tmptopmass
   TopDecays = tmptopdec
   RETURN
END SUBROUTINE



SUBROUTINE InitProcess_TTBH()
implicit none
integer :: NumQuarks,NumGluons,NumBoson
integer :: NumTrees
integer :: iTree,NumParticles


! gg->ttbar+H
  NumQuarks=2; NumGluons=2; NumBoson=1;
  NumParticles=NumQuarks+NumGluons+NumBoson
  NumTrees=2
  do iTree=1,NumTrees
      TheTreeAmps_GG_TTBH(iTree)%NumPart=NumParticles
      TheTreeAmps_GG_TTBH(iTree)%NumQua=NumQuarks
      allocate( TheTreeAmps_GG_TTBH(iTree)%NumGlu(0:NumQuarks+1) )
      allocate( TheTreeAmps_GG_TTBH(iTree)%PartRef(1:NumParticles) )
      allocate( TheTreeAmps_GG_TTBH(iTree)%PartType(1:NumParticles) )
      allocate( TheTreeAmps_GG_TTBH(iTree)%Quarks(1:NumQuarks) )
      allocate( TheTreeAmps_GG_TTBH(iTree)%Gluons(1:NumGluons) )
      TheTreeAmps_GG_TTBH(iTree)%NumGlu(0) = NumGluons
      TheTreeAmps_GG_TTBH(iTree)%NumSca = 0
  enddo

! qqbar->ttbar+H
  NumQuarks=4; NumGluons=0; NumBoson=1;
  NumParticles=NumQuarks+NumGluons+NumBoson
  NumTrees=1
  do iTree=1,NumTrees
      TheTreeAmps_QQB_TTBH(iTree)%NumPart=NumParticles
      TheTreeAmps_QQB_TTBH(iTree)%NumQua=NumQuarks
      allocate( TheTreeAmps_QQB_TTBH(iTree)%NumGlu(0:NumQuarks+1) )
      allocate( TheTreeAmps_QQB_TTBH(iTree)%PartRef(1:NumParticles) )
      allocate( TheTreeAmps_QQB_TTBH(iTree)%PartType(1:NumParticles) )
      allocate( TheTreeAmps_QQB_TTBH(iTree)%Quarks(1:NumQuarks) )
      allocate( TheTreeAmps_QQB_TTBH(iTree)%Gluons(1:NumGluons) )
      TheTreeAmps_QQB_TTBH(iTree)%NumGlu(0) = NumGluons
      TheTreeAmps_QQB_TTBH(iTree)%NumSca = 0
  enddo

  ExtParticles(1)%PartType = ATop_
  ExtParticles(1)%ExtRef   = 1
  ExtParticles(1)%Mass = m_Top
  ExtParticles(1)%Mass2= ExtParticles(1)%Mass**2
  ExtParticles(1)%Helicity = 0

  ExtParticles(2)%PartType = Top_
  ExtParticles(2)%ExtRef   = 2
  ExtParticles(2)%Mass = m_Top
  ExtParticles(2)%Mass2= ExtParticles(2)%Mass**2
  ExtParticles(2)%Helicity = 0

  ExtParticles(3)%PartType = Glu_
  ExtParticles(3)%ExtRef   = 3
  ExtParticles(3)%Mass = 0d0
  ExtParticles(3)%Mass2= 0d0
  ExtParticles(3)%Helicity = 0

  ExtParticles(4)%PartType = Glu_
  ExtParticles(4)%ExtRef   = 4
  ExtParticles(4)%Mass = 0d0
  ExtParticles(4)%Mass2= 0d0
  ExtParticles(4)%Helicity = 0

  ExtParticles(5)%PartType = AStr_
  ExtParticles(5)%ExtRef   = 5
  ExtParticles(5)%Mass = 0d0
  ExtParticles(5)%Mass2= 0d0
  ExtParticles(5)%Helicity = 0

  ExtParticles(6)%PartType = Str_
  ExtParticles(6)%ExtRef   = 6
  ExtParticles(6)%Mass = 0d0
  ExtParticles(6)%Mass2= 0d0
  ExtParticles(6)%Helicity = 0

  ExtParticles(7)%PartType = Hig_
  ExtParticles(7)%ExtRef   = 7
  ExtParticles(7)%Mass = m_Reso
  ExtParticles(7)%Mass2= ExtParticles(5)%Mass**2
  ExtParticles(7)%Helicity = 0


  TheTreeAmps_GG_TTBH(1)%PartRef(1:5) = (/3,4,1,7,2/)! (/1,7,2,3,4/)
  TheTreeAmps_GG_TTBH(2)%PartRef(1:5) = (/3,1,7,2,4/)! (/1,7,2,4,3/)
  do iTree=1,2
      call LinkTreeParticles(TheTreeAmps_GG_TTBH(iTree),ExtParticles(1:7))
  enddo

  TheTreeAmps_QQB_TTBH(1)%PartRef(1:5) = (/5,6,1,7,2/)! (/1,7,2,5,6/)
  call LinkTreeParticles(TheTreeAmps_QQB_TTBH(1),ExtParticles(1:7))


  couplHTT_right_dyn = m_top/vev/2d0 * ( kappa + (0d0,1d0)*kappa_tilde )
  couplHTT_left_dyn  = m_top/vev/2d0 * ( kappa - (0d0,1d0)*kappa_tilde )


RETURN
END SUBROUTINE


SUBROUTINE ExitProcess_TTBH()
implicit none
integer :: NumTrees, iTree


! gg->ttbar+H
  NumTrees=2
  do iTree=1,NumTrees
      !call UnLinkTreeParticles(TheTreeAmps_GG_TTBH(iTree))

      if(allocated( TheTreeAmps_GG_TTBH(iTree)%NumGlu )) then
         deallocate( TheTreeAmps_GG_TTBH(iTree)%NumGlu )
      endif
      if(allocated( TheTreeAmps_GG_TTBH(iTree)%PartRef )) then
         deallocate( TheTreeAmps_GG_TTBH(iTree)%PartRef )
      endif
      if(allocated( TheTreeAmps_GG_TTBH(iTree)%PartType )) then
         deallocate( TheTreeAmps_GG_TTBH(iTree)%PartType )
      endif
      if(allocated( TheTreeAmps_GG_TTBH(iTree)%Quarks )) then
         deallocate( TheTreeAmps_GG_TTBH(iTree)%Quarks )
      endif
      if(allocated( TheTreeAmps_GG_TTBH(iTree)%Gluons )) then
         deallocate( TheTreeAmps_GG_TTBH(iTree)%Gluons )
      endif
      if(allocated( TheTreeAmps_GG_TTBH(iTree)%Scalars )) then
         deallocate( TheTreeAmps_GG_TTBH(iTree)%Scalars )
      endif
  enddo

  NumTrees=1
  do iTree=1,NumTrees
      !call UnLinkTreeParticles(TheTreeAmps_QQB_TTBH(iTree))

      if(allocated( TheTreeAmps_QQB_TTBH(iTree)%NumGlu )) then
         deallocate( TheTreeAmps_QQB_TTBH(iTree)%NumGlu )
      endif
      if(allocated( TheTreeAmps_QQB_TTBH(iTree)%PartRef )) then
         deallocate( TheTreeAmps_QQB_TTBH(iTree)%PartRef )
      endif
      if(allocated( TheTreeAmps_QQB_TTBH(iTree)%PartType )) then
         deallocate( TheTreeAmps_QQB_TTBH(iTree)%PartType )
      endif
      if(allocated( TheTreeAmps_QQB_TTBH(iTree)%Quarks )) then
         deallocate( TheTreeAmps_QQB_TTBH(iTree)%Quarks )
      endif
      if(allocated( TheTreeAmps_QQB_TTBH(iTree)%Gluons )) then
         deallocate( TheTreeAmps_QQB_TTBH(iTree)%Gluons )
      endif
      if(allocated( TheTreeAmps_QQB_TTBH(iTree)%Scalars )) then
         deallocate( TheTreeAmps_QQB_TTBH(iTree)%Scalars )
      endif
  enddo


RETURN
END SUBROUTINE



SUBROUTINE EvalAmp_GG_TTBH(Mom,SqAmp)
use ModTopDecay
implicit none
real(8) :: Mom(1:4,1:13),SqAmp
complex(8) :: ResOffSh(1:4,1:2),Res(1:2,1:2)
complex(8) :: GluPol(1:4,1:2,1:2)
integer :: hel4,TopHel1,TopHel2,nhel
real(8),parameter :: c_aa=64.D0/3.D0, c_ab=-8.D0/3.D0
integer, parameter :: inLeft=1,inRight=2,Hbos=3,tbar=4,t=5,  bbar=6,Wm=7,lepM=8,nubar=9,  b=10,Wp=11,lepP=12,nu=13
SqAmp = 0d0


     ExtParticles(1)%Mom(1:4) = Mom(1:4,tbar)
     ExtParticles(2)%Mom(1:4) = Mom(1:4,t)
     ExtParticles(3)%Mom(1:4) =-Mom(1:4,inLeft)
     ExtParticles(4)%Mom(1:4) =-Mom(1:4,inRight)
     ExtParticles(7)%Mom(1:4) = Mom(1:4,Hbos)

     if( TOPDECAYS.ne.0 ) then
        call TopDecay(ATop_,(/Mom(1:4,bbar),Mom(1:4,lepM),Mom(1:4,nubar)/),ExtParticles(1)%Pol(1:4))
        call TopDecay(Top_,(/Mom(1:4,b),Mom(1:4,lepP),Mom(1:4,nu)/),ExtParticles(2)%Pol(1:4))
     endif
     ExtParticles(7)%Pol(1:4) = 1d0
!    call HDecay(ExtParticles(7),DK_LO,MomExt(1:4,12:13))
     GluPol(1:4,1,1) = pol_mless(ExtParticles(3)%Mom(1:4),+1,outgoing=.true.)
     GluPol(1:4,1,2) = pol_mless(ExtParticles(3)%Mom(1:4),-1,outgoing=.true.)
     GluPol(1:4,2,1) = pol_mless(ExtParticles(4)%Mom(1:4),+1,outgoing=.true.)
     GluPol(1:4,2,2) = pol_mless(ExtParticles(4)%Mom(1:4),-1,outgoing=.true.)
!      GluPol(1:4,1,1) = ExtParticles(3)%Mom(1:4);  GluPol(1:4,1,2) = ExtParticles(3)%Mom(1:4); print *, "checking gauge invariance"


     nhel=-1
     if( TOPDECAYS.EQ.0 ) nhel=+1
     do TopHel1=-1,nhel,2
     do TopHel2=-1,nhel,2
     if( TOPDECAYS.eq.0 ) then
             call ubarSpi_Dirac(ExtParticles(2)%Mom(1:4),M_Top,TopHel1,ExtParticles(2)%Pol(1:4))
             call    vSpi_Dirac(ExtParticles(1)%Mom(1:4),M_Top,TopHel2,ExtParticles(1)%Pol(1:4))
     endif
     do hel4=1,2

        ExtParticles(4)%Pol(1:4) = GluPol(1:4,2,hel4)
!         ExtParticles(4)%Pol(1:4) = ExtParticles(4)%Mom(1:4); print *, "checking gauge invariance"
        call new_calc_ampl(0,0,TheTreeAmps_GG_TTBH(1),ResOffSh(1:4,1))
        call new_calc_ampl(0,0,TheTreeAmps_GG_TTBH(2),ResOffSh(1:4,2))

        Res(1,1) = (ResOffSh(1:4,1).Ndot.GluPol(1:4,1,1))! col1 hel+
        Res(2,1) = (ResOffSh(1:4,2).Ndot.GluPol(1:4,1,1))! col2 hel+
        Res(1,2) = (ResOffSh(1:4,1).Ndot.GluPol(1:4,1,2))! col1 hel-
        Res(2,2) = (ResOffSh(1:4,2).Ndot.GluPol(1:4,1,2))! col2 hel-

        SqAmp = SqAmp   &
              + c_aa * dreal( Res(1,1)*dconjg(Res(1,1)) + Res(1,2)*dconjg(Res(1,2)) )  &
              + c_ab * dreal( Res(1,1)*dconjg(Res(2,1)) + Res(1,2)*dconjg(Res(2,2)) )  &
              + c_ab * dreal( Res(2,1)*dconjg(Res(1,1)) + Res(2,2)*dconjg(Res(1,2)) )  &
              + c_aa * dreal( Res(2,1)*dconjg(Res(2,1)) + Res(2,2)*dconjg(Res(2,2)) )
    enddo
    enddo
    enddo

    SqAmp = SqAmp * SpinAvg * GluonColAvg**2 * (4d0*Pi*alphas)**2  !* (4d0*pi*alpha_QED) * (m_top/(2d0*sitW*M_W))**2


RETURN
END SUBROUTINE







SUBROUTINE EvalAmp_QQB_TTBH(Mom,SqAmp)
use ModTopDecay
implicit none
real(8) :: Mom(1:4,1:13),SqAmp
complex(8) :: ResOffSh(1:4),Res(1:2)
complex(8) :: QuaPol(1:4,1:2,1:2)
integer :: hel4,TopHel1,TopHel2,nhel
real(8),parameter :: c_aa=8.0D0
integer, parameter :: inLeft=1,inRight=2,Hbos=3,tbar=4,t=5,  bbar=6,Wm=7,lepM=8,nubar=9,  b=10,Wp=11,lepP=12,nu=13
SqAmp = 0d0

     ExtParticles(1)%Mom(1:4) = Mom(1:4,tbar)
     ExtParticles(2)%Mom(1:4) = Mom(1:4,t)
     ExtParticles(5)%Mom(1:4) =-Mom(1:4,inLeft)
     ExtParticles(6)%Mom(1:4) =-Mom(1:4,inRight)
     ExtParticles(7)%Mom(1:4) = Mom(1:4,Hbos)

     if( TOPDECAYS.ne.0 ) then
        call TopDecay(ATop_,(/Mom(1:4,bbar),Mom(1:4,lepM),Mom(1:4,nubar)/),ExtParticles(1)%Pol(1:4))
        call TopDecay(Top_,(/Mom(1:4,b),Mom(1:4,lepP),Mom(1:4,nu)/),ExtParticles(2)%Pol(1:4))
     endif
     ExtParticles(7)%Pol(1:4) = 1d0
!    call HDecay(ExtParticles(7),DK_LO,MomExt(1:4,12:13))
     call ubarSpi_Dirac(ExtParticles(6)%Mom(1:4),0d0,-1,QuaPol(1:4,1,1))
     call ubarSpi_Dirac(ExtParticles(6)%Mom(1:4),0d0,+1,QuaPol(1:4,1,2))
     call    vSpi_Dirac(ExtParticles(5)%Mom(1:4),0d0,-1,QuaPol(1:4,2,1))
     call    vSpi_Dirac(ExtParticles(5)%Mom(1:4),0d0,+1,QuaPol(1:4,2,2))



     nhel=-1
     if( TOPDECAYS.EQ.0 ) nhel=+1
     do TopHel1=-1,nhel,2
     do TopHel2=-1,nhel,2
     if( TOPDECAYS.eq.0 ) then
             call ubarSpi_Dirac(ExtParticles(2)%Mom(1:4),M_Top,TopHel1,ExtParticles(2)%Pol(1:4))
             call    vSpi_Dirac(ExtParticles(1)%Mom(1:4),M_Top,TopHel2,ExtParticles(1)%Pol(1:4))
     endif
     do hel4=1,2

        ExtParticles(6)%Pol(1:4) = QuaPol(1:4,2,hel4)
!         ExtParticles(4)%Pol(1:4) = ExtParticles(4)%Mom(1:4); print *, "checking gauge invariance"
        call new_calc_ampl(0,0,TheTreeAmps_QQB_TTBH(1),ResOffSh(1:4))

        Res(1) = psp1_(ResOffSh(1:4),QuaPol(1:4,1,1))! hel+
        Res(2) = psp1_(ResOffSh(1:4),QuaPol(1:4,1,2))! hel-

        SqAmp = SqAmp   &
              + c_aa * dreal( Res(1)*dconjg(Res(1)) + Res(2)*dconjg(Res(2)) )
    enddo
    enddo
    enddo
    SqAmp = SqAmp * SpinAvg * QuarkColAvg**2 * (4d0*Pi*alphas)**2  !* (4d0*pi*alpha_QED) * (m_top/(2d0*sitW*M_W))**2

RETURN
END SUBROUTINE








SUBROUTINE new_calc_ampl(tag_f,tag_Z,TreeProc,Res)
implicit none
integer :: tag_f,tag_Z,n
complex(8) :: Res(1:Ds)
type(TreeProcess) :: TreeProc
logical,parameter :: Boson=.true.
integer :: i,j,Order(1:6)


      if( TreeProc%NumQua.eq.2 .and. TreeProc%NumSca.eq.0 ) then!  2 quarks and no scalars
          if ( TreeProc%PartType(1).eq.Glu_ .and. Boson ) then
             Res(1:Dv) = cur_g_2fV( TreeProc%Gluons(2:TreeProc%NumGlu(0)),TreeProc%Quarks(1:TreeProc%NumQua),TreeProc%Boson,TreeProc%NumGlu(0:3) )
          elseif( IsAQuark(TreeProc%PartType(1)) .and. Boson ) then
             if( TreeProc%NumV.eq.1 ) then
                Res(1:Ds) = cur_f_2fV(TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Quarks(2:2),TreeProc%Quarks(1)%PartType,TreeProc%Boson,TreeProc%NumGlu(0:2))
             else
                print *, "requested current with a boson is not available"
                stop
             endif
          else
             print *, "requested current is not available 2q"
             stop
          endif

      elseif( TreeProc%NumQua.eq.4  .and. TreeProc%NumSca.eq.0) then!  4 quarks, no scalars
          if( TreeProc%PartType(1).eq.Glu_ ) then
                 Res(1:Dv) = cur_g_4fV( TreeProc%Gluons(2:TreeProc%NumGlu(0)),TreeProc%Quarks(1:TreeProc%NumQua),TreeProc%Boson,TreeProc%BosonVertex,TreeProc%NumGlu(0:5) )
          elseif( IsAQuark(TreeProc%PartType(1)) ) then
                 Res(1:Ds) = cur_f_4fV( TreeProc%Gluons(1:TreeProc%NumGlu(0)),TreeProc%Quarks(2:4),TreeProc%Quarks(1)%PartType,TreeProc%Boson,TreeProc%BosonVertex,TreeProc%NumGlu(0:4),tag_f,tag_Z )
          else
             print *, "requested current is not available 4q"
             stop
          endif
      else
           print *, "requested current is not available xx"
           stop
      endif

return
END SUBROUTINE




SUBROUTINE LinkTreeParticles(TheTreeAmp,TheParticles)
implicit none
type(TreeProcess) :: TheTreeAmp
type(Particle),target :: TheParticles(:)
integer :: iPart,PartRef,PartType,ig,iq,ib,NPart,counterQ,counterG,LastQuark,QuarkPos(1:6)


            TheTreeAmp%NumW = 0
            TheTreeAmp%NumV = 0
            counterQ = 0
            counterG = 0

            do NPart=1,TheTreeAmp%NumPart
                  TheTreeAmp%PartType(NPart) = TheParticles( TheTreeAmp%PartRef(NPart) )%PartType
                  if( IsAQuark(TheTreeAmp%PartType(NPart)) ) then
                     !TheTreeAmp%NumQua = TheTreeAmp%NumQua + 1   ! this is suppoed to be done outside this subroutine
                     counterQ = counterQ + 1
                     QuarkPos(counterQ) = counterQ + counterG
                     LastQuark = counterQ! only required for BosonVertex below
!                   elseif( IsAScalar(TheTreeAmp%PartType(NPart)) ) then
! !                      TheTreeAmp%NumSca = TheTreeAmp%NumSca + 1   ! this is suppoed to be done outside this subroutine
!                      counterQ = counterQ + 1!     treat the scalar like a quark here because this is only to determine NumGlu
!                      QuarkPos(counterQ) = counterQ + counterG
                  elseif( TheTreeAmp%PartType(NPart).eq.Glu_ ) then
                     counterG = counterG + 1
                  elseif( IsABoson(TheTreeAmp%PartType(NPart)) ) then! careful: bosons should only be places *between* same flavor quark lines
                     if( NPart.eq.1 ) print *, "Vector boson should not be the first particle."
                     if( abs(TheTreeAmp%PartType(NPart)).eq.abs(Wp_) ) TheTreeAmp%NumW = TheTreeAmp%NumW + 1
                     if( abs(TheTreeAmp%PartType(NPart)).eq.abs(Z0_) ) TheTreeAmp%NumV = TheTreeAmp%NumV + 1
                     if( abs(TheTreeAmp%PartType(NPart)).eq.abs(Pho_)) TheTreeAmp%NumV = TheTreeAmp%NumV + 1
!                      TheTreeAmp%BosonVertex = TheTreeAmp%PartType(LastQuark)! this variable specifies the position of the boson wrt. to the quark lines
                     TheTreeAmp%BosonVertex = LastQuark
                  endif
            enddo


!           set number of gluons between quark lines
!            if( IsAQuark( TheTreeAmp%PartType(1) ) .or. IsAScalar(TheTreeAmp%PartType(1)) ) then ! not a gluon and not a boson
           if( IsAQuark( TheTreeAmp%PartType(1) ) ) then ! not a gluon and not a boson
             if( TheTreeAmp%NumQua+TheTreeAmp%NumSca .eq. 2 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(2) = TheTreeAmp%NumSca+TheTreeAmp%NumQua+TheTreeAmp%NumGlu(0) - QuarkPos(2)
             endif
             if( TheTreeAmp%NumQua .eq. 4 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(2) = QuarkPos(3) - QuarkPos(2) - 1
                   TheTreeAmp%NumGlu(3) = QuarkPos(4) - QuarkPos(3) - 1
                   TheTreeAmp%NumGlu(4) = TheTreeAmp%NumSca+TheTreeAmp%NumQua+TheTreeAmp%NumGlu(0) - QuarkPos(4)
             endif
             if( TheTreeAmp%NumQua .eq. 6 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(2) = QuarkPos(3) - QuarkPos(2) - 1
                   TheTreeAmp%NumGlu(3) = QuarkPos(4) - QuarkPos(3) - 1
                   TheTreeAmp%NumGlu(4) = QuarkPos(5) - QuarkPos(4) - 1
                   TheTreeAmp%NumGlu(5) = QuarkPos(6) - QuarkPos(5) - 1
                   TheTreeAmp%NumGlu(6) = TheTreeAmp%NumSca+TheTreeAmp%NumQua+TheTreeAmp%NumGlu(0) - QuarkPos(6)
             endif

           elseif( TheTreeAmp%PartType(1).eq.10 ) then ! is a gluon
             if( TheTreeAmp%NumQua .eq. 2 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(1) - 2
                   TheTreeAmp%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(3) = TheTreeAmp%NumSca+TheTreeAmp%NumQua+TheTreeAmp%NumGlu(0) - QuarkPos(2)
             endif
             if( TheTreeAmp%NumQua .eq. 4 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(1) - 2
                   TheTreeAmp%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(3) = QuarkPos(3) - QuarkPos(2) - 1
                   TheTreeAmp%NumGlu(4) = QuarkPos(4) - QuarkPos(3) - 1
                   TheTreeAmp%NumGlu(5) = TheTreeAmp%NumSca+TheTreeAmp%NumQua+TheTreeAmp%NumGlu(0) - QuarkPos(4)
             endif
             if( TheTreeAmp%NumQua .eq. 6 ) then
                   TheTreeAmp%NumGlu(1) = QuarkPos(1) - 2
                   TheTreeAmp%NumGlu(2) = QuarkPos(2) - QuarkPos(1) - 1
                   TheTreeAmp%NumGlu(3) = QuarkPos(3) - QuarkPos(2) - 1
                   TheTreeAmp%NumGlu(4) = QuarkPos(4) - QuarkPos(3) - 1
                   TheTreeAmp%NumGlu(5) = QuarkPos(5) - QuarkPos(4) - 1
                   TheTreeAmp%NumGlu(6) = QuarkPos(6) - QuarkPos(5) - 1
                   TheTreeAmp%NumGlu(7) = TheTreeAmp%NumSca+TheTreeAmp%NumQua+TheTreeAmp%NumGlu(0) - QuarkPos(6)
             endif
           endif


  ig=0; iq=0; ib=0;
  do iPart=1,TheTreeAmp%NumPart
     PartRef = TheTreeAmp%PartRef(iPart)
     PartType= TheParticles(PartRef)%PartType
     if( PartType.eq.Glu_ ) then
           ig=ig+1
           TheTreeAmp%Gluons(ig)%PartType => TheParticles(PartRef)%PartType
           TheTreeAmp%Gluons(ig)%ExtRef   => TheParticles(PartRef)%ExtRef
           TheTreeAmp%Gluons(ig)%Mass     => TheParticles(PartRef)%Mass
           TheTreeAmp%Gluons(ig)%Mass2    => TheParticles(PartRef)%Mass2
           TheTreeAmp%Gluons(ig)%Helicity => TheParticles(PartRef)%Helicity
           TheTreeAmp%Gluons(ig)%Mom      => TheParticles(PartRef)%Mom
           TheTreeAmp%Gluons(ig)%Pol      => TheParticles(PartRef)%Pol
           if( PartType.ne.TheParticles(PartRef)%PartType ) print *,"Error1 in LinkTreeParticles"
     elseif( IsAQuark(PartType) ) then ! PartType==Quark
           iq=iq+1
           TheTreeAmp%Quarks(iq)%PartType => TheParticles(PartRef)%PartType
           TheTreeAmp%Quarks(iq)%ExtRef   => TheParticles(PartRef)%ExtRef
           TheTreeAmp%Quarks(iq)%Mass     => TheParticles(PartRef)%Mass
           TheTreeAmp%Quarks(iq)%Mass2    => TheParticles(PartRef)%Mass2
           TheTreeAmp%Quarks(iq)%Helicity => TheParticles(PartRef)%Helicity
           TheTreeAmp%Quarks(iq)%Mom      => TheParticles(PartRef)%Mom
           TheTreeAmp%Quarks(iq)%Pol      => TheParticles(PartRef)%Pol
           if( PartType.ne.TheParticles(PartRef)%PartType ) print *,"Error2 in LinkTreeParticles"
     elseif( IsABoson(PartType) ) then  ! PartType==Boson
           ib=ib+1
           if( ib.ge.2 ) print *, "Too many bosons in LinkTreeParticles"
           TheTreeAmp%Boson%PartType => TheParticles(PartRef)%PartType
           TheTreeAmp%Boson%ExtRef   => TheParticles(PartRef)%ExtRef
           TheTreeAmp%Boson%Mass     => TheParticles(PartRef)%Mass
           TheTreeAmp%Boson%Mass2    => TheParticles(PartRef)%Mass2
           TheTreeAmp%Boson%Helicity => TheParticles(PartRef)%Helicity
           TheTreeAmp%Boson%Mom      => TheParticles(PartRef)%Mom
           TheTreeAmp%Boson%Pol      => TheParticles(PartRef)%Pol
           if( PartType.ne.TheParticles(PartRef)%PartType ) print *,"Error2 in LinkTreeParticles"
     endif
  enddo
  if( ig.ne.TheTreeAmp%NumGlu(0) .OR. iq.ne.TheTreeAmp%NumQua+TheTreeAmp%NumSca .OR. ib.ne.TheTreeAmp%NumPart-TheTreeAmp%NumGlu(0)-TheTreeAmp%NumQua-TheTreeAmp%NumSca) print *,"Error3 in LinkTreeParticles"


RETURN
END SUBROUTINE



!SUBROUTINE UnLinkTreeParticles(TheTreeAmp)
!implicit none
!type(TreeProcess) :: TheTreeAmp
!integer :: iPart,PartType,PartRef,ig,iq,ib

!  ig=0; iq=0; ib=0;
!  do iPart=1,TheTreeAmp%NumPart
!     PartRef = TheTreeAmp%PartRef(iPart)
!     PartType = TheParticles(PartRef)%PartType

!     if( PartType.eq.Glu_ ) then
!           ig=ig+1
!           nullify(TheTreeAmp%Gluons(ig)%PartType)
!           nullify(TheTreeAmp%Gluons(ig)%ExtRef)
!           nullify(TheTreeAmp%Gluons(ig)%Mass)
!           nullify(TheTreeAmp%Gluons(ig)%Mass2)
!           nullify(TheTreeAmp%Gluons(ig)%Helicity)
!           nullify(TheTreeAmp%Gluons(ig)%Mom)
!           nullify(TheTreeAmp%Gluons(ig)%Pol)
!     elseif( IsAQuark(PartType) ) then ! PartType==Quark
!           iq=iq+1
!           nullify(TheTreeAmp%Quarks(iq)%PartType)
!           nullify(TheTreeAmp%Quarks(iq)%ExtRef)
!           nullify(TheTreeAmp%Quarks(iq)%Mass)
!           nullify(TheTreeAmp%Quarks(iq)%Mass2)
!           nullify(TheTreeAmp%Quarks(iq)%Helicity)
!           nullify(TheTreeAmp%Quarks(iq)%Mom)
!           nullify(TheTreeAmp%Quarks(iq)%Pol)
!     elseif( IsABoson(PartType) ) then  ! PartType==Boson
!           ib=ib+1
!           nullify(TheTreeAmp%Boson%PartType)
!           nullify(TheTreeAmp%Boson%ExtRef)
!           nullify(TheTreeAmp%Boson%Mass)
!           nullify(TheTreeAmp%Boson%Mass2)
!           nullify(TheTreeAmp%Boson%Helicity)
!           nullify(TheTreeAmp%Boson%Mom)
!           nullify(TheTreeAmp%Boson%Pol)
!     endif
!  enddo

!RETURN
!END SUBROUTINE




! ----------------------------------------------------


FUNCTION cur_f_2fV(Gluons,Quark,Quark1PartType,Boson,NumGlu) result(Res)           ! Quarks(:) DOES include the OFF-shell quark, in contrast to all other routines!
implicit none
complex(8) :: Res(1:Ds)
integer :: NumGlu(0:2),i,rIn,rOut,Quark1PartType
type(PtrToParticle) :: Gluons(1:),Boson,Quark(2:2)
complex(8) :: GluMom(1:Dv,NumGlu(0)), QuarkMom(1:Dv)
complex(8) :: GluPol(1:Dv,NumGlu(0)), QuarkPol(1:Ds)



   do i=1,NumGlu(0)
    GluMom(1:Dv,i) = Gluons(i)%Mom(1:Dv)
    GluPol(1:Dv,i) = Gluons(i)%Pol(1:Dv)
   enddo
   QuarkMom(1:Dv) = Quark(2)%Mom(1:Dv)
   QuarkPol(1:Ds) = Quark(2)%Pol(1:Ds)

   if( Quark1PartType.ne.-Quark(2)%PartType ) print *, "Wrong quark flavors in cur_f_2fV"
   if( NumGlu(0)-NumGlu(1)-NumGlu(2).ne.0 ) print *, "Wrong NumGlu in cur_f_2fV",NumGlu(0)-NumGlu(1)-NumGlu(2)

   rIn =1
   rOut=NumGlu(0)
   if( Quark(2)%PartType .gt.0 ) then      !    X----->----
      Res(:) = fV(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),Quark(2)%Mass,Quark1PartType,Boson%Pol(1:Dv),Boson%Mom(1:Dv),NumGlu(1))
   else
      Res(:) = bfV(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),Quark(2)%Mass,Quark1PartType,Boson%Pol(1:Dv),Boson%Mom(1:Dv),NumGlu(1))
   endif


return
END FUNCTION







      recursive function fV(e,k,sp,p,mass,QuarkFlavor,eV,kV,ms) result(res)
      implicit none
      complex(8), intent(in) :: e(:,:), k(:,:)
      complex(8), intent(in) :: sp(:), p(:)
      complex(8), intent(in) :: eV(:), kV(:)
      integer, intent(in) ::  ms,QuarkFlavor
      integer             :: ms1,m,ng1, ng2
      integer :: ngluon
      complex(8)             :: res(size(sp))
      complex(8)             :: tmp(size(sp))
      complex(8)             :: k1(size(p))
      complex(8)             :: k2(size(p))
      complex(8)             :: sp2(size(sp))
      complex(8)             :: sp3(size(sp))
      complex(8)             :: e1(size(e,dim=1))
      complex(8)             :: e2(size(e,dim=1))
      complex(8)  :: k1sq,k2sq,k3sq
      real(8) :: mass
      complex(8) :: couplVQQ_left,couplVQQ_right,couplVQQ_left2,couplVQQ_right2
      character,parameter :: FerFla*3="dum" ! dummy, only used for check of flavor consistency inside the functions f,bf


      ngluon = size(e,dim=2)
      ng1 = ms   !#gluons to the left of a f-line
      ng2 = ngluon - ms  !#gluons to the right of the f-line
      if (ng2 < 0) write(*,*) 'WRONG DEFINITION OF CURRENT fV'

      if( abs(QuarkFlavor).eq.Top_ .or. abs(QuarkFlavor).eq.Bot_ ) then!   note that Bot_ is treated as top quark in TOPAZ!
         couplVQQ_left  = couplHTT_left_dyn
         couplVQQ_right = couplHTT_right_dyn
      else
         couplVQQ_left=0d0
         couplVQQ_right=0d0
         print *, "this should not happen for the Higgs",QuarkFlavor,Top_
      endif



if (ngluon == 0) then
         res = vbqV(sp,eV,couplVQQ_left,couplVQQ_right)
else

       res = (0d0,0d0)
       do m=0,ng2-1
           k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
           e1=g(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon))
           k1sq=sc_(k1,k1)

           k2 = sum(k(:,1:ng1+m),dim=2)
           k2 = k2 + p + kV
           k2sq = sc_(k2,k2)-mass**2
           sp2 = fV(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,QuarkFlavor,eV,kV,ng1)
           sp2 = spb2_(sp2,k2)+mass*sp2

           tmp = vqg(sp2,e1)
           if (m < ng2-1)  then
              if(abs(k1sq) > propcut) then
                  tmp = -(0d0,1d0)/k1sq*tmp
              else
                  tmp = (0d0,0d0)
              endif
           endif

           if (abs(k2sq) > propcut) then
                  tmp =  (0d0,1d0)/k2sq*tmp
           else
                  tmp = (0d0,0d0)
           endif
           res = res + tmp
        enddo




        do m=1,ng1
           k1 = sum(k(:,1:m),dim=2)
           e1=g(e(:,1:m),k(:,1:m))
           k1sq = sc_(k1,k1)

           k2 = sum(k(:,m+1:ngluon),dim=2)
           k2 = k2 + p + kV
           k2sq = sc_(k2,k2) - mass**2
           ms1 = ng1 - m
           sp2=fV(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,QuarkFlavor,eV,kV,ms1)

           sp2 = spb2_(sp2,k2)+mass*sp2

           tmp = vgq(e1,sp2)
           if (m > 1) then
              if (abs(k1sq) > propcut) then
                  tmp=-(0d0,1d0)/k1sq*tmp
              else
                  tmp = (0d0,0d0)
              endif
           endif

           if (abs(k2sq) > propcut) then
              tmp=(0d0,1d0)/k2sq*tmp
           else
              tmp = (0d0,0d0)
           endif

           res = res + tmp
        enddo


        sp2 = f(e,k,sp,p,mass,FerFla,FerFla,ms)
        k2 = sum(k(:,1:ngluon),dim=2)
        k2 = k2 + p
        k2sq = sc_(k2,k2)  - mass**2

        sp2 = spb2_(sp2,k2)+ mass*sp2

        tmp = vbqV(sp2,eV,couplVQQ_left,couplVQQ_right)
        if (abs(k2sq) > propcut) then
               tmp = (0d0,1d0)/k2sq*tmp
        else
                tmp = (0d0,0d0)
        endif
        res = res + tmp

endif

end function fV




      recursive function bfV(e,k,sp,p,mass,QuarkFlavor,eV,kV,ms) result(res)
      implicit none
      complex(8), intent(in) :: e(:,:), k(:,:)
      complex(8), intent(in) :: sp(:), p(:)
      complex(8), intent(in) :: eV(:), kV(:)
      integer, intent(in) ::  ms,QuarkFlavor
      integer             :: ms1,m,ng1, ng2
      integer :: ngluon
      complex(8)             :: res(size(sp))
      complex(8)             :: tmp(size(sp))
      complex(8)             :: k1(size(p))
      complex(8)             :: k2(size(p))
      complex(8)             :: sp2(size(sp))
      complex(8)             :: sp3(size(sp))
      complex(8)             :: e1(size(e,dim=1))
      complex(8)             :: e2(size(e,dim=1))
      complex(8)  :: k1sq,k2sq,k3sq
      real(8) :: mass
      complex(8) :: couplVQQ_left,couplVQQ_right,couplVQQ_left2,couplVQQ_right2
      character,parameter :: FerFla*3="dum" ! dummy, only used for check of flavor consistency inside the functions f,bf


      ngluon = size(e,dim=2)
      ng1 = ms   !#gluons to the left of a f-line
      ng2 = ngluon - ms  !#gluons to the right of the f-line

      if (ng2 < 0) write(*,*) 'WRONG DEFINITION OF CURRENT fbV'

      if( abs(QuarkFlavor).eq.Top_ .or. abs(QuarkFlavor).eq.Bot_ ) then!   note that Bot_ is treated as top quark in TOPAZ!
          couplVQQ_left  = couplHTT_left_dyn
          couplVQQ_right = couplHTT_right_dyn
      else
          couplVQQ_left=0d0
          couplVQQ_right=0d0
          print *, "this should not happen for the Higgs",QuarkFlavor,Top_
      endif



if (ngluon == 0) then
         res = vVq(eV,sp,couplVQQ_left,couplVQQ_right)
else

       res = (0d0,0d0)
       do m=0,ng2-1
           k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
           e1=g(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon))
           k1sq=sc_(k1,k1)

           k2 = sum(k(:,1:ng1+m),dim=2)
           k2 = -k2 - p - kV
           k2sq = sc_(k2,k2)-mass**2

           sp2 = bfV(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,QuarkFlavor,eV,kV,ng1)
           sp2 = spi2_(k2,sp2)+mass*sp2

           tmp = vbqg(sp2,e1)

           if (m < ng2-1) then
            if (abs(k1sq) > propcut) then
                tmp = -(0d0,1d0)/k1sq*tmp
            else
                tmp = (0d0,0d0)
            endif
           endif

           if (abs(k2sq) > propcut) then
              tmp =  (0d0,1d0)/k2sq*tmp
           else
               tmp = (0d0,0d0)
           endif

           res = res + tmp
        enddo


        do m=1,ng1

           k1 = sum(k(:,1:m),dim=2)
           e1=g(e(:,1:m),k(:,1:m))
           k1sq = sc_(k1,k1)

           k2 = sum(k(:,m+1:ngluon),dim=2)
           k2 = -k2 - p - kV
           k2sq = sc_(k2,k2) - mass**2
           ms1 = ng1 - m
           sp2=bfV(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,QuarkFlavor,eV,kV,ms1)

           sp2 = spi2_(k2,sp2)+mass*sp2

           tmp = vgbq(e1,sp2)

           if (m > 1) then
              if (abs(k1sq) > propcut) then
                  tmp=-(0d0,1d0)/k1sq*tmp
              else
                  tmp = (0d0,0d0)
              endif
           endif

           if (abs(k2sq) > propcut) then
                 tmp=(0d0,1d0)/k2sq*tmp
              else
                 tmp = (0d0,0d0)
           endif

           res = res + tmp

        enddo



        sp2 = bf(e,k,sp,p,mass,FerFla,FerFla,ms)
        k2 = sum(k(:,1:ngluon),dim=2)
        k2 = -k2 - p
        k2sq = sc_(k2,k2)   -mass**2

!        sp2 = spb2_(sp2,k2) +mass*sp2
        sp2 = spi2_(k2,sp2) +mass*sp2

        tmp = vVq(eV,sp2,couplVQQ_left,couplVQQ_right)
        if (abs(k2sq) > propcut) then
               tmp = (0d0,1d0)/k2sq*tmp
        else
               tmp = (0d0,0d0)
        endif
        res = res + tmp

endif




      end function bfV





FUNCTION cur_f_4fV(Gluons,Quarks,Quark1PartType,Boson,BosonVertex,NumGlu,tag_f,tag_Z) result(res)           ! Quarks(:) does not include the OFF-shell quark
implicit none
integer :: NumGlu(0:4),Quark1PartType
type(PtrToParticle) :: Gluons(1:),Quarks(2:4),Boson
integer :: tag_f,BosonVertex,tag_Z
integer,target :: TmpExtRef
complex(8) :: res(1:Ds),tmp(1:Ds)
complex(8) :: ubar1(1:Ds)
complex(8),target :: ubar0(1:Ds)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpQuark(1:1)
complex(8) :: PropFac1,PropFac2
complex(8),target :: pmom1(1:Dv)
complex(8) :: pmom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4).ne.0 ) print *, "wrong number of gluons in cur_f_4f"
    if(Quarks(3)%PartType.eq.-Quarks(4)%PartType .and. Quark1PartType.ne.-Quarks(2)%PartType ) print *,"wrong flavor in cur_f_4f (1)"
    if(Quarks(2)%PartType.eq.-Quarks(3)%PartType .and. Quark1PartType.ne.-Quarks(4)%PartType ) print *,"wrong flavor in cur_f_4f (2)"
!DEC$ ENDIF
    Res(:)=(0d0,0d0)
   if ( Quark1PartType.eq.-Quarks(2)%PartType .and. Quarks(3)%PartType.eq.-Quarks(4)%PartType ) then
      do n2a=0,NumGlu(2)
         do n4a=0,NumGlu(4)
            n2b = NumGlu(2)-n2a
            n4b = NumGlu(4)-n4a
            rIn =NumGlu(1)+n2a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            if (BosonVertex .eq. 1 .or. BosonVertex .eq. 2 .or. BosonVertex .eq. 4) then
               Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(3:4),(/1+n2b+NumGlu(3)+n4a,n2b,NumGlu(3),n4a/))
               PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom
            elseif (BosonVertex .eq. 3) then
               Eps2 = cur_g_2fV(Gluons(rIn:rOut),Quarks(3:4),Boson,(/1+n2b+NumGlu(3)+n4a,n2b,NumGlu(3),n4a/))
               PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom(:) + Quarks(4)%Mom(:) + Boson%Mom(:)
            endif
            PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
            if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1
            do n1a=0,NumGlu(1)
               n1b = NumGlu(1)-n1a
               ! Fer2
               if (BosonVertex.eq.1 .or. BosonVertex.eq.2) then
                  rIn =n1a+1
                  rOut=NumGlu(1)+n2a
                  ubar1(:) = cur_f_2fV(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,Boson,(/n2a+n1b,n1b,n2a/) )
                  PMom2(:) = Quarks(2)%Mom(:) + SumMom(Gluons,rIn,rOut) + Boson%Mom(:)
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)

                  if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) then
                     PropFac2=(0d0,0d0)
                  endif

                  if( Quarks(2)%PartType.lt.0 ) then
                     ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(2)%Mass*ubar1(:))*PropFac2
                  else
                     ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(2)%Mass*ubar1(:))*PropFac2
                  endif

                  if( Quarks(2)%PartType.lt.0 ) then
                     ubar0(:) = vbqg(ubar1,eps2)
                  else
                     ubar0(:) = vqg(ubar1,eps2)
                  endif

                  PMom1 = Quarks(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a) + Boson%Mom(:)
                  if(n1a.ge.1 .or. n4b.ge.1) then
                     PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
                     if( abs(sc_(PMom1,PMom1)-Quarks(2)%Mass2).lt.PropCut ) then
                        PropFac1=(0d0,0d0)
                     endif
                     if( Quarks(2)%PartType.lt.0 ) then
                        ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(2)%Mass*ubar0(:))*PropFac1
                     else
                        ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(2)%Mass*ubar0(:))*PropFac1
                     endif
                  endif

                  TmpQuark(1)%Mom  => PMom1(:)
                  TmpQuark(1)%Pol  => ubar0(:)
                  TmpQuark(1)%Mass => Quarks(2)%Mass
                  TmpQuark(1)%Mass2=> Quarks(2)%Mass2
                  TmpExtRef = -1
                  TmpQuark(1)%ExtRef => TmpExtRef
                  TmpQuark(1)%PartType => Quarks(2)%PartType
                  counter=1
                  rIn =1
                  rOut=n1a
                  do i=rIn,rOut
                     call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                     counter=counter+1
                  enddo
                  rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut=NumGlu(0)
                  do i=rIn,rOut
                     call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                     counter=counter+1
                  enddo
                  tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
                  Res(:) = Res(:) + tmp(:)
                endif

                ! Fer 3
                if (BosonVertex .eq. 1 .or. BosonVertex .eq. 4) then
                   rIn =n1a+1
                   rOut=NumGlu(1)+n2a
                   ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n2a+n1b,n1b,n2a/) )
                   if(n1b.ge.1 .or. n2a.ge.1) then
                      PMom2(:) = Quarks(2)%Mom + SumMom(Gluons,rIn,rOut)  ! can be moved outside the n1a-loop
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
                      if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) cycle
                      if( Quarks(2)%PartType.lt.0 ) then
                         ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(2)%Mass*ubar1(:))*PropFac2
                      else
                         ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(2)%Mass*ubar1(:))*PropFac2
                      endif
                   endif
                   if( Quarks(2)%PartType.lt.0 ) then
                      ubar0(:) = vbqg(ubar1,eps2)
                   else
                      ubar0(:) = vqg(ubar1,eps2)
                   endif

                   PMom1 = Quarks(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a)
                   PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
                   if( abs(sc_(PMom1,PMom1)-Quarks(2)%Mass2).lt.PropCut ) then
                      PropFac1=(0d0,0d0)
                   endif
                   if( Quarks(2)%PartType.lt.0 ) then
                      ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(2)%Mass*ubar0(:))*PropFac1
                   else
                      ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(2)%Mass*ubar0(:))*PropFac1
                   endif

                   TmpQuark(1)%Mom  => PMom1(:)
                   TmpQuark(1)%Pol  => ubar0(:)
                   TmpQuark(1)%Mass => Quarks(2)%Mass
                   TmpQuark(1)%Mass2=> Quarks(2)%Mass2
                   TmpExtRef = -1
                   TmpQuark(1)%ExtRef => TmpExtRef
                   TmpQuark(1)%PartType => Quarks(2)%PartType
                   counter=1
                   rIn =1
                   rOut=n1a
                   do i=rIn,rOut
                      call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                      counter=counter+1
                   enddo
                   rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                   rOut=NumGlu(0)
                   do i=rIn,rOut
                      call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                      counter=counter+1
                   enddo
                   tmp(:) = cur_f_2fV(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,Boson,(/counter-1,n1a,n4b/) )
                   Res(:) = Res(:) + tmp(:)
                endif

                if (BosonVertex .eq. 3) then
                   rIn =n1a+1
                   rOut=NumGlu(1)+n2a
                   ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n2a+n1b,n1b,n2a/) )
                   if(n1b.ge.1 .or. n2a.ge.1) then
                      PMom2(:) = Quarks(2)%Mom + SumMom(Gluons,rIn,rOut)  ! can be moved outside the n1a-loop
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
                      if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) then
                         PropFac2=(0d0,0d0)
                      endif
                      if( Quarks(2)%PartType.lt.0 ) then
                         ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(2)%Mass*ubar1(:))*PropFac2
                      else
                         ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(2)%Mass*ubar1(:))*PropFac2
                      endif
                   endif
                   if( Quarks(2)%PartType.lt.0 ) then
                      ubar0(:) = vbqg(ubar1,eps2)
                   else
                      ubar0(:) = vqg(ubar1,eps2)
                   endif

                   PMom1 = Quarks(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a)+Boson%Mom
                   if(n1a.ge.1 .or. n4b.ge.1) then
                      PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
                      if( abs(sc_(PMom1,PMom1)-Quarks(2)%Mass2).lt.PropCut ) cycle
                      if( Quarks(2)%PartType.lt.0 ) then
                         ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(2)%Mass*ubar0(:))*PropFac1
                      else
                         ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(2)%Mass*ubar0(:))*PropFac1
                      endif
                   endif

                   TmpQuark(1)%Mom  => PMom1(:)
                   TmpQuark(1)%Pol  => ubar0(:)
                   TmpQuark(1)%Mass => Quarks(2)%Mass
                   TmpQuark(1)%Mass2=> Quarks(2)%Mass2
                   TmpExtRef = -1
                   TmpQuark(1)%ExtRef => TmpExtRef
                   TmpQuark(1)%PartType => Quarks(2)%PartType
                   counter=1
                   rIn =1
                   rOut=n1a
                   do i=rIn,rOut
                      call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                      counter=counter+1
                   enddo
                   rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                   rOut=NumGlu(0)
                   do i=rIn,rOut
                      call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                      counter=counter+1
                   enddo
                   tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
                   Res(:) = Res(:) + tmp(:)
                endif
             enddo
          enddo
       enddo
    endif


    if( (Quark1PartType.eq.-Quarks(4)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType)      .AND.  &
         (Quarks(4)%ExtRef.ne.-1 .or. tag_f.ne.1 .or. abs(Quark1PartType).ne.abs(Quarks(2)%PartType)) ) then

       do n1a=0,NumGlu(1)
          do n3a=0,NumGlu(3)
             n1b = NumGlu(1)-n1a
             n3b = NumGlu(3)-n3a

             rIn =n1a+1
             rOut=NumGlu(1)+NumGlu(2)+n3a

! This means that all the ext gluons are on this lines,
! and we must remove this to prevent color issues with a Z on the quark loop, see RR notes

             if (BosonVertex.eq.1 .or. BosonVertex.eq.3 .or. BosonVertex.eq.4)  then
                if ( tag_Z .eq. 1 .and. n1b+NumGlu(2)+n3a .eq. NumGlu(0) ) then
!                     print  * , "cycle for tag_Z=1 in cur_f_4fV"
                    cycle
                endif
                Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(2:3),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
                PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
             elseif (BosonVertex .eq. 2) then
                Eps2 = cur_g_2fV(Gluons(rIn:rOut),Quarks(2:3),Boson,(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
                PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Boson%Mom
             endif
             PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
             if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
             Eps2 = Eps2*PropFac1

             do n4a=0,NumGlu(4)
                n4b = NumGlu(4)-n4a
                ! radiate V off Fer4
                if (BosonVertex .eq. 3 .or. BosonVertex .eq. 4) then
                   rIn =NumGlu(1)+NumGlu(2)+n3a+1
                   rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a

                   ubar1(:) = cur_f_2fV(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,Boson,(/n3b+n4a,n3b,n4a/) )
                   PMom2(:) = Quarks(4)%Mom(:) + SumMom(Gluons,rIn,rOut) + Boson%Mom(:)
                   PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(4)%Mass2)
                   if( abs(sc_(PMom2,PMom2)-Quarks(4)%Mass2).lt.PropCut ) then
                      PropFac2=(0d0,0d0)
                   endif

                   if( Quarks(4)%PartType.lt.0 ) then
                      ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(4)%Mass*ubar1(:))*PropFac2
                   else
                      ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(4)%Mass*ubar1(:))*PropFac2
                   endif
                   if( Quarks(4)%PartType.lt.0 ) then
                      ubar0(:) = vgbq(eps2,ubar1)
                   else
                      ubar0(:) = vgq(eps2,ubar1)
                   endif

                   PMom1 = Quarks(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a) + Boson%Mom(:)
                   if(n1a.ge.1 .or. n4b.ge.1) then
                      PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(4)%Mass2)
                      if( abs(sc_(PMom1,PMom1)-Quarks(4)%Mass2).lt.PropCut ) then
                         PropFac1=(0d0,0d0)
                      endif
                      if( Quarks(4)%PartType.lt.0 ) then
                         ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(4)%Mass*ubar0(:))*PropFac1
                      else
                         ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(4)%Mass*ubar0(:))*PropFac1
                      endif
                   endif

                   TmpQuark(1)%Mom  => PMom1(:)
                   TmpQuark(1)%Pol  => ubar0(:)
                   TmpQuark(1)%Mass => Quarks(4)%Mass
                   TmpQuark(1)%Mass2=> Quarks(4)%Mass2
                   TmpExtRef = -1
                   TmpQuark(1)%ExtRef => TmpExtRef
                   TmpQuark(1)%PartType => Quarks(4)%PartType
                   counter=1
                   rIn =1
                   rOut=n1a
                   do i=rIn,rOut
                      call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                      counter=counter+1
                   enddo
                   rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                   rOut=NumGlu(0)
                   do i=rIn,rOut
                      call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                      counter=counter+1
                   enddo
                   tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
                   Res(:) = Res(:) + tmp(:)
                endif

                if (BosonVertex .eq. 1 .or. BosonVertex .eq. 4) then
                   ! radiate V off Fer1
                   rIn =NumGlu(1)+NumGlu(2)+n3a+1
                   rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
                   ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/n3b+n4a,n3b,n4a/) )
                   if(n3b.ge.1 .or. n4a.ge.1) then
                      PMom2(:) = Quarks(4)%Mom + SumMom(Gluons,rIn,rOut)  ! can be moved outside the n1a-loop
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2)-Quarks(4)%Mass2).lt.PropCut ) then
                         PropFac2=(0d0,0d0)
                      endif

                      if( Quarks(4)%PartType.lt.0 ) then
                         ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(4)%Mass*ubar1(:))*PropFac2
                      else
                         ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(4)%Mass*ubar1(:))*PropFac2
                      endif
                   endif
                   if( Quarks(4)%PartType.lt.0 ) then
                      ubar0(:) = vgbq(eps2,ubar1)
                   else
                      ubar0(:) = vgq(eps2,ubar1)
                   endif

                   PMom1 = Quarks(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a)
                   PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(4)%Mass2)
                   if( abs(sc_(PMom1,PMom1)-Quarks(4)%Mass2).lt.PropCut ) then
                      PropFac1=(0d0,0d0)
                   endif
                   if( Quarks(4)%PartType.lt.0 ) then
                      ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(4)%Mass*ubar0(:))*PropFac1
                   else
                      ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(4)%Mass*ubar0(:))*PropFac1
                   endif

                   TmpQuark(1)%Mom  => PMom1(:)
                   TmpQuark(1)%Pol  => ubar0(:)
                   TmpQuark(1)%Mass => Quarks(4)%Mass
                   TmpQuark(1)%Mass2=> Quarks(4)%Mass2
                   TmpExtRef = -1
                   TmpQuark(1)%ExtRef => TmpExtRef
                   TmpQuark(1)%PartType => Quarks(4)%PartType
                   counter=1
                   rIn =1
                   rOut=n1a
                   do i=rIn,rOut
                      call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                      counter=counter+1
                   enddo
                   rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                   rOut=NumGlu(0)
                   do i=rIn,rOut
                      call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                      counter=counter+1
                   enddo
                   tmp(:) = cur_f_2fV(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,Boson,(/counter-1,n1a,n4b/) )
                   Res(:) = Res(:) + tmp(:)
                endif


                if (BosonVertex .eq. 2) then
                   rIn =NumGlu(1)+NumGlu(2)+n3a+1
                   rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
                   ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/n3b+n4a,n3b,n4a/) )
                   if ( n3b .ge. 1 .or. n4a .ge. 1) then
                      PMom2(:) = Quarks(4)%Mom(:) + SumMom(Gluons,rIn,rOut)
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(4)%Mass2)

                      if( abs(sc_(PMom2,PMom2)-Quarks(4)%Mass2).lt.PropCut ) then
                         PropFac2=(0d0,0d0)
                      endif

                      if( Quarks(4)%PartType.lt.0 ) then
                         ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(4)%Mass*ubar1(:))*PropFac2
                      else
                         ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(4)%Mass*ubar1(:))*PropFac2
                      endif
                   endif
                   if( Quarks(4)%PartType.lt.0 ) then
                      ubar0(:) = vgbq(eps2,ubar1)
                   else
                      ubar0(:) = vgq(eps2,ubar1)
                   endif
                   PMom1 = Quarks(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a)+Boson%Mom
                   if(n1a.ge.1 .or. n4b.ge.1) then
                      PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(4)%Mass2)
                      if( abs(sc_(PMom1,PMom1)-Quarks(4)%Mass2).lt.PropCut ) then
                         cycle
                      endif
                      if( Quarks(4)%PartType.lt.0 ) then
                         ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(4)%Mass*ubar0(:))*PropFac1
                      else
                         ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(4)%Mass*ubar0(:))*PropFac1
                      endif
                   endif

                   TmpQuark(1)%Mom  => PMom1(:)
                   TmpQuark(1)%Pol  => ubar0(:)
                   TmpQuark(1)%Mass => Quarks(4)%Mass
                   TmpQuark(1)%Mass2=> Quarks(4)%Mass2
                   TmpExtRef = -1
                   TmpQuark(1)%ExtRef => TmpExtRef
                   TmpQuark(1)%PartType => Quarks(4)%PartType
                   counter=1
                   rIn =1
                   rOut=n1a
                   do i=rIn,rOut
                      call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                      counter=counter+1
                   enddo
                   rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                   rOut=NumGlu(0)
                   do i=rIn,rOut
                      call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                      counter=counter+1
                   enddo
                   tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
                   Res(:) = Res(:) + tmp(:)
                endif
             enddo
          enddo
       enddo

    endif

    return
  END FUNCTION cur_f_4fV






FUNCTION cur_f_4f(Gluons,Quarks,Quark1PartType,NumGlu,tag_f,tag_Z_arg) result(res)           ! Quarks(:) does not include the OFF-shell quark
implicit none
integer :: NumGlu(0:4),Quark1PartType
type(PtrToParticle) :: Gluons(1:),Quarks(2:4)
integer :: tag_f,tag_Z
integer, optional :: tag_Z_arg
integer,target :: TmpExtRef
complex(8) :: res(1:Ds),tmp(1:Ds)
complex(8) :: ubar1(1:Ds)
complex(8),target :: ubar0(1:Ds)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(4)),TmpQuark(1:1)
complex(8) :: PropFac1,PropFac2
complex(8),target :: pmom1(1:Dv)
complex(8) :: pmom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4).ne.0 ) print *, "wrong number of gluons in cur_f_4f"
    if(Quarks(3)%PartType.eq.-Quarks(4)%PartType .and. Quark1PartType.ne.-Quarks(2)%PartType ) print *,"wrong flavor in cur_f_4f (1)"
    if(Quarks(2)%PartType.eq.-Quarks(3)%PartType .and. Quark1PartType.ne.-Quarks(4)%PartType ) print *,"wrong flavor in cur_f_4f (2)"
!DEC$ ENDIF

   Res(:)=(0d0,0d0)

   tag_Z=0
   if( present(tag_Z_arg) ) then
      tag_Z=tag_Z_arg
   endif

   if( Quark1PartType.eq.-Quarks(2)%PartType .and. Quarks(3)%PartType.eq.-Quarks(4)%PartType ) then
!     (I)
      do n2a=0,NumGlu(2)
      do n4a=0,NumGlu(4)
         n2b = NumGlu(2)-n2a
         n4b = NumGlu(4)-n4a

         rIn =NumGlu(1)+n2a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
         Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(3:4),(/1+n2b+NumGlu(3)+n4a,n2b,NumGlu(3),n4a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1
         do n1a=0,NumGlu(1)
            n1b = NumGlu(1)-n1a
            ! Fer2
            rIn =n1a+1
            rOut=NumGlu(1)+n2a
            ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n2a+n1b,n1b,n2a/) )
            if(n1b.ge.1 .or. n2a.ge.1) then
               PMom2(:) = Quarks(2)%Mom + SumMom(Gluons,rIn,rOut)  ! can be moved outside the n1a-loop
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(2)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(2)%Mass*ubar1(:))*PropFac2
               endif
            endif
            if( Quarks(2)%PartType.lt.0 ) then
               ubar0(:) = vbqg(ubar1,eps2)       ! re-checked
            else
               ubar0(:) = vqg(ubar1,eps2)        ! re-checked
            endif

            PMom1 = Quarks(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a)
            if(n1a.ge.1 .or. n4b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(2)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(2)%Mass*ubar0(:))*PropFac1
               endif
            endif

            TmpQuark(1)%Mom  => PMom1(:)
            TmpQuark(1)%Pol  => ubar0(:)
            TmpQuark(1)%Mass => Quarks(2)%Mass
            TmpQuark(1)%Mass2=> Quarks(2)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(2)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
            Res(:) = Res(:) + tmp(:)
         enddo
      enddo
      enddo
   endif


!    if( (Quark1PartType.eq.-Quarks(4)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType) .AND.  &
!        ((abs(Quark1PartType).ne.abs(Quarks(2)%PartType).and.(tag_f.ne.3)) .or. &
!         (abs(Quark1PartType).eq.abs(Quarks(2)%PartType).and.(tag_f.ne.1.and.tag_f.ne.3))) ) then
!    if( (Quark1PartType.eq.-Quarks(4)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType) .AND.  &
!        .not.(Quarks(4)%ExtRef.eq.-1 .and. tag_f.eq.1 .and. abs(Quark1PartType).eq.(Quarks(2)%PartType))  ) then
   if( (Quark1PartType.eq.-Quarks(4)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType) .AND.  &
       (Quarks(4)%ExtRef.ne.-1 .or. tag_f.ne.1 .or. abs(Quark1PartType).ne.abs(Quarks(2)%PartType)) &
     ) then
!     (II)
      do n1a=0,NumGlu(1)
      do n3b=0,NumGlu(3)
         n1b = NumGlu(1)-n1a
         n3a = NumGlu(3)-n3b
         rIn =n1a+1
         rOut=NumGlu(1)+NumGlu(2)+n3a


! This prevents color issues with a Z on the quark loop, see RR notes
         if ( tag_Z .eq. 1 .and. (n1b+NumGlu(2)+n3a .eq. NumGlu(0)) .and. Quarks(4)%ExtRef .eq. -1) then
!             print  * , "cycle for tag_Z=1 in cur_f_4f",n1a,n1b
            cycle
         endif




         Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(2:3),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1

         do n4a=0,NumGlu(4)
            n4b = NumGlu(4)-n4a
            ! Fer4
            rIn =NumGlu(1)+NumGlu(2)+n3a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/NumGlu(3)+n4a-n3a,n3b,n4a/) )
            if(n3b.ge.1 .or. n4a.ge.1) then
               PMom2(:) = Quarks(4)%Mom + SumMom(Gluons,rIn,rOut)
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(4)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(4)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(4)%Mass*ubar1(:))*PropFac2
               endif
            endif
            if( Quarks(4)%PartType.lt.0 ) then
               ubar0(:) = vgbq(eps2,ubar1)   !! changed from vqg(ubar1,eps2)       ! re-checked
            else
               ubar0(:) = vgq(eps2,ubar1)   !! changed from vbqg(ubar1,eps2)       ! re-checked
            endif

            PMom1 = Quarks(2)%Mom+Quarks(3)%Mom+Quarks(4)%Mom+SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a)
            if(n1a.ge.1 .or. n4b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(4)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(4)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(4)%Mass*ubar0(:))*PropFac1
               endif
            endif
            TmpQuark(1)%Mom  => PMom1(:)
            TmpQuark(1)%Pol  => ubar0(:)
            TmpQuark(1)%Mass => Quarks(4)%Mass
            TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(4)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n4b/) )
            Res(:) = Res(:) + tmp(:)
         enddo
      enddo
      enddo
   endif
return
END FUNCTION




      recursive function g(e,k) result(res)
      implicit none
      complex(8), intent(in) :: e(:,:),k(:,:)
      complex(8)             :: e1(size(e,dim=1))
      complex(8)             :: res(size(e,dim=1))
      complex(8) :: k1(size(e,dim=1))
      complex(8) :: k2(size(e,dim=1))
      complex(8) :: k3(size(e,dim=1))
      complex(8) :: e2(size(e,dim=1))
      complex(8) :: e3(size(e,dim=1))
      complex(8) :: tmp(size(e,dim=1))
      complex(8) :: k1sq, k2sq, k3sq
      integer :: npart, m, m1

      npart = size(e,dim=2)

      if (npart == 1) then
         res = e(:,1)

      elseif (npart == 2) then
         res = vggg(e(:,1),k(:,1),e(:,2),k(:,2))
        else

           res = (0d0,0d0)

           do m=1,npart-1
              k1=sum(k(:,1:m),dim=2)
              k2=sum(k(:,m+1:npart),dim=2)

              e1 = g(e(:,1:m),k(:,1:m))
              e2 = g(e(:,m+1:npart),k(:,m+1:npart))

              tmp = vggg(e1,k1,e2,k2)

              if (m > 1) then
              k1sq = sc_(k1,k1)
               if (abs(k1sq) > propcut) then
              tmp = -(0d0,1d0)*tmp/k1sq
               else
                  tmp = (0d0,0d0)
               endif
              endif

              if (m + 1 < npart) then
              k2sq = sc_(k2,k2)
                if (abs(k2sq) > propcut) then
              tmp = -(0d0,1d0)*tmp/k2sq
                else
              tmp = (0d0,0d0)
                endif
              endif

              res = res + tmp

      if (m <= npart-2) then

        do m1=m+1,npart-1
!        e1 = g(e(:,1:m),k(:,1:m))  ! e1 is already computed
        e2=g(e(:,m+1:m1),k(:,m+1:m1))
        e3=g(e(:,m1+1:npart),k(:,m1+1:npart))
        k2 = sum(k(:,m+1:m1),dim=2)
        k3 = sum(k(:,m1+1:npart),dim=2)
        tmp = vgggg(e1,e2,e3)
        if (m > 1) then
        k1sq = sc_(k1,k1)
          if(abs(k1sq) > propcut) then
        tmp = -(0d0,1d0)*tmp/k1sq
         else
            tmp = (0d0,0d0)
          endif
        endif
        if (m+1 < m1) then
         k2sq = sc_(k2,k2)
           if(abs(k2sq) > propcut) then
           tmp = -(0d0,1d0)*tmp/k2sq
           else
            tmp = (0d0,0d0)
           endif
        endif
        if (m1+1 < npart) then
           k3sq = sc_(k3,k3)
           if(abs(k3sq) > propcut) then
           tmp = -(0d0,1d0)*tmp/k3sq
           else
           tmp = (0d0,0d0)
           endif
        endif
        res = res + tmp
        enddo

        endif
        enddo
        endif

      end function





FUNCTION cur_g(Gluons,NumGlu)!  note the off-shell gluon has be counted in NumGlu
implicit none
type(PtrToParticle) :: Gluons(1:)
complex(8) :: cur_g(1:Dv)
integer :: Ngluons,NumGlu
complex(8) :: glu_subcur(1:Dv,1:36)  ! max. 8 gluons allowed
complex(8) :: mom_sum(1:Dv,1:36), PropFactor,PropDenom
integer :: ind0,ind1,ind2,ind3,j,l,mu
integer :: a,b,i1,i2

!DEC$ IF (_DebugWriteCurrents==1)
character :: parts(20)*4
integer :: parti(20)

    do i1=1,NumGlu-1
       if( Gluons(i1)%ExtRef.eq.-1 ) then
          exit
       else
          parti(i1)=Gluons(i1)%ExtRef
       endif
       if(i1.eq.NumGlu-1) print *, parti(1:NumGlu-1)
!       write (parts(1:20), '(I20)') parti(1:NumGlu-1)
    enddo
!DEC$ ENDIF

   Ngluons = NumGlu-1

   do a=0,Ngluons-1
      do b=1,Ngluons-a

         i1 = b
         i2 = a+b
         ind0 = linear_map(i1,i2,Ngluons)

         if (i1.eq.i2) then
            glu_subcur(1:Dv,ind0) = Gluons(i1)%Pol(1:Dv)
            mom_sum(1:Dv,ind0)    = Gluons(i1)%Mom(1:Dv)
         else
            mom_sum(1:Dv,ind0) = mom_sum(1:Dv,ind0-Ngluons+i2-i1-1) + Gluons(i2)%Mom(1:Dv)
            if ( i1 .ne. 1 .or. i2 .ne. Ngluons ) then
               PropDenom = mom_sum(1:Dv,ind0).Ndot.mom_sum(1:Dv,ind0)
               if( abs(PropDenom).lt.PropCut ) cycle
               PropFactor = (0d0,-1d0)/PropDenom
            else
               PropFactor = 1d0
            endif
            do mu=1,Dv
               glu_subcur(mu,ind0) = 0d0
            enddo
            do j=i1,i2-1
               ind1 = linear_map(i1,j,Ngluons)
               ind2 = linear_map(j+1,i2,Ngluons)
               glu_subcur(1:Dv,ind0) = glu_subcur(1:Dv,ind0) +            &
                   eval_TripVert( mom_sum(1:Dv,ind1),mom_sum(1:Dv,ind2),  &
                                   glu_subcur(1:Dv,ind1),glu_subcur(1:Dv,ind2)) * PropFactor
            enddo
            do j=i1,i2-2
               do l=j+1,i2-1
                  ind1 = linear_map(i1,j,Ngluons)
                  ind2 = linear_map(j+1,l,Ngluons)
                  ind3 = linear_map(l+1,i2,Ngluons)
                  glu_subcur(1:Dv,ind0) = glu_subcur(1:Dv,ind0) +  &
                      eval_QuadVert(                                       &
                        glu_subcur(1:Dv,ind1),glu_subcur(1:Dv,ind2),glu_subcur(1:Dv,ind3)) * PropFactor
               enddo
            enddo
         endif
      enddo
   enddo
   cur_g(1:Dv) = glu_subcur(1:Dv,ind0)

return
END FUNCTION





      recursive function f(e,k,sp,p,mass,flout,flin,ms) result(res)
      implicit none
      complex(8), intent(in) :: e(:,:), k(:,:)
      complex(8), intent(in) :: sp(:), p(:)
      character,  intent(in) :: flin*3  ! flavor of off-shell f-line
      character,  intent(in)  :: flout*3 ! flavor of on-shell f-line
      integer, intent(in) ::  ms
      integer             :: ms1,m,ng1, ng2
      integer :: ngluon
      complex(8)             :: res(size(sp))
      complex(8)             :: tmp(size(sp))
      complex(8)             :: k1(size(p))
      complex(8)             :: k2(size(p))
      complex(8)             :: sp2(size(sp))
      complex(8)             :: sp3(size(sp))
      complex(8)             :: e1(size(e,dim=1))
      complex(8)             :: e2(size(e,dim=1))
      complex(8)  :: k1sq,k2sq,k3sq
      real(8) :: mass

      ngluon = size(e,dim=2)
      ng1 = ms   !#gluons to the left of a f-line
      ng2 = ngluon - ms  !#gluons to the right of the f-line

      if (flout.ne.flin) then
         res = (0d0,0d0)
      else

      if (ng2 < 0) write(*,*) 'WRONG DEFINITION OF CURRENT A'

      if (ngluon == 0) then
         res = sp

      else

       res = (0d0,0d0)
       do m=0,ng2-1
           k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
           e1=g(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon))
           k1sq=sc_(k1,k1)

           k2 = sum(k(:,1:ng1+m),dim=2)
           k2 = k2 + p
           k2sq = sc_(k2,k2)-mass**2
           sp2 = f(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,flout,flin,ng1)
           if (ng1 >0.or.m>0) sp2 = spb2_(sp2,k2)+mass*sp2
           tmp = vqg(sp2,e1)

! print *, m,e1
! pause
           if (m < ng2-1)  then
              if(abs(k1sq) > propcut) then
                  tmp = -(0d0,1d0)/k1sq*tmp
              else
                  tmp = (0d0,0d0)
              endif
           endif

           if (ng1>0.or.m>0) then
              if (abs(k2sq) > propcut) then
                  tmp =  (0d0,1d0)/k2sq*tmp
              else
                  tmp = (0d0,0d0)
               endif
           endif
           res = res + tmp

        enddo

        do m=1,ng1

           k1 = sum(k(:,1:m),dim=2)
           e1=g(e(:,1:m),k(:,1:m))
           k1sq = sc_(k1,k1)

           k2 = sum(k(:,m+1:ngluon),dim=2)
           k2 = k2 + p
           k2sq = sc_(k2,k2) - mass**2
           ms1 = ng1 - m
           sp2=f(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,flout,flin,ms1)

           if (ng2 > 0.or.m < ng1) sp2 = spb2_(sp2,k2)+mass*sp2
           tmp = vgq(e1,sp2)
           if (m > 1) then
              if (abs(k1sq) > propcut) then
                  tmp=-(0d0,1d0)/k1sq*tmp
              else
                  tmp = (0d0,0d0)
              endif
           endif

           if (ng2 > 0.or. m < ng1) then
            if (abs(k2sq) > propcut) then
              tmp=(0d0,1d0)/k2sq*tmp
            else
              tmp = (0d0,0d0)
            endif
           endif

           res = res + tmp
           enddo

         endif
         endif ! endif for flavor consistency condition

      end function f



      recursive function bf(e,k,sp,p,mass,flout,flin,ms) result(res)
      implicit none
      complex(8), intent(in) :: e(:,:), k(:,:)
      complex(8), intent(in) :: sp(:), p(:)
      integer, intent(in) ::  ms
      character, intent(in) :: flout*3  ! flavor of on-shell f-line
      character, intent(in) :: flin*3   ! flavor of off-shell f-line
      integer             :: ms1,m,ng1, ng2
      integer :: ngluon
      complex(8)             :: res(size(sp))
      complex(8)             :: tmp(size(sp))
      complex(8)             :: k1(size(p))
      complex(8)             :: k2(size(p))
      complex(8)             :: sp2(size(sp))
      complex(8)             :: sp3(size(sp))
      complex(8)             :: e1(size(e,dim=1))
      complex(8)             :: e2(size(e,dim=1))
      complex(8)  :: k1sq,k2sq,k3sq
      real(8) :: mass


       if (flout.ne.flin) then
          res = (0d0,0d0)
       else

       ngluon = size(e,dim=2)
      ng1 = ms   !#gluons to the left of a f-line
      ng2 = ngluon - ms  !#gluons to the right of the f-line

      if (ng2 < 0) write(*,*) 'WRONG DEFINITION OF CURRENT B'
      if (ngluon == 0) then
         res = sp
      else

       res = (0d0,0d0)
       do m=0,ng2-1
           k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
           e1=g(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon))
           k1sq=sc_(k1,k1)

           k2 = sum(k(:,1:ng1+m),dim=2)
           k2 = -k2 - p
           k2sq = sc_(k2,k2)-mass**2
           sp2 = bf(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,mass,flout,flin,ng1)
           if (ng1>0.or.m>0) sp2 = spi2_(k2,sp2)+mass*sp2

            tmp = vbqg(sp2,e1)

           if (m < ng2-1) then
            if (abs(k1sq) > propcut) then
           tmp = -(0d0,1d0)/k1sq*tmp
            else
           tmp = (0d0,0d0)
            endif
            endif
           if (ng1>0.or.m>0) then
            if (abs(k2sq) > propcut) then
           tmp =  (0d0,1d0)/k2sq*tmp
            else
            tmp = (0d0,0d0)
             endif
             endif

           res = res + tmp


        enddo


        do m=1,ng1

           k1 = sum(k(:,1:m),dim=2)
           e1=g(e(:,1:m),k(:,1:m))
           k1sq = sc_(k1,k1)

           k2 = sum(k(:,m+1:ngluon),dim=2)
           k2 = -k2 - p
           k2sq = sc_(k2,k2) - mass**2
           ms1 = ng1 - m
           sp2=bf(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,mass,flout,flin,ms1)

           if (ng2 > 0.or.m < ng1) sp2 = spi2_(k2,sp2)+mass*sp2

           tmp = vgbq(e1,sp2)

           if (m > 1) then
               if (abs(k1sq) > propcut) then
              tmp=-(0d0,1d0)/k1sq*tmp
              else
              tmp = (0d0,0d0)
              endif
           endif

           if (ng2 > 0.or. m < ng1) then
              if (abs(k2sq) > propcut) then
              tmp=(0d0,1d0)/k2sq*tmp
              else
              tmp = (0d0,0d0)
              endif
           endif

           res = res + tmp

           enddo

          endif
          endif ! endif for flavor consisency condition
      end function bf






FUNCTION cur_g_2f(Gluons,Quarks,NumGlu) result(Res)           ! Gluons(:) does not include the OFF-shell gluon,  however NumGlu(0) is the number of all gluons
implicit none
integer :: NumGlu(0:3),i,counter
type(PtrToParticle) :: Gluons(1:),Quarks(1:2)
integer :: rIn,rOut,n1a,n1b,n2a,n2b,n3a,n3b
integer,target :: TmpExtRef
complex(8) :: Res(1:Dv)
complex(8) :: u1(1:Ds),ubar2(1:Ds)
complex(8),target :: Eps1(1:Dv)
complex(8) :: Eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(3)+1)
complex(8) :: PMom1(1:Dv),PMom2(1:Dv),PMom4(1:Dv)
complex(8),target :: PMom3(1:Dv)
complex(8) :: PropFac1,PropFac2,PropFac3,PropFac4
integer :: PartKey,HelKey,CurrKey,Hel_Tmp




   res = (0d0,0d0)
   if( Quarks(1)%PartType .ne. -Quarks(2)%PartType ) return
   do n1a=0,NumGlu(1)
   do n3a=0,NumGlu(3)
   do n2a=0,NumGlu(2)
      n1b=NumGlu(1)-n1a
      n2b=NumGlu(2)-n2a
      n3b=NumGlu(3)-n3a
      ! Fer1
      rIn=n1a+1
      rOut=NumGlu(1)+n2a
      PMom1(:) = Quarks(1)%Mom(:)+ SumMom(Gluons,rIn,rOut)
      u1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,(/n1b+n2a,n1b,n2a/))
      if(n1b.ge.1 .or. n2a.ge.1) then
         PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(1)%Mass2)
         if( abs(sc_(PMom1,PMom1)-Quarks(1)%Mass2).lt.PropCut ) cycle
         if( Quarks(1)%PartType.lt.0 ) then
            u1(:) = (-spi2_(PMom1,u1)+Quarks(1)%Mass*u1(:) )*PropFac1
         else
            u1(:) = (+spb2_(u1,PMom1)+Quarks(1)%Mass*u1(:) )*PropFac1
         endif
      endif

      ! Fer2
      rIn=NumGlu(1)+n2a+1
      rOut=NumGlu(1)+NumGlu(2)+n3a
      PMom2(:) = Quarks(2)%Mom(:)+ SumMom(Gluons,rIn,rOut)
      ubar2(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n2b+n3a,n2b,n3a/))
      if(n2b.ge.1 .or. n3a.ge.1) then
         PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
         if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) cycle
         if( Quarks(2)%PartType.lt.0 ) then
            ubar2(:) = (-spi2_(PMom2,ubar2)+Quarks(2)%Mass*ubar2(:) )*PropFac2
         else
            ubar2(:) = (+spb2_(ubar2,PMom2)+Quarks(2)%Mass*ubar2(:) )*PropFac2
         endif
      endif

      if( Quarks(1)%PartType.lt.0 ) then
         Eps1(:)= -vbqq(Dv,ubar2,u1)       ! re-checked
      else
         Eps1(:)= +vbqq(Dv,u1,ubar2)       ! re-checked
      endif
      PMom3(:) = Quarks(1)%Mom(:)+Quarks(2)%Mom(:) + SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+n3a)
      counter=1
      rIn =1
      rOut=n1a
      do i=rIn,rOut
         call CopyParticlePtr(Gluons(i),TmpGluons(counter))
         counter=counter+1
      enddo
      TmpGluons(counter)%Mom => PMom3(:)
      TmpGluons(counter)%Pol => Eps1(:)
      TmpExtRef = -1
      TmpGluons(counter)%ExtRef => TmpExtRef
      counter=counter+1
      rIn =NumGlu(1)+NumGlu(2)+n3a+1
      rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)
      do i=rIn,rOut
         call CopyParticlePtr(Gluons(i),TmpGluons(counter))
         counter=counter+1
      enddo
      Eps2(:) = cur_g(TmpGluons(1:counter-1),1+n1a+n3b+1)


      if(n1a.ge.1 .or. n3b.ge.1) then
         PropFac3 = (0d0,-1d0)/sc_(PMom3,PMom3)
         if( abs(sc_(PMom3,PMom3)).lt.PropCut ) cycle
         Eps2(:) = Eps2(:)*PropFac3
      endif

      Res(:) = Res(:) + Eps2(:)
   enddo
   enddo
   enddo



return
END FUNCTION





FUNCTION cur_g_2fV(Gluons,Quarks,Boson,NumGlu) result(Res)           ! Gluons(:) does not include the OFF-shell gluon,  however NumGlu(0) is the number of all gluons
implicit none
integer :: NumGlu(0:3),i,counter
type(PtrToParticle) :: Gluons(1:),Quarks(1:2),Boson
integer :: rIn,rOut,n1a,n1b,n2a,n2b,n3a,n3b
integer,target :: TmpExtRef
complex(8) :: Res(1:Dv)
complex(8) :: u1(1:Ds),ubar2(1:Ds)
complex(8),target :: Eps1(1:Dv)
complex(8) :: Eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(3)+1)
complex(8) :: PMom1(1:Dv),PMom2(1:Dv),PMom4(1:Dv)
complex(8),target :: PMom3(1:Dv)
complex(8) :: PropFac1,PropFac2,PropFac3,PropFac4
integer :: PartKey,HelKey,CurrKey,Hel_Tmp




!DEC$ IF (_DebugCheckMyImpl1==1)
   if(Quarks(1)%PartType*Quarks(2)%PartType.ge.0) print *,"Error in cur_g_2f: wrong PartTypes"
!DEC$ ENDIF
!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-1-NumGlu(1)-NumGlu(2)-NumGlu(3).ne.0 ) print *, "wrong number of gluons in cur_g_2f"
!DEC$ ENDIF

   res = (0d0,0d0)
   if( Quarks(1)%PartType .ne. -Quarks(2)%PartType ) return
   do n1a=0,NumGlu(1)
   do n3a=0,NumGlu(3)
   do n2a=0,NumGlu(2)
      n1b=NumGlu(1)-n1a
      n2b=NumGlu(2)-n2a
      n3b=NumGlu(3)-n3a

      ! Fer1 and V coupling
      rIn=n1a+1
      rOut=NumGlu(1)+n2a
      PMom1(:) = Quarks(1)%Mom(:) + SumMom(Gluons,rIn,rOut) + Boson%Mom(:)
      u1(:) = cur_f_2fV(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,Boson,(/n1b+n2a,n1b,n2a/))
      !if(n1b.ge.1 .or. n2a.ge.1) then
         PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(1)%Mass2)
         if( abs(sc_(PMom1,PMom1)-Quarks(1)%Mass2).lt.PropCut ) then
            PropFac1=(0d0,0d0)
         endif
         if( Quarks(1)%PartType.lt.0 ) then
            u1(:) = (-spi2_(PMom1,u1)+Quarks(1)%Mass*u1(:) )*PropFac1
         else
            u1(:) = (+spb2_(u1,PMom1)+Quarks(1)%Mass*u1(:) )*PropFac1
         endif
      !endif

      ! Fer2
      rIn=NumGlu(1)+n2a+1
      rOut=NumGlu(1)+NumGlu(2)+n3a
      PMom2(:) = Quarks(2)%Mom(:)+ SumMom(Gluons,rIn,rOut)
      ubar2(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n2b+n3a,n2b,n3a/))
      if(n2b.ge.1 .or. n3a.ge.1) then
         PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
         if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) then
            PropFac2=(0d0,0d0)
         endif
         if( Quarks(2)%PartType.lt.0 ) then
            ubar2(:) = (-spi2_(PMom2,ubar2)+Quarks(2)%Mass*ubar2(:) )*PropFac2
         else
            ubar2(:) = (+spb2_(ubar2,PMom2)+Quarks(2)%Mass*ubar2(:) )*PropFac2
         endif
      endif


      if( Quarks(1)%PartType.lt.0 ) then
         Eps1(:)= -vbqq(Dv,ubar2,u1)
      else
         Eps1(:)= +vbqq(Dv,u1,ubar2)
      endif





      ! Fer1
      rIn=n1a+1
      rOut=NumGlu(1)+n2a
      PMom1(:) = Quarks(1)%Mom(:) + SumMom(Gluons,rIn,rOut)
      u1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,(/n1b+n2a,n1b,n2a/))
      if(n1b.ge.1 .or. n2a.ge.1) then
         PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(1)%Mass2)
         if( abs(sc_(PMom1,PMom1)-Quarks(1)%Mass2).lt.PropCut ) cycle
         if( Quarks(1)%PartType.lt.0 ) then
            u1(:) = (-spi2_(PMom1,u1)+Quarks(1)%Mass*u1(:) )*PropFac1
         else
            u1(:) = (+spb2_(u1,PMom1)+Quarks(1)%Mass*u1(:) )*PropFac1
         endif
      endif

      ! Fer2 and V coupling
      rIn=NumGlu(1)+n2a+1
      rOut=NumGlu(1)+NumGlu(2)+n3a
      PMom2(:) = Quarks(2)%Mom(:)+ SumMom(Gluons,rIn,rOut) + Boson%Mom(:)
      ubar2(:) = cur_f_2fV(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,Boson,(/n2b+n3a,n2b,n3a/))
      !if(n2b.ge.1 .or. n3a.ge.1) then
         PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
         if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) cycle
         if( Quarks(2)%PartType.lt.0 ) then
            ubar2(:) = (-spi2_(PMom2,ubar2)+Quarks(2)%Mass*ubar2(:) )*PropFac2
         else
            ubar2(:) = (+spb2_(ubar2,PMom2)+Quarks(2)%Mass*ubar2(:) )*PropFac2
         endif
      !endif

      if( Quarks(1)%PartType.lt.0 ) then
         Eps1(:)= Eps1(:) - vbqq(Dv,ubar2,u1)
      else
         Eps1(:)= Eps1(:) + vbqq(Dv,u1,ubar2)
      endif
!

      PMom3(:) = Quarks(1)%Mom(:)+Quarks(2)%Mom(:) + SumMom(Gluons,n1a+1,NumGlu(1)+NumGlu(2)+n3a) + Boson%Mom(:)
      counter=1
      rIn =1
      rOut=n1a
      do i=rIn,rOut
         call CopyParticlePtr(Gluons(i),TmpGluons(counter))
         counter=counter+1
      enddo
      TmpGluons(counter)%Mom => PMom3(:)
      TmpGluons(counter)%Pol => Eps1(:)
      TmpExtRef = -1
      TmpGluons(counter)%ExtRef => TmpExtRef
      counter=counter+1
      rIn =NumGlu(1)+NumGlu(2)+n3a+1
      rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)
      do i=rIn,rOut
         call CopyParticlePtr(Gluons(i),TmpGluons(counter))
         counter=counter+1
      enddo
      Eps2(:) = cur_g(TmpGluons(1:counter-1),1+n1a+n3b+1)


      if(n1a.ge.1 .or. n3b.ge.1) then
         PropFac3 = (0d0,-1d0)/sc_(PMom3,PMom3)
         if( abs(sc_(PMom3,PMom3)).lt.PropCut ) cycle
         Eps2(:) = Eps2(:)*PropFac3
      endif

      Res(:) = Res(:) + Eps2(:)

   enddo
   enddo
   enddo



return
END FUNCTION








FUNCTION cur_f_2f(Gluons,Quarks,Quark1PartType,NumGlu) result(Res)           ! Quarks(:) does not include the OFF-shell quark
implicit none
complex(8) :: Res(1:Ds)
integer :: NumGlu(0:2),i,rIn,rOut,Quark1PartType
type(PtrToParticle) :: Gluons(1:),Quarks(2:2)      ! off-shell quark is not included in Quarks(:)
complex(8) :: GluMom(1:Dv,NumGlu(0)), QuarkMom(1:Dv)
complex(8) :: GluPol(1:Dv,NumGlu(0)), QuarkPol(1:Ds)
character :: FerFla1*3,FerFla2*3
integer :: PartKey,HelKey,CurrKey,Hel_Tmp






!DEC$ IF (_DebugGeneralChecks==1)
   if( Quarks(2)%PartType .eq.0 .or. .not.IsAQuark(Quarks(2)%PartType)) then
      print *, "Error in cur_f_2f"
      stop
   endif
!DEC$ ENDIF
!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2).ne.0 ) print *, "wrong number of gluons in cur_f_2f"
    if( Quarks(2)%PartType.ne.-Quark1PartType ) print *, "unequal flavors in cur_f_2f"
!DEC$ ENDIF


   do i=1,NumGlu(0)
    GluMom(1:Dv,i) = Gluons(i)%Mom(1:Dv)
    GluPol(1:Dv,i) = Gluons(i)%Pol(1:Dv)
   enddo
   QuarkMom(1:Dv) = Quarks(2)%Mom(1:Dv)
   QuarkPol(1:Ds) = Quarks(2)%Pol(1:Ds)

   if( abs(Quark1PartType).eq.5 ) then
    FerFla1="top"
   elseif( abs(Quark1PartType).eq.6 ) then
    FerFla1="bot"
   elseif( abs(Quark1PartType).eq.3 ) then
    FerFla1="chm"
   else
    FerFla1="str"
   endif

   if( abs(Quarks(2)%PartType).eq.5 ) then
    FerFla2="top"
   elseif( abs(Quarks(2)%PartType).eq.6 ) then
    FerFla2="bot"
   elseif( abs(Quarks(2)%PartType).eq.3 ) then
    FerFla2="chm"
   else
    FerFla2="str"
   endif

   rIn =1
   rOut=NumGlu(0)
   if( Quarks(2)%PartType .gt.0 ) then      !    X----->----
      Res(:) = f(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),Quarks(2)%Mass,FerFla2,FerFla1,NumGlu(1))
   else                                     !    X-----<----
      Res(:) = bf(GluPol(1:Dv,rIn:rOut),GluMom(1:Dv,rIn:rOut),QuarkPol(1:Ds),QuarkMom(1:Dv),Quarks(2)%Mass,FerFla2,FerFla1,NumGlu(1))
   endif


return
END FUNCTION






FUNCTION cur_g_4f(Gluons,Quarks,NumGlu) result(res)           ! Gluons(:) does not include the OFF-shell gluon,  however NumGlu is the number of all gluons
implicit none
integer,intent(in) :: NumGlu(0:5)
type(PtrToParticle) :: Gluons(1:),Quarks(1:)
integer :: na,nb,nc,nd,ne,nf,ng,nh,ni,nj,nk
integer :: rIn,rOut
integer :: tag_f,counter,i
complex(8) :: res(Dv)
type(PtrToParticle) :: TmpGluons(1:2+NumGlu(1)+NumGlu(3)+NumGlu(5))
complex(8),target :: TmpMom1(1:Dv),TmpMom2(1:Dv)
integer,target :: TmpExtRef1,TmpExtRef2
complex(8),target :: Eps1(1:Dv)
complex(8),target :: Eps2(1:Dv)
complex(8) :: Eps3(1:Dv)
complex(8) :: u1(1:Ds)
complex(8) :: ubar2(1:Ds)
complex(8) :: PropFac1,PropFac2,PropFac3,PropFac4
complex(8) :: PMom1(1:Dv)
complex(8) :: PMom2(1:Dv)
complex(8) :: PMom3(1:Dv)
complex(8) :: PMom4(1:Dv)

!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-1-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5).ne.0 ) print *, "wrong number of gluons in cur_g_4f"
!DEC$ ENDIF

      res = (0d0,0d0)
      if( (Quarks(1)%PartType.eq.-Quarks(2)%PartType .and. Quarks(3)%PartType.eq.-Quarks(4)%PartType) &
     .OR. (Quarks(1)%PartType.eq.-Quarks(4)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType) ) then
!        type (1)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            nf=NumGlu(5)-ne
            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(2:4),Quarks(1)%PartType,(/nd+NumGlu(3)+NumGlu(4)+ne,nd,NumGlu(3),NumGlu(4),ne/),0,0)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(1)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Quarks(1)%Mass2).lt.PropCut ) cycle
            if( Quarks(1)%PartType.lt.0 ) then
               u1 = ( spb2_(u1,PMom2) + Quarks(1)%Mass*u1 )*PropFac2
            else
               u1 = (-spi2_(PMom2,u1) + Quarks(1)%Mass*u1 )*PropFac2
            endif

            rIn = na+1
            rOut= NumGlu(1)+nc
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,(/nb+nc,nb,nc/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom
            if( nb.ge.1 .or. nc.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(1)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Quarks(1)%Mass2).lt.PropCut ) cycle
               if( Quarks(1)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom3,ubar2) + Quarks(1)%Mass*ubar2 )*PropFac3
               else
                  ubar2 = (+spb2_(ubar2,PMom3) + Quarks(1)%Mass*ubar2 )*PropFac3
               endif
            endif

            if( Quarks(1)%PartType.lt.0 ) then
                Eps1 = -vbqq(Dv,u1,ubar2)       ! re-checked
            else
                Eps1 = +vbqq(Dv,ubar2,u1)       ! re-checked
            endif

            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2
         enddo
         enddo
         enddo


!        type (2)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(4)
         do ne=0,NumGlu(5)        ! can be replaced by above ne-loop
            nb=NumGlu(1)-na
            nd=NumGlu(4)-nc
            nf=NumGlu(5)-ne

            rIn = na+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nc
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(1:3),Quarks(4)%PartType,(/nb+NumGlu(2)+NumGlu(3)+nc,nb,NumGlu(2),NumGlu(3),nc/),0,0)
            PMom2 =  SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom + Quarks(2)%Mom + Quarks(3)%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
            if( Quarks(4)%PartType.lt.0 ) then
              u1 = (+spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
            else
              u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
            endif
            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/nd+ne,nd,ne/))
            PMom3  = SumMom(Gluons,rIn,rOut) + Quarks(4)%Mom
            if( nd.ge.1 .or. ne.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(4)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom3,ubar2) + Quarks(4)%Mass*ubar2 )*PropFac3
               else
                  ubar2 = (+spb2_(ubar2,PMom3) + Quarks(4)%Mass*ubar2 )*PropFac3
               endif
            endif

            if( Quarks(4)%PartType.lt.0 ) then
              Eps1 = +vbqq(Dv,u1,ubar2)       ! re-checked
            else
              Eps1 = -vbqq(Dv,ubar2,u1)       ! re-checked
            endif

            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2
         enddo
         enddo
         enddo
      endif



      if( Quarks(1)%PartType.eq.-Quarks(2)%PartType .and. Quarks(3)%PartType.eq.-Quarks(4)%PartType) then
!        type(3)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(3)
         do nf=0,NumGlu(3)-ne
         do nh=0,NumGlu(4)   ! this loop can be placed after Eps1 has been calculated
         do nj=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            ng=NumGlu(3)-ne-nf
            ni=NumGlu(4)-nh
            nk=NumGlu(5)-nj

            rIn = na+1
            rOut= NumGlu(1)+nc
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,(/nb+nc,nb,nc/))
            PMom1 = SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom
            if( nb.ge.1 .or. nc.ge.1 ) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1) - Quarks(1)%Mass2)
               if( abs(sc_(PMom1,PMom1) - Quarks(1)%Mass2).lt.PropCut ) cycle
               if( Quarks(1)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom1,ubar2) + Quarks(1)%Mass*ubar2 )*PropFac1
               else
                  ubar2 = (+spb2_(ubar2,PMom1) + Quarks(1)%Mass*ubar2 )*PropFac1
               endif
            endif
            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+ne
            u1 = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/nd+ne,nd,ne/))
            PMom2 = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom
            if( nd.ge.1 .or. ne.ge.1 ) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
               if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  u1 = (-spi2_(PMom2,u1) + Quarks(2)%Mass*u1 )*PropFac2
               else
                  u1 = (+spb2_(u1,PMom2) + Quarks(2)%Mass*u1 )*PropFac2
               endif
            endif

            if( Quarks(2)%PartType.lt.0 ) then
              Eps1 = +vbqq(Dv,ubar2,u1)       ! re-checked
            else
              Eps1 = -vbqq(Dv,u1,ubar2)       ! re-checked
            endif
            TmpMom1 = PMom1 + PMom2
            PropFac3 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
            if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
            Eps1 = Eps1*PropFac3


            rIn = NumGlu(1)+NumGlu(2)+ne+nf+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nh
            u1 = cur_f_2f(Gluons(rIn:rOut),Quarks(3:3),-Quarks(3)%PartType,(/ng+nh,ng,nh/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom
            if( ng.ge.1 .or. nh.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(3)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Quarks(3)%Mass2).lt.PropCut ) cycle
               if( Quarks(3)%PartType.lt.0 ) then
                  u1 = (-spi2_(PMom3,u1) + Quarks(3)%Mass*u1 )*PropFac3
               else
                  u1 = (+spb2_(u1,PMom3) + Quarks(3)%Mass*u1 )*PropFac3
               endif
            endif
            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nh+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/ni+nj,ni,nj/))
            PMom4 = SumMom(Gluons,rIn,rOut) + Quarks(4)%Mom
            if( ni.ge.1 .or. nj.ge.1 ) then
               PropFac4 = (0d0,1d0)/(sc_(PMom4,PMom4) - Quarks(4)%Mass2)
               if( abs(sc_(PMom4,PMom4) - Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom4,ubar2) + Quarks(4)%Mass*ubar2 )*PropFac4
               else
                  ubar2 = (+spb2_(ubar2,PMom4) + Quarks(4)%Mass*ubar2 )*PropFac4
               endif
            endif

            if( Quarks(4)%PartType.lt.0 ) then
               Eps2 = +vbqq(Dv,u1,ubar2)       ! re-checked
            else
               Eps2 = -vbqq(Dv,ubar2,u1)       ! re-checked
            endif
            TmpMom2 = PMom3 + PMom4
            PropFac1 = (0d0,-1d0)/sc_(TmpMom2,TmpMom2)
            if( abs(sc_(TmpMom2,TmpMom2)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1


            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+ne+1
            rOut=NumGlu(1)+NumGlu(2)+ne+nf
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpGluons(counter)%Mom => TmpMom2(:)
            TmpGluons(counter)%Pol => Eps2(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps3(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+nk+2)
            Res = Res + Eps3
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
      endif


return
END FUNCTION




FUNCTION cur_f_6f(Gluons,Quarks,Quark1PartType,NumGlu,tag_f,tag_Z) result(res)           ! Quarks(:) does not include the OFF-shell quark
implicit none
integer :: NumGlu(0:6),Quark1PartType,tag_f,tag_Z
type(PtrToParticle) :: Gluons(1:),Quarks(2:6)
integer,target :: TmpPartType,TmpExtRef
complex(8) :: Res(1:Ds),tmp(1:Ds)
! complex(8) :: Res1(1:Ds),Res2(1:Ds),Res3(1:Ds),Res4(1:Ds)
complex(8) :: u1(1:Ds),ubar1(1:Ds)
complex(8),target :: ubar0(1:Ds)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(6)),TmpQuark(1:1)
complex(8) :: PropFac1,PropFac2
complex(8),target :: PMom1(1:Dv)
complex(8),target :: PMom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,n5a,n5b,n6a,n6b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5)-NumGlu(6).ne.0 ) print *, "wrong number of gluons in cur_f_6f"
!DEC$ ENDIF

    Res = (0d0,0d0)
!     Res1=(0d0,0d0); Res2=(0d0,0d0); Res3=(0d0,0d0); Res4=(0d0,0d0);


!   (A)
    if( Quark1PartType.eq.-Quarks(2)%PartType .AND. (Quarks(3)%PartType.eq.-Quarks(4)%PartType .or. Quarks(3)%PartType.eq.-Quarks(6)%PartType) &
        .AND. (Quarks(2)%ExtRef.ne.-1 .or. tag_f.ne.1) &
      ) then

      do n2a=0,NumGlu(2)
      do n6a=0,NumGlu(6)
         n2b = NumGlu(2)-n2a
         n6b = NumGlu(6)-n6a

         rIn =NumGlu(1)+n2a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
         Eps2 = cur_g_4f(Gluons(rIn:rOut),Quarks(3:6),(/1+n2b+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a,n2b,NumGlu(3),NumGlu(4),NumGlu(5),n6a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut) cycle
         Eps2 = Eps2*PropFac1
         do n1a=0,NumGlu(1)
            n1b = NumGlu(1)-n1a
            ! Fer2
            rIn =n1a+1
            rOut=NumGlu(1)+n2a
            ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n1b+n2a,n1b,n2a/) )
            if(n1b.ge.1 .or. n2a.ge.1) then
               PMom2(:) = Quarks(2)%Mom + SumMom(Gluons,rIn,rOut)
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(2)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(2)%Mass*ubar1(:))*PropFac2
               endif
            endif
            if( Quarks(2)%PartType.lt.0 ) then
               ubar0(:) = vbqg(ubar1,eps2)       ! re-checked
            else
               ubar0(:) = vqg(ubar1,eps2)       ! re-checked
            endif

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            PMom1 = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom  ! can be simplified with PMom1(:)
            if(n1a.ge.1 .or. n6b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(2)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(2)%Mass*ubar0(:))*PropFac1
               endif
            endif

            TmpQuark(1)%Mom  => PMom1(:)
            TmpQuark(1)%Pol  => ubar0(:)
            TmpQuark(1)%Mass => Quarks(2)%Mass
            TmpQuark(1)%Mass2=> Quarks(2)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(2)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n6b/) )

            Res(:) = Res(:) + tmp(:)
!             Res1(:) = Res1(:) + tmp(:)
!             print *, "1",tmp(:)
         enddo
      enddo
      enddo
    endif



!   (B)
    if( Quark1PartType.eq.-Quarks(6)%PartType .AND. (Quarks(2)%PartType.eq.-Quarks(5)%PartType .or. Quarks(2)%PartType.eq.-Quarks(3)%PartType) &
        .AND. (Quarks(6)%ExtRef.ne.-1 .or. tag_f.ne.1) &
      ) then
      do n1a=0,NumGlu(1)
      do n5b=0,NumGlu(5)
         n1b = NumGlu(1)-n1a
         n5a = NumGlu(5)-n5b

         rIn =n1a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a

! preventing color issues with Z on the loop, see RR notes 05-22-2013
!note: only tested for zero gluons
         if ( tag_Z .eq. 1 .and. (n1b+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a .eq. NumGlu(0)) ) then
!             print  * , "cycle for tag_Z=1 in cur_f_6f"
            cycle
         endif

         Eps2 = cur_g_4f(Gluons(rIn:rOut),Quarks(2:5),(/1+n1b+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a,n1b,NumGlu(2),NumGlu(3),NumGlu(4),n5a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Quarks(5)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1

         do n6a=0,NumGlu(6)
            n6b = NumGlu(6)-n6a
            ! Fer6
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+n5a+1
            rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(6:6),-Quarks(6)%PartType,(/n5b+n6a,n5b,n6a/) )
            if(n5b.ge.1 .or. n6a.ge.1) then
               PMom2(:) = SumMom(Gluons,rIn,rOut) + Quarks(6)%Mom
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(6)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(6)%Mass2).lt.PropCut ) cycle
               if( Quarks(6)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(6)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(6)%Mass*ubar1(:))*PropFac2
               endif
            endif
            if( Quarks(6)%PartType.lt.0 ) then
               ubar0(:) = vgbq(Eps2,ubar1)       ! re-checked
            else
               ubar0(:) = vgq(Eps2,ubar1)       ! re-checked
            endif

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            PMom1 = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom
            if(n1a.ge.1 .or. n6b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(6)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(6)%Mass2).lt.PropCut ) cycle
               if( Quarks(6)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(6)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(6)%Mass*ubar0(:))*PropFac1
               endif
            endif
            TmpQuark(1)%Mom => PMom1(:)
            TmpQuark(1)%Pol => ubar0(:)
            TmpQuark(1)%Mass => Quarks(6)%Mass
            TmpQuark(1)%Mass2=> Quarks(6)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(6)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),Quark1PartType,(/counter-1,n1a,n6b/) )

            Res(:) = Res(:) + tmp(:)

!             Res2(:) = Res2(:) + tmp(:)

!             print *, "2",tmp(:)
         enddo
      enddo
      enddo
    endif




!   (C)
!     if( Quarks(5)%PartType.eq.-Quarks(6)%PartType .AND. (Quark1PartType.eq.-Quarks(2)%PartType .or. Quark1PartType.eq.-Quarks(4)%PartType) &
!         .AND. .not.(Quarks(6)%ExtRef.eq.-1 .and. tag_f.eq.1) &
!       ) then
!     if( Quark1PartType.eq.-Quarks(2)%PartType .and. Quarks(5)%PartType.eq.-Quarks(6)%PartType .and..not.(Quarks(2)%ExtRef.eq.-1 .and. tag_f.eq.1) &
!    .OR. Quark1PartType.eq.-Quarks(4)%PartType .and. Quarks(5)%PartType.eq.-Quarks(6)%PartType .and..not.(Quarks(4)%ExtRef.eq.-1 .and. tag_f.eq.1) &
!       ) then

    if( Quarks(5)%PartType.eq.-Quarks(6)%PartType .AND. &
        ((Quark1PartType.eq.-Quarks(2)%PartType .and. (Quarks(2)%ExtRef.ne.-1.or.tag_f.ne.1) ) &
    .OR. (Quark1PartType.eq.-Quarks(4)%PartType .and. (Quarks(4)%ExtRef.ne.-1.or.tag_f.ne.1) ))&
      ) then

      do n1a=0,NumGlu(1)
      do n4a=0,NumGlu(4)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n4b = NumGlu(4)-n4a
         n6b = NumGlu(6)-n6a

            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(5:6),(/1+n4b+NumGlu(5)+n6a,n4b,NumGlu(5),n6a/))
            PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
            PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
            if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(2:4),Quark1PartType,(/n1b+NumGlu(2)+NumGlu(3)+n4a,n1b,NumGlu(2),NumGlu(3),n4a/),0,0)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom

            if( Quark1PartType.eq.-Quarks(2)%PartType) then
                if( Quarks(2)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(2)%Mass*u1 )*PropFac2
                  ubar0 = vbqg(u1,eps2)       ! re-checked
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom   ! this PMom2 will be re-used below
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(2)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(2)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(2)%Mass*u1 )*PropFac2
                  ubar0 = vqg(u1,eps2)       ! re-checked
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(2)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(2)%Mass
                TmpQuark(1)%Mass2=> Quarks(2)%Mass2
            elseif( Quark1PartType.eq.-Quarks(4)%PartType) then
                if( Quarks(4)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vbqg(u1,eps2)       ! re-checked
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(4)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vqg(u1,eps2)       ! re-checked
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(4)%Mass
                TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            endif
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpPartType = -Quark1PartType
            TmpQuark(1)%PartType => TmpPartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),Quark1PartType,(/counter-1,n1a,n6b/) )

            Res(:) = Res(:) + tmp(:)
!             Res3(:) = Res3(:) + tmp(:)
!             print *, "3",tmp(:)
      enddo
      enddo
      enddo
    endif


!   (D)
    if( Quarks(2)%PartType.eq.-Quarks(3)%PartType .AND. ( &
        (Quark1PartType.eq.-Quarks(4)%PartType .and. (Quarks(4)%ExtRef.ne.-1.or.tag_f.ne.1)) &
   .OR. (Quark1PartType.eq.-Quarks(6)%PartType .and. (Quarks(6)%ExtRef.ne.-1.or.tag_f.ne.1))  )) then
!     if( Quark1PartType.eq.-Quarks(4)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType .and..not.(Quarks(4)%ExtRef.eq.-1 .and. tag_f.eq.1) &
!     .OR.Quark1PartType.eq.-Quarks(6)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType .and..not.(Quarks(6)%ExtRef.eq.-1 .and. tag_f.eq.1) &
!       ) then

      do n1a=0,NumGlu(1)
      do n3a=0,NumGlu(3)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n3b = NumGlu(3)-n3a
         n6b = NumGlu(6)-n6a

            counter=1
            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+n3a


            Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(2:3),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
            PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
            PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
            if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1

            rIn = NumGlu(1)+NumGlu(2)+n3a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(4:6),Quark1PartType,(/n3b+NumGlu(4)+NumGlu(5)+n6a,n3b,NumGlu(4),NumGlu(5),n6a/),tag_f,0)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom

            if( Quark1PartType.eq.-Quarks(4)%PartType ) then
                if( Quarks(4)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgbq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom   !this PMom2 will be re-used below  ! CAN be written as PMom2=PMom2+PMom1
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(4)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(4)%Mass
                TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            elseif(Quark1PartType.eq.-Quarks(6)%PartType) then
                if( Quarks(6)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(6)%Mass*u1 )*PropFac2
                  ubar0 = vgbq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(6)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(6)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(6)%Mass*u1 )*PropFac2
                  ubar0 = vgq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(6)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(6)%Mass
                TmpQuark(1)%Mass2=> Quarks(6)%Mass2
            endif
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpPartType = -Quark1PartType
            TmpQuark(1)%PartType => TmpPartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),Quark1PartType,(/counter-1,n1a,n6b/))

            Res(:)  = Res(:) + Tmp(:)
!             Res4(:) = Res4(:) + Tmp(:)
!             print *, "4",tmp(:)
      enddo
      enddo
      enddo
    endif
!    print *, 'res', res
!    pause
!     Res(:) = Res1(:)+Res2(:)+Res3(:)+Res4(:)
!     print *, "Res1:",Res1(1:4)
!     print *, "Res2:",Res2(1:4)
!     print *, "Res3:",Res3(1:4)
!     print *, "Res4:",Res4(1:4)
return
END FUNCTION








FUNCTION cur_g_4fV(Gluons,Quarks,Boson,BosonVertex,NumGlu) result(res)           ! Gluons(:) does not include the OFF-shell gluon,  however NumGlu is the number of all gluons
implicit none
integer,intent(in) :: NumGlu(0:5),BosonVertex
type(PtrToParticle) :: Gluons(1:),Quarks(1:),Boson
integer :: na,nb,nc,nd,ne,nf,ng,nh,ni,nj,nk,BosonVertex_mod
integer :: rIn,rOut
integer :: tag_f,counter,i
complex(8) :: res(Dv)
type(PtrToParticle) :: TmpGluons(1:2+NumGlu(1)+NumGlu(3)+NumGlu(5))
complex(8),target :: TmpMom1(1:Dv),TmpMom2(1:Dv)
integer,target :: TmpExtRef1,TmpExtRef2
complex(8),target :: Eps1(1:Dv)
complex(8),target :: Eps2(1:Dv)
complex(8) :: Eps3(1:Dv)
complex(8) :: u1(1:Ds)
complex(8) :: ubar2(1:Ds)
complex(8) :: PropFac1,PropFac2,PropFac3,PropFac4
complex(8) :: PMom1(1:Dv)
complex(8) :: PMom2(1:Dv)
complex(8) :: PMom3(1:Dv)
complex(8) :: PMom4(1:Dv)

!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-1-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5).ne.0 ) print *, "wrong number of gluons in cur_g_4f"
!DEC$ ENDIF

    res = (0d0,0d0)
    if (Quarks(1)%PartType.eq.-Quarks(4)%PartType .and. Quarks(2)%PartType.eq.-Quarks(3)%PartType) then
       ! for this flavor structure, the BosonVertex implies the following:
       ! BV = 1 : Boson on quark 1 only
       ! BV = 2 : Boson on both quarks 2 and 3
       ! BV = 3 : Boson on quark 4 only
       ! See RR notes

       if( BosonVertex.eq.1)  then
          !        type (1)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            nf=NumGlu(5)-ne
            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(2:4),Quarks(1)%PartType,(/nd+NumGlu(3)+NumGlu(4)+ne,nd,NumGlu(3),NumGlu(4),ne/),0,0)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(1)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Quarks(1)%Mass2).lt.PropCut ) cycle
            if( Quarks(1)%PartType.lt.0 ) then
               u1 = ( spb2_(u1,PMom2) + Quarks(1)%Mass*u1 )*PropFac2
            else
               u1 = (-spi2_(PMom2,u1) + Quarks(1)%Mass*u1 )*PropFac2
            endif

            rIn = na+1
            rOut= NumGlu(1)+nc
            ubar2 = cur_f_2fV(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,Boson,(/nb+nc,nb,nc/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom + Boson%Mom
            PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(1)%Mass2)
            if( abs(sc_(PMom3,PMom3) - Quarks(1)%Mass2).lt.PropCut ) cycle
            if( Quarks(1)%PartType.lt.0 ) then
               ubar2 = (-spi2_(PMom3,ubar2) + Quarks(1)%Mass*ubar2 )*PropFac3
            else
               ubar2 = (+spb2_(ubar2,PMom3) + Quarks(1)%Mass*ubar2 )*PropFac3
            endif

            if( Quarks(1)%PartType.lt.0 ) then
                Eps1 = -vbqq(Dv,u1,ubar2)       ! re-checked
            else
                Eps1 = +vbqq(Dv,ubar2,u1)       ! re-checked
            endif

            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2
         enddo
         enddo
         enddo
      endif! Bosonvertex

      if( BosonVertex.eq.1 .or. BosonVertex.eq.2 .or. BosonVertex .eq. 3)  then
!        type (2)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(4)
         do ne=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(4)-nc
            nf=NumGlu(5)-ne

            rIn = na+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nc
            BosonVertex_mod = BosonVertex  + 1
!            if (BosonVertex .eq. 1) BosonVertex_mod=BosonVertex
            u1 = cur_f_4fV(Gluons(rIn:rOut),Quarks(1:3),Quarks(4)%PartType,Boson,BosonVertex_mod,(/nb+NumGlu(2)+NumGlu(3)+nc,nb,NumGlu(2),NumGlu(3),nc/),0,0)
            PMom2 =  SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom + Quarks(2)%Mom + Quarks(3)%Mom + Boson%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
            if( Quarks(4)%PartType.lt.0 ) then
              u1 = (+spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
            else
              u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
            endif

            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/nd+ne,nd,ne/))
            PMom3  = SumMom(Gluons,rIn,rOut) + Quarks(4)%Mom
            if( nd.ge.1 .or. ne.ge.1 ) then
               PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(4)%Mass2)
               if( abs(sc_(PMom3,PMom3) - Quarks(4)%Mass2).lt.PropCut ) cycle
               if( Quarks(4)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom3,ubar2) + Quarks(4)%Mass*ubar2 )*PropFac3
               else
                  ubar2 = (+spb2_(ubar2,PMom3) + Quarks(4)%Mass*ubar2 )*PropFac3
               endif
            endif

            if( Quarks(4)%PartType.lt.0 ) then
              Eps1 = +vbqq(Dv,u1,ubar2)
            else
              Eps1 = -vbqq(Dv,ubar2,u1)
            endif

            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2

         enddo
         enddo
         enddo

      endif ! BosonVertex


      if (BosonVertex .eq. 1 .or. BosonVertex .eq. 2 .or. BosonVertex .eq. 3 ) then

         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            nf=NumGlu(5)-ne
            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            BosonVertex_mod=BosonVertex
            u1 = cur_f_4fV(Gluons(rIn:rOut),Quarks(2:4),Quarks(1)%PartType,Boson,BosonVertex_mod,(/nd+NumGlu(3)+NumGlu(4)+ne,nd,NumGlu(3),NumGlu(4),ne/),0,0)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Boson%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(1)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Quarks(1)%Mass2).lt.PropCut ) cycle
            if( Quarks(1)%PartType.lt.0 ) then
               u1 = ( spb2_(u1,PMom2) + Quarks(1)%Mass*u1 )*PropFac2
            else
               u1 = (-spi2_(PMom2,u1) + Quarks(1)%Mass*u1 )*PropFac2
            endif

            rIn = na+1
            rOut= NumGlu(1)+nc
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,(/nb+nc,nb,nc/))
            PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom
            PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(1)%Mass2)
            if ( nb .ge. 1 .or. nc .ge. 1) then
               if( abs(sc_(PMom3,PMom3) - Quarks(1)%Mass2).lt.PropCut ) cycle
               if( Quarks(1)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom3,ubar2) + Quarks(1)%Mass*ubar2 )*PropFac3
               else
                  ubar2 = (+spb2_(ubar2,PMom3) + Quarks(1)%Mass*ubar2 )*PropFac3
               endif
            endif

            if( Quarks(1)%PartType.lt.0 ) then
                Eps1 = -vbqq(Dv,u1,ubar2)       ! re-checked
            else
                Eps1 = +vbqq(Dv,ubar2,u1)       ! re-checked
            endif

            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2
         enddo
         enddo
         enddo
      endif! Bosonvertex


      if( BosonVertex.eq.3)  then
!        type (2)
         do na=0,NumGlu(1)
         do nc=0,NumGlu(4)
         do ne=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(4)-nc
            nf=NumGlu(5)-ne

            rIn = na+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nc
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(1:3),Quarks(4)%PartType,(/nb+NumGlu(2)+NumGlu(3)+nc,nb,NumGlu(2),NumGlu(3),nc/),0,0)
            PMom2 =  SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom + Quarks(2)%Mom + Quarks(3)%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
            if( Quarks(4)%PartType.lt.0 ) then
              u1 = (+spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
            else
              u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
            endif

            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            BosonVertex_mod = BosonVertex - 3
            ubar2 = cur_f_2fV(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,Boson,(/nd+ne,nd,ne/))
            PMom3  = SumMom(Gluons,rIn,rOut) + Quarks(4)%Mom + Boson%Mom
            PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(4)%Mass2)
            if( abs(sc_(PMom3,PMom3) - Quarks(4)%Mass2).lt.PropCut ) cycle
            if( Quarks(4)%PartType.lt.0 ) then
               ubar2 = (-spi2_(PMom3,ubar2) + Quarks(4)%Mass*ubar2 )*PropFac3
            else
               ubar2 = (+spb2_(ubar2,PMom3) + Quarks(4)%Mass*ubar2 )*PropFac3
            endif

            if( Quarks(4)%PartType.lt.0 ) then
              Eps1 = +vbqq(Dv,u1,ubar2)       ! re-checked
            else
              Eps1 = -vbqq(Dv,ubar2,u1)       ! re-checked
            endif

            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2
         enddo
         enddo
         enddo
      endif! Bosonvertex


      if (BosonVertex .ne. 1 .and. BosonVertex .ne. 2 .and. BosonVertex .ne. 3) then
         print *, 'cur_g_4fV not implemented for this flavor choice and BosonVertex=', BosonVertex
         stop
      endif

      endif! Flavor check

      if( Quarks(1)%PartType.eq.-Quarks(2)%PartType .and. Quarks(3)%PartType.eq.-Quarks(4)%PartType) then

         if (BosonVertex .eq. 1 .or. BosonVertex .eq. 3) then
             do na=0,NumGlu(1)
             do nc=0,NumGlu(2)
             do ne=0,NumGlu(5)

                nb=NumGlu(1)-na
                nd=NumGlu(2)-nc
                nf=NumGlu(5)-ne
                rIn = NumGlu(1)+nc+1
                rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
                u1 = cur_f_4fV(Gluons(rIn:rOut),Quarks(2:4),Quarks(1)%PartType,Boson, BosonVertex,(/nd+NumGlu(3)+NumGlu(4)+ne,nd,NumGlu(3),NumGlu(4),ne/),0,0)
                PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Boson%Mom
                PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(1)%Mass2)
                if( abs(sc_(PMom2,PMom2) - Quarks(1)%Mass2).lt.PropCut ) cycle
                if( Quarks(1)%PartType.lt.0 ) then
                   u1 = ( spb2_(u1,PMom2) + Quarks(1)%Mass*u1 )*PropFac2
                else
                   u1 = (-spi2_(PMom2,u1) + Quarks(1)%Mass*u1 )*PropFac2
                endif

                rIn = na+1
                rOut= NumGlu(1)+nc
                ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,(/nb+nc,nb,nc/))
                PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom
                if( nb.ge.1 .or. nc.ge.1 ) then
                   PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(1)%Mass2)
                   if( abs(sc_(PMom3,PMom3) - Quarks(1)%Mass2).lt.PropCut ) cycle
                   if( Quarks(1)%PartType.lt.0 ) then
                      ubar2 = (-spi2_(PMom3,ubar2) + Quarks(1)%Mass*ubar2 )*PropFac3
                   else
                      ubar2 = (+spb2_(ubar2,PMom3) + Quarks(1)%Mass*ubar2 )*PropFac3
                   endif
                endif

                if( Quarks(1)%PartType.lt.0 ) then
                   Eps1 = -vbqq(Dv,u1,ubar2)       ! re-checked
                else
                   Eps1 = +vbqq(Dv,ubar2,u1)       ! re-checked
                endif

                counter=1
                rIn =1
                rOut=na
                do i=rIn,rOut
                   call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                   counter=counter+1
                enddo
                TmpMom1(:) = PMom2(:)+PMom3(:)
                TmpExtRef1 = -1
                TmpGluons(counter)%Mom => TmpMom1(:)
                TmpGluons(counter)%Pol => Eps1(:)
                TmpGluons(counter)%ExtRef => TmpExtRef1
                counter=counter+1
                rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
                rOut=NumGlu(0)-1
                do i=rIn,rOut
                   call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                   counter=counter+1
                enddo
                Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

                if( na.ge.1 .or. nf.ge.1 ) then
                   PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
                   if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
                   Eps2 = Eps2*PropFac1
                endif

                Res = Res + Eps2
             enddo
             enddo
          enddo
       endif

       if (BosonVertex .eq. 1) then
          do na=0,NumGlu(1)
          do nc=0,NumGlu(2)
          do ne=0,NumGlu(5)
             nb=NumGlu(1)-na
             nd=NumGlu(2)-nc
             nf=NumGlu(5)-ne
             rIn = NumGlu(1)+nc+1
             rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
             u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(2:4),Quarks(1)%PartType,(/nd+NumGlu(3)+NumGlu(4)+ne,nd,NumGlu(3),NumGlu(4),ne/),0,0)
             PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom
             PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(1)%Mass2)
             if( abs(sc_(PMom2,PMom2) - Quarks(1)%Mass2).lt.PropCut ) cycle
             if( Quarks(1)%PartType.lt.0 ) then
                u1 = ( spb2_(u1,PMom2) + Quarks(1)%Mass*u1 )*PropFac2
             else
                u1 = (-spi2_(PMom2,u1) + Quarks(1)%Mass*u1 )*PropFac2
             endif

             rIn = na+1
             rOut= NumGlu(1)+nc
             ubar2 = cur_f_2fV(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,Boson,(/nb+nc,nb,nc/))
             PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom +  Boson%Mom
             PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(1)%Mass2)
             if( abs(sc_(PMom3,PMom3) - Quarks(1)%Mass2).lt.PropCut ) cycle
             if( Quarks(1)%PartType.lt.0 ) then
                ubar2 = (-spi2_(PMom3,ubar2) + Quarks(1)%Mass*ubar2 )*PropFac3
             else
                ubar2 = (+spb2_(ubar2,PMom3) + Quarks(1)%Mass*ubar2 )*PropFac3
             endif

             if( Quarks(1)%PartType.lt.0 ) then
                Eps1 = -vbqq(Dv,u1,ubar2)       ! re-checked
             else
                Eps1 = +vbqq(Dv,ubar2,u1)       ! re-checked
             endif

             counter=1
             rIn =1
             rOut=na
             do i=rIn,rOut
                call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                counter=counter+1
             enddo
             TmpMom1(:) = PMom2(:)+PMom3(:)
             TmpExtRef1 = -1
             TmpGluons(counter)%Mom => TmpMom1(:)
             TmpGluons(counter)%Pol => Eps1(:)
             TmpGluons(counter)%ExtRef => TmpExtRef1
             counter=counter+1
             rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
             rOut=NumGlu(0)-1
             do i=rIn,rOut
                call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                counter=counter+1
             enddo
             Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

             if( na.ge.1 .or. nf.ge.1 ) then
                PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
                if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
                Eps2 = Eps2*PropFac1
             endif

             Res = Res + Eps2
          enddo
       enddo
    enddo
 endif


 if (BosonVertex .eq. 1 .or. BosonVertex .eq. 3) then
!        type (3)
    do na=0,NumGlu(1)
    do nc=0,NumGlu(4)
    do ne=0,NumGlu(5)        ! can be replaced by above ne-loop
       nb=NumGlu(1)-na
       nd=NumGlu(4)-nc
       nf=NumGlu(5)-ne

       rIn = na+1
       rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nc
       BosonVertex_mod=BosonVertex+1
       u1 = cur_f_4fV(Gluons(rIn:rOut),Quarks(1:3),Quarks(4)%PartType,Boson,BosonVertex_mod,(/nb+NumGlu(2)+NumGlu(3)+nc,nb,NumGlu(2),NumGlu(3),nc/),0,0)
       PMom2 =  SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom + Quarks(2)%Mom + Quarks(3)%Mom + Boson%Mom
       PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
       if ( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) then
          cycle
       endif
       if( Quarks(4)%PartType.lt.0 ) then
          u1 = (+spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
       else
          u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
       endif
       rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nc+1
       rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
       ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/nd+ne,nd,ne/))
       PMom3  = SumMom(Gluons,rIn,rOut) + Quarks(4)%Mom
       if( nd.ge.1 .or. ne.ge.1 ) then
          PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(4)%Mass2)
          if( abs(sc_(PMom3,PMom3) - Quarks(4)%Mass2).lt.PropCut ) cycle
          if( Quarks(4)%PartType.lt.0 ) then
             ubar2 = (-spi2_(PMom3,ubar2) + Quarks(4)%Mass*ubar2 )*PropFac3
          else
             ubar2 = (+spb2_(ubar2,PMom3) + Quarks(4)%Mass*ubar2 )*PropFac3
          endif
       endif

       if( Quarks(4)%PartType.lt.0 ) then
          Eps1 = +vbqq(Dv,u1,ubar2)       ! re-checked
       else
          Eps1 = -vbqq(Dv,ubar2,u1)       ! re-checked
       endif

       counter=1
       rIn =1
       rOut=na
       do i=rIn,rOut
          call CopyParticlePtr(Gluons(i),TmpGluons(counter))
          counter=counter+1
       enddo
       TmpMom1(:) = PMom2(:)+PMom3(:)
       TmpExtRef1 = -1
       TmpGluons(counter)%Mom => TmpMom1(:)
       TmpGluons(counter)%Pol => Eps1(:)
       TmpGluons(counter)%ExtRef => TmpExtRef1
       counter=counter+1
       rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
       rOut=NumGlu(0)-1
       do i=rIn,rOut
          call CopyParticlePtr(Gluons(i),TmpGluons(counter))
          counter=counter+1
       enddo
       Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

       if( na.ge.1 .or. nf.ge.1 ) then
          PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
          if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
          Eps2 = Eps2*PropFac1
       endif

       Res = Res + Eps2
    enddo
    enddo
    enddo
   endif

   if (BosonVertex .eq. 1) then
! type(4)
      do na=0,NumGlu(1)
      do nc=0,NumGlu(2)
      do ne=0,NumGlu(3)
      do nf=0,NumGlu(3)-ne
      do nh=0,NumGlu(4)   ! this loop can be placed after Eps1 has been calculated
      do nj=0,NumGlu(5)
         nb=NumGlu(1)-na
         nd=NumGlu(2)-nc
         ng=NumGlu(3)-ne-nf
         ni=NumGlu(4)-nh
         nk=NumGlu(5)-nj

         rIn = na+1
         rOut= NumGlu(1)+NumGlu(2)+ne
         eps1 = cur_g_2fV(Gluons(rIn:rOut),Quarks(1:2),Boson,(/1+nb+nc+nd+ne,nc,nc+nd,ne/) )
         TmpMom1 =  SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom + Quarks(2)%Mom + Boson%Mom
         PropFac3 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
         if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
         Eps1 = Eps1*PropFac3

         rIn = NumGlu(1)+NumGlu(2)+ne+nf+1
         rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nh
         u1 = cur_f_2f(Gluons(rIn:rOut),Quarks(3:3),-Quarks(3)%PartType,(/ng+nh,ng,nh/))
         PMom3 = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom
         if( ng.ge.1 .or. nh.ge.1 ) then
            PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(3)%Mass2)
            if( abs(sc_(PMom3,PMom3) - Quarks(3)%Mass2).lt.PropCut ) cycle
            if( Quarks(3)%PartType.lt.0 ) then
               u1 = (-spi2_(PMom3,u1) + Quarks(3)%Mass*u1 )*PropFac3
            else
               u1 = (+spb2_(u1,PMom3) + Quarks(3)%Mass*u1 )*PropFac3
            endif
         endif
         rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nh+1
         rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj
         ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,(/ni+nj,ni,nj/))
         PMom4 = SumMom(Gluons,rIn,rOut) + Quarks(4)%Mom
         if( ni.ge.1 .or. nj.ge.1 ) then
            PropFac4 = (0d0,1d0)/(sc_(PMom4,PMom4) - Quarks(4)%Mass2)
            if( abs(sc_(PMom4,PMom4) - Quarks(4)%Mass2).lt.PropCut ) cycle
            if( Quarks(4)%PartType.lt.0 ) then
               ubar2 = (-spi2_(PMom4,ubar2) + Quarks(4)%Mass*ubar2 )*PropFac4
            else
               ubar2 = (+spb2_(ubar2,PMom4) + Quarks(4)%Mass*ubar2 )*PropFac4
            endif
         endif

         if( Quarks(4)%PartType.lt.0 ) then
            Eps2 = +vbqq(Dv,u1,ubar2)       ! re-checked
         else
            Eps2 = -vbqq(Dv,ubar2,u1)       ! re-checked
         endif
         TmpMom2 = PMom3 + PMom4
         PropFac1 = (0d0,-1d0)/sc_(TmpMom2,TmpMom2)
         if( abs(sc_(TmpMom2,TmpMom2)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1


         counter=1
         rIn =1
         rOut=na
         do i=rIn,rOut
            call CopyParticlePtr(Gluons(i),TmpGluons(counter))
            counter=counter+1
         enddo
         TmpExtRef1 = -1
         TmpGluons(counter)%Mom => TmpMom1(:)
         TmpGluons(counter)%Pol => Eps1(:)
         TmpGluons(counter)%ExtRef => TmpExtRef1
         counter=counter+1
         rIn =NumGlu(1)+NumGlu(2)+ne+1
         rOut=NumGlu(1)+NumGlu(2)+ne+nf
         do i=rIn,rOut
            call CopyParticlePtr(Gluons(i),TmpGluons(counter))
            counter=counter+1
         enddo
         TmpGluons(counter)%Mom => TmpMom2(:)
         TmpGluons(counter)%Pol => Eps2(:)
         TmpGluons(counter)%ExtRef => TmpExtRef1
         counter=counter+1
         rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj+1
         rOut=NumGlu(0)-1
         do i=rIn,rOut
            call CopyParticlePtr(Gluons(i),TmpGluons(counter))
            counter=counter+1
         enddo
         Eps3(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+nk+2)
         Res = Res + Eps3
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
      endif

      if (BosonVertex .eq. 3) then

         do na=0,NumGlu(1)
         do nc=0,NumGlu(4)
         do ne=0,NumGlu(5)        ! can be replaced by above ne-loop
            nb=NumGlu(1)-na
            nd=NumGlu(4)-nc
            nf=NumGlu(5)-ne

            rIn = na+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+nc
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(1:3),Quarks(4)%PartType,(/nb+NumGlu(2)+NumGlu(3)+nc,nb,NumGlu(2),NumGlu(3),nc/),0,0)
            PMom2 =  SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom + Quarks(2)%Mom + Quarks(3)%Mom
            PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
            if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
            if( Quarks(4)%PartType.lt.0 ) then
              u1 = (+spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
            else
              u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
            endif
            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne
            ubar2 = cur_f_2fV(Gluons(rIn:rOut),Quarks(4:4),-Quarks(4)%PartType,Boson,(/nd+ne,nd,ne/))
            PMom3  = SumMom(Gluons,rIn,rOut) + Quarks(4)%Mom + Boson%Mom
            PropFac3 = (0d0,1d0)/(sc_(PMom3,PMom3) - Quarks(4)%Mass2)
            if( abs(sc_(PMom3,PMom3) - Quarks(4)%Mass2).lt.PropCut ) cycle
            if( Quarks(4)%PartType.lt.0 ) then
               ubar2 = (-spi2_(PMom3,ubar2) + Quarks(4)%Mass*ubar2 )*PropFac3
            else
               ubar2 = (+spb2_(ubar2,PMom3) + Quarks(4)%Mass*ubar2 )*PropFac3
            endif

            if( Quarks(4)%PartType.lt.0 ) then
              Eps1 = +vbqq(Dv,u1,ubar2)       ! re-checked
            else
              Eps1 = -vbqq(Dv,ubar2,u1)       ! re-checked
            endif

            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpMom1(:) = PMom2(:)+PMom3(:)
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+ne+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps2(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+1)

            if( na.ge.1 .or. nf.ge.1 ) then
               PropFac1 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
               if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
               Eps2 = Eps2*PropFac1
            endif

            Res = Res + Eps2
         enddo
         enddo
         enddo
      endif

      if (BosonVertex .eq. 3) then

         do na=0,NumGlu(1)
         do nc=0,NumGlu(2)
         do ne=0,NumGlu(3)
         do nf=0,NumGlu(3)-ne
         do nh=0,NumGlu(4)   ! this loop can be placed after Eps1 has been calculated
         do nj=0,NumGlu(5)
            nb=NumGlu(1)-na
            nd=NumGlu(2)-nc
            ng=NumGlu(3)-ne-nf
            ni=NumGlu(4)-nh
            nk=NumGlu(5)-nj

            rIn = na+1
            rOut= NumGlu(1)+nc
            ubar2 = cur_f_2f(Gluons(rIn:rOut),Quarks(1:1),-Quarks(1)%PartType,(/nb+nc,nb,nc/))
            PMom1 = SumMom(Gluons,rIn,rOut) + Quarks(1)%Mom
            if( nb.ge.1 .or. nc.ge.1 ) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1) - Quarks(1)%Mass2)
               if( abs(sc_(PMom1,PMom1) - Quarks(1)%Mass2).lt.PropCut ) cycle
               if( Quarks(1)%PartType.lt.0 ) then
                  ubar2 = (-spi2_(PMom1,ubar2) + Quarks(1)%Mass*ubar2 )*PropFac1
               else
                  ubar2 = (+spb2_(ubar2,PMom1) + Quarks(1)%Mass*ubar2 )*PropFac1
               endif
            endif
            rIn = NumGlu(1)+nc+1
            rOut= NumGlu(1)+NumGlu(2)+ne
            u1 = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/nd+ne,nd,ne/))
            PMom2 = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom
            if( nd.ge.1 .or. ne.ge.1 ) then
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
               if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  u1 = (-spi2_(PMom2,u1) + Quarks(2)%Mass*u1 )*PropFac2
               else
                  u1 = (+spb2_(u1,PMom2) + Quarks(2)%Mass*u1 )*PropFac2
               endif
            endif

            if( Quarks(2)%PartType.lt.0 ) then
              Eps1 = +vbqq(Dv,ubar2,u1)       ! re-checked
            else
              Eps1 = -vbqq(Dv,u1,ubar2)       ! re-checked
            endif
            TmpMom1 = PMom1 + PMom2
            PropFac3 = (0d0,-1d0)/sc_(TmpMom1,TmpMom1)
            if( abs(sc_(TmpMom1,TmpMom1)).lt.PropCut ) cycle
            Eps1 = Eps1*PropFac3

            rIn = NumGlu(1)+NumGlu(2)+ne+nf+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj
            Eps2=cur_g_2fV(Gluons(rIn:rOut),Quarks(3:4),Boson,(/1+ng+nh+ni+nj,ng,nh,ni,nj /) )
            TmpMom2 = SumMom(Gluons,rIn,rOut) +Quarks(3)%Mom + Quarks(4)%Mom + Boson%Mom
            PropFac1 = (0d0,-1d0)/sc_(TmpMom2,TmpMom2)
            if( abs(sc_(TmpMom2,TmpMom2)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1

            counter=1
            rIn =1
            rOut=na
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpExtRef1 = -1
            TmpGluons(counter)%Mom => TmpMom1(:)
            TmpGluons(counter)%Pol => Eps1(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+ne+1
            rOut=NumGlu(1)+NumGlu(2)+ne+nf
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            TmpGluons(counter)%Mom => TmpMom2(:)
            TmpGluons(counter)%Pol => Eps2(:)
            TmpGluons(counter)%ExtRef => TmpExtRef1
            counter=counter+1
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+nj+1
            rOut=NumGlu(0)-1
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Eps3(:) = cur_g(TmpGluons(1:counter-1),1+na+nf+nk+2)
            Res = Res + Eps3

         enddo
      enddo
   enddo
enddo
enddo
enddo

endif   ! BosonVertex


if (BosonVertex .ne. 1 .and. BosonVertex .ne. 3) then
   print *, 'Not implemented for for cur_g_4fV for this flavor structure and BosonVertex=', BosonVertex
   stop
endif
endif

return
END FUNCTION











FUNCTION cur_f_6fV(Gluons,Quarks,Quark1PartType,Boson,BosonVertex,NumGlu,tag_f) result(res)           ! Quarks(:) does not include the OFF-shell quark
implicit none
integer :: NumGlu(0:6),Quark1PartType,tag_f,BosonVertex,BosonVertex_mod
type(PtrToParticle) :: Gluons(1:),Quarks(2:6),Boson
integer,target :: TmpPartType,TmpExtRef
complex(8) :: Res(1:Ds),tmp(1:Ds)
! complex(8) :: Res1(1:Ds),Res2(1:Ds),Res3(1:Ds),Res4(1:Ds)
complex(8) :: u1(1:Ds),ubar1(1:Ds)
complex(8),target :: ubar0(1:Ds)
complex(8) :: eps1(1:Dv)
complex(8) :: eps2(1:Dv)
type(PtrToParticle) :: TmpGluons(1:NumGlu(1)+NumGlu(6)),TmpQuark(1:1)
complex(8) :: PropFac1,PropFac2
complex(8),target :: PMom1(1:Dv)
complex(8),target :: PMom2(1:Dv)
integer :: n1a,n1b,n2a,n2b,n3a,n3b,n4a,n4b,n5a,n5b,n6a,n6b
integer :: rIn,rOut,i,counter


!DEC$ IF (_DebugCheckMyImpl1==1)
    if( NumGlu(0)-NumGlu(1)-NumGlu(2)-NumGlu(3)-NumGlu(4)-NumGlu(5)-NumGlu(6).ne.0 ) print *, "wrong number of gluons in cur_f_6fV"
!DEC$ ENDIF

    Res = (0d0,0d0)


! -- not all possible choices of flavors and BosonVertex values are used in ttb+Z calculation. This warns you if you try to use something new and potentially buggy.
    if (Quark1PartType.eq.-Quarks(2)%PartType) then
       if ( ((Quarks(3)%PartType.eq.-Quarks(4)%PartType) .AND. (BosonVertex .eq. 2 .OR.BosonVertex .eq. 4 .OR.BosonVertex .eq. 6)) .OR. &
         ((Quarks(3)%PartType.eq.-Quarks(6)%PartType) .AND. (BosonVertex .eq. 2 .OR.BosonVertex .eq. 6)) ) then
          print *, 'WARNING :  cur_f_6fV with this flavor structure and this choice of BosonVertex is implemented, but not checked!'
          print *, 'Quark flavors: ',Quark1PartType, Quarks(2)%PartType,Quarks(3)%PartType,Quarks(4)%PartType,Quarks(5)%PartType,Quarks(6)%PartType
          print *, 'BosonVertex =', BosonVertex
       endif
    elseif (Quark1PartType.eq.-Quarks(6)%PartType ) then
       if ( ((Quark1PartType.eq.-Quarks(2)%PartType) .AND. (BosonVertex .eq. 2 .or. BosonVertex .eq. 4 .or. BosonVertex .eq. 6)) .OR. &
            & ((Quark1PartType.eq.-Quarks(4)%PartType) .AND. (BosonVertex .eq. 4 .or. BosonVertex .eq. 6)) ) then
          print *, 'WARNING :  cur_f_6fV with this flavor structure and this choice of BosonVertex is implemented, but not checked!'
          print *, 'Quark flavors: ',Quark1PartType, Quarks(2)%PartType,Quarks(3)%PartType,Quarks(4)%PartType,Quarks(5)%PartType,Quarks(6)%PartType
          print *, 'BosonVertex =', BosonVertex
       endif
    elseif (Quarks(2)%PartType.eq.-Quarks(3)%PartType) then
       if ( ((Quark1PartType.eq.-Quarks(4)%PartType) .AND.  (BosonVertex .eq. 4 .or. BosonVertex .eq. 6) ) .OR. &
            ((Quark1PartType.eq.-Quarks(6)%PartType) .AND.  (BosonVertex .eq. 4 .or. BosonVertex .eq. 6) ) ) then
          print *, 'WARNING :  cur_f_6fV with this flavor structure and this choice of BosonVertex is implemented, but not checked!'
          print *, 'Quark flavors: ',Quark1PartType, Quarks(2)%PartType,Quarks(3)%PartType,Quarks(4)%PartType,Quarks(5)%PartType,Quarks(6)%PartType
          print *, 'BosonVertex =', BosonVertex
       endif
    endif


! probably some tag_Z checks are needed here: not yet implemented


!   (A)
    if( Quark1PartType.eq.-Quarks(2)%PartType .AND. (Quarks(3)%PartType.eq.-Quarks(4)%PartType .or. Quarks(3)%PartType.eq.-Quarks(6)%PartType) &
         .AND. (Quarks(2)%ExtRef.ne.-1 .or. tag_f.ne.1) &
         .AND. ( BosonVertex.eq.1 .OR. BosonVertex.eq.2) &
         ) then
!       print *, 'A-1'
      do n2a=0,NumGlu(2)
      do n6a=0,NumGlu(6)
         n2b = NumGlu(2)-n2a
         n6b = NumGlu(6)-n6a

         rIn =NumGlu(1)+n2a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
         Eps2 = cur_g_4f(Gluons(rIn:rOut),Quarks(3:6),(/1+n2b+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a,n2b,NumGlu(3),NumGlu(4),NumGlu(5),n6a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut) cycle
         Eps2 = Eps2*PropFac1
         do n1a=0,NumGlu(1)
            n1b = NumGlu(1)-n1a
            ! Fer2 couple V on the top
            rIn =n1a+1
            rOut=NumGlu(1)+n2a
            if (BosonVertex .eq. 1 .or. BosonVertex .eq. 2) then
               ubar1(:) = cur_f_2fV(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,Boson,(/n1b+n2a,n1b,n2a/) )
               PMom2(:) = Quarks(2)%Mom + SumMom(Gluons,rIn,rOut) + Boson%Mom
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) then
                  PropFac2=(0d0,0d0)
               endif

               if( Quarks(2)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(2)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(2)%Mass*ubar1(:))*PropFac2
               endif
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar0(:) = vbqg(ubar1,eps2)
               else
                  ubar0(:) = vqg(ubar1,eps2)
               endif

               rIn = n1a+1
               rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
               PMom1 = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom + Boson%Mom  ! can be simplified with PMom1(:)
               if(n1a.ge.1 .or. n6b.ge.1) then
                  PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
                  if( abs(sc_(PMom1,PMom1)-Quarks(2)%Mass2).lt.PropCut ) then
                     PropFac1=(0d0,0d0)
                  endif

                  if( Quarks(2)%PartType.lt.0 ) then
                     ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(2)%Mass*ubar0(:))*PropFac1
                  else
                     ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(2)%Mass*ubar0(:))*PropFac1
                  endif
               endif

               TmpQuark(1)%Mom  => PMom1(:)
               TmpQuark(1)%Pol  => ubar0(:)
               TmpQuark(1)%Mass => Quarks(2)%Mass
               TmpQuark(1)%Mass2=> Quarks(2)%Mass2
               TmpExtRef = -1
               TmpQuark(1)%ExtRef => TmpExtRef
               TmpQuark(1)%PartType => Quarks(2)%PartType
               counter=1
               rIn =1
               rOut=n1a
               do i=rIn,rOut
                  call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                  counter=counter+1
               enddo
               rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
               rOut=NumGlu(0)
               do i=rIn,rOut
                  call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                  counter=counter+1
               enddo
               tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n6b/) )
               Res(:) = Res(:) + tmp(:)
            endif
            if (BosonVertex .eq. 1 .or. BosonVertex .eq. 6) then
               ! Fer2  couple V on the bottom
               rIn =n1a+1
               rOut=NumGlu(1)+n2a
               ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n1b+n2a,n1b,n2a/) )
               if(n1b.ge.1 .or. n2a.ge.1) then
                  PMom2(:) = Quarks(2)%Mom + SumMom(Gluons,rIn,rOut)
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
                  if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) cycle
                  if( Quarks(2)%PartType.lt.0 ) then
                     ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(2)%Mass*ubar1(:))*PropFac2
                  else
                     ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(2)%Mass*ubar1(:))*PropFac2
                  endif
               endif
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar0(:) = vbqg(ubar1,eps2)       ! re-checked
               else
                  ubar0(:) = vqg(ubar1,eps2)       ! re-checked
               endif

               rIn = n1a+1
               rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
               PMom1 = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom

               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(2)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(2)%Mass*ubar0(:))*PropFac1
                  endif

                  TmpQuark(1)%Mom  => PMom1(:)
                  TmpQuark(1)%Pol  => ubar0(:)
                  TmpQuark(1)%Mass => Quarks(2)%Mass
                  TmpQuark(1)%Mass2=> Quarks(2)%Mass2
                  TmpExtRef = -1
                  TmpQuark(1)%ExtRef => TmpExtRef
                  TmpQuark(1)%PartType => Quarks(2)%PartType
                  counter=1
                  rIn =1
                  rOut=n1a
                  do i=rIn,rOut
                     call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                     counter=counter+1
                  enddo
                  rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
                  rOut=NumGlu(0)
                  do i=rIn,rOut
                     call CopyParticlePtr(Gluons(i),TmpGluons(counter))
                     counter=counter+1
                  enddo
                  tmp(:) = cur_f_2fV(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,Boson,(/counter-1,n1a,n6b/) )
                  Res(:) = Res(:) + tmp(:)
               endif
            enddo
         enddo
      enddo
   endif



!   (A)
    if( Quark1PartType.eq.-Quarks(2)%PartType .AND. (Quarks(3)%PartType.eq.-Quarks(4)%PartType .or. Quarks(3)%PartType.eq.-Quarks(6)%PartType) &
        .AND. (Quarks(2)%ExtRef.ne.-1 .or. tag_f.ne.1) &
        .AND. ( BosonVertex.eq.3 .OR. BosonVertex .eq. 4 .OR. BosonVertex.eq.5 ) &
      ) then
!       print *, 'A-2'
      do n2a=0,NumGlu(2)
      do n6a=0,NumGlu(6)
         n2b = NumGlu(2)-n2a
         n6b = NumGlu(6)-n6a

         rIn =NumGlu(1)+n2a+1
         rOut=NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
         BosonVertex_mod = BosonVertex - 2
         Eps2 = cur_g_4fV(Gluons(rIn:rOut),Quarks(3:6),Boson,BosonVertex_mod,(/1+n2b+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a,n2b,NumGlu(3),NumGlu(4),NumGlu(5),n6a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(3)%Mom + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom + Boson%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut) cycle
         Eps2 = Eps2*PropFac1
         do n1a=0,NumGlu(1)
            n1b = NumGlu(1)-n1a
            ! Fer2
            rIn =n1a+1
            rOut=NumGlu(1)+n2a
            ubar1(:) = cur_f_2f(Gluons(rIn:rOut),Quarks(2:2),-Quarks(2)%PartType,(/n1b+n2a,n1b,n2a/) )
            if(n1b.ge.1 .or. n2a.ge.1) then
               PMom2(:) = Quarks(2)%Mom + SumMom(Gluons,rIn,rOut)
               PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2)-Quarks(2)%Mass2)
               if( abs(sc_(PMom2,PMom2)-Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar1(:) = (-spi2_(PMom2,ubar1)+Quarks(2)%Mass*ubar1(:))*PropFac2
               else
                  ubar1(:) = (+spb2_(ubar1,PMom2)+Quarks(2)%Mass*ubar1(:))*PropFac2
               endif
            endif
            if( Quarks(2)%PartType.lt.0 ) then
               ubar0(:) = vbqg(ubar1,eps2)
            else
               ubar0(:) = vqg(ubar1,eps2)
            endif

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            PMom1 = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom + Boson%Mom
            if(n1a.ge.1 .or. n6b.ge.1) then
               PropFac1 = (0d0,1d0)/(sc_(PMom1,PMom1)-Quarks(2)%Mass2)
               if( abs(sc_(PMom1,PMom1)-Quarks(2)%Mass2).lt.PropCut ) cycle
               if( Quarks(2)%PartType.lt.0 ) then
                  ubar0(:) = (-spi2_(PMom1,ubar0)+Quarks(2)%Mass*ubar0(:))*PropFac1
               else
                  ubar0(:) = (+spb2_(ubar0,PMom1)+Quarks(2)%Mass*ubar0(:))*PropFac1
               endif
            endif

            TmpQuark(1)%Mom  => PMom1(:)
            TmpQuark(1)%Pol  => ubar0(:)
            TmpQuark(1)%Mass => Quarks(2)%Mass
            TmpQuark(1)%Mass2=> Quarks(2)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpQuark(1)%PartType => Quarks(2)%PartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),-TmpQuark(1)%PartType,(/counter-1,n1a,n6b/) )
            Res(:) = Res(:) + tmp(:)

         enddo
      enddo
      enddo
   endif



!   (B)
!! RR -- I think that only times this will be called in ttbZ will be color-zero in any case (tadpoles), so I'm going to just comment it out.
!! I imagine that if implemented correctly, then the it would return zero for ttbZ for kinematic reasons, but there doesn't seem to be any poin in coding it up for now...

!    if( Quark1PartType.eq.-Quarks(6)%PartType .AND. (Quarks(2)%PartType.eq.-Quarks(5)%PartType .or. Quarks(2)%PartType.eq.-Quarks(3)%PartType) &
!        .AND. (Quarks(6)%ExtRef.ne.-1 .or. tag_f.ne.1) &
!      ) then
!      call Error("This 6fV(B) current is not yet implemented")
!
!    endif

!   (C)
    if( Quarks(5)%PartType.eq.-Quarks(6)%PartType .AND. &
        ((Quark1PartType.eq.-Quarks(2)%PartType .and. (Quarks(2)%ExtRef.ne.-1.or.tag_f.ne.1) ) &
    .OR. (Quark1PartType.eq.-Quarks(4)%PartType .and. (Quarks(4)%ExtRef.ne.-1.or.tag_f.ne.1) ))&
    .AND. ( BosonVertex.eq.1 .OR. BosonVertex.eq.6 )  &
      ) then
!       print *, 'C-1'
      do n1a=0,NumGlu(1)
      do n4a=0,NumGlu(4)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n4b = NumGlu(4)-n4a
         n6b = NumGlu(6)-n6a

            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(5:6),(/1+n4b+NumGlu(5)+n6a,n4b,NumGlu(5),n6a/))
            PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
            PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
            if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(2:4),Quark1PartType,(/n1b+NumGlu(2)+NumGlu(3)+n4a,n1b,NumGlu(2),NumGlu(3),n4a/),0,0)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom

            if( Quark1PartType.eq.-Quarks(2)%PartType) then
                 if( Quarks(2)%PartType.lt.0 ) then
                   PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                   if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                   u1 = (-spi2_(PMom2,u1) + Quarks(2)%Mass*u1 )*PropFac2
                   ubar0 = vbqg(u1,eps2)       ! re-checked
                   rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                   rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                   PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom   ! this PMom2 will be re-used below

                   if( n1a.ge.1 .or. n6b.ge.1 .or. BosonVertex.eq.1 .or. BosonVertex.eq.6 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(2)%Mass*ubar0 )*PropFac2
                   endif
                elseif( Quarks(2)%PartType.gt.0 ) then
                   PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                   if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                   u1 = ( spb2_(u1,PMom2) + Quarks(2)%Mass*u1 )*PropFac2
                   ubar0 = vqg(u1,eps2)       ! re-checked
                   rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                   rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                   PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom

                   if( n1a.ge.1 .or. n6b.ge.1 .or. BosonVertex.eq.1 .or. BosonVertex.eq.6 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                       if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                       ubar0 = ( spb2_(ubar0,PMom2) + Quarks(2)%Mass*ubar0 )*PropFac2
                    endif
                 endif
                 TmpQuark(1)%Mom  => PMom2(:)
                 TmpQuark(1)%Pol  => ubar0(:)
                 TmpQuark(1)%Mass => Quarks(2)%Mass
                 TmpQuark(1)%Mass2=> Quarks(2)%Mass2
            elseif( Quark1PartType.eq.-Quarks(4)%PartType) then
                if( Quarks(4)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vbqg(u1,eps2)
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 .or. BosonVertex.eq.1 .or. BosonVertex.eq.6 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(4)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vqg(u1,eps2)
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 .or. BosonVertex.eq.1 .or. BosonVertex.eq.6 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(4)%Mass
                TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            endif
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpPartType = -Quark1PartType
            TmpQuark(1)%PartType => TmpPartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2fV(TmpGluons(1:counter-1),TmpQuark(1:1),Quark1PartType,Boson,(/counter-1,n1a,n6b/) )

            Res(:) = Res(:) + tmp(:)

      enddo
      enddo
      enddo
    endif



    if( Quarks(5)%PartType.eq.-Quarks(6)%PartType .AND. &
        ((Quark1PartType.eq.-Quarks(2)%PartType .and. (Quarks(2)%ExtRef.ne.-1.or.tag_f.ne.1) ) &
    .OR. (Quark1PartType.eq.-Quarks(4)%PartType .and. (Quarks(4)%ExtRef.ne.-1.or.tag_f.ne.1) ))&
    .AND. ( BosonVertex.eq.1 .OR. BosonVertex.eq.2 .OR. BosonVertex.eq.3 .OR. BosonVertex.eq.4  )  &
      ) then
!       print *, 'C-2'
      do n1a=0,NumGlu(1)
      do n4a=0,NumGlu(4)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n4b = NumGlu(4)-n4a
         n6b = NumGlu(6)-n6a

            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(5:6),(/1+n4b+NumGlu(5)+n6a,n4b,NumGlu(5),n6a/))
            PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
            PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
            if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            BosonVertex_mod = BosonVertex
            u1 = cur_f_4fV(Gluons(rIn:rOut),Quarks(2:4),Quark1PartType,Boson,BosonVertex_mod,(/n1b+NumGlu(2)+NumGlu(3)+n4a,n1b,NumGlu(2),NumGlu(3),n4a/),0,0)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom + Boson%Mom

            if( Quark1PartType.eq.-Quarks(2)%PartType) then
                 if( Quarks(2)%PartType.lt.0 ) then
                   PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                   if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                   u1 = (-spi2_(PMom2,u1) + Quarks(2)%Mass*u1 )*PropFac2
                   ubar0 = vbqg(u1,eps2)       ! re-checked
                   rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                   rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                   PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom   ! this PMom2 will be re-used below
                   if( n1a.ge.1 .or. n6b.ge.1 ) then
                       PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                       if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                       ubar0 = (-spi2_(PMom2,ubar0) + Quarks(2)%Mass*ubar0 )*PropFac2
                   endif
                 elseif( Quarks(2)%PartType.gt.0 ) then
                   PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                   if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                   u1 = ( spb2_(u1,PMom2) + Quarks(2)%Mass*u1 )*PropFac2
                   ubar0 = vqg(u1,eps2)       ! re-checked
                   rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                   rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                   PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
                   if( n1a.ge.1 .or. n6b.ge.1 ) then
                       PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                       if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                       ubar0 = ( spb2_(ubar0,PMom2) + Quarks(2)%Mass*ubar0 )*PropFac2
                   endif
                 endif
                 TmpQuark(1)%Mom  => PMom2(:)
                 TmpQuark(1)%Pol  => ubar0(:)
                 TmpQuark(1)%Mass => Quarks(2)%Mass
                 TmpQuark(1)%Mass2=> Quarks(2)%Mass2
            elseif( Quark1PartType.eq.-Quarks(4)%PartType) then
                if( Quarks(4)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vbqg(u1,eps2)
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(4)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vqg(u1,eps2)
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(4)%Mass
                TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            endif
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpPartType = -Quark1PartType
            TmpQuark(1)%PartType => TmpPartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),Quark1PartType,(/counter-1,n1a,n6b/) )

            Res(:) = Res(:) + tmp(:)

      enddo
      enddo
      enddo
    endif


!   (C)
    if( Quarks(5)%PartType.eq.-Quarks(6)%PartType .AND. &
        ((Quark1PartType.eq.-Quarks(2)%PartType .and. (Quarks(2)%ExtRef.ne.-1.or.tag_f.ne.1) ) &
    .OR. (Quark1PartType.eq.-Quarks(4)%PartType .and. (Quarks(4)%ExtRef.ne.-1.or.tag_f.ne.1) ))&
    .AND. ( BosonVertex.eq.5 )  &
      ) then
!       print *, 'C-3'
       do n1a=0,NumGlu(1)
       do n4a=0,NumGlu(4)
       do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n4b = NumGlu(4)-n4a
         n6b = NumGlu(6)-n6a

            rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            Eps2 = cur_g_2fV(Gluons(rIn:rOut),Quarks(5:6),Boson,(/1+n4b+NumGlu(5)+n6a,n4b,NumGlu(5),n6a/))
            PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom + Boson%Mom
            PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
            if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1

            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(2:4),Quark1PartType,(/n1b+NumGlu(2)+NumGlu(3)+n4a,n1b,NumGlu(2),NumGlu(3),n4a/),0,0)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(2)%Mom + Quarks(3)%Mom + Quarks(4)%Mom

            if( Quark1PartType.eq.-Quarks(2)%PartType) then

                 if( Quarks(2)%PartType.lt.0 ) then
                   PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                   if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                   u1 = (-spi2_(PMom2,u1) + Quarks(2)%Mass*u1 )*PropFac2
                   ubar0 = vbqg(u1,eps2)       ! re-checked
                   rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                   rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                   PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom   + Boson%Mom! this PMom2 will be re-used below
                   if( n1a.ge.1 .or. n6b.ge.1 ) then
                       PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                       if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                       ubar0 = (-spi2_(PMom2,ubar0) + Quarks(2)%Mass*ubar0 )*PropFac2
                   endif
                 elseif( Quarks(2)%PartType.gt.0 ) then
                    PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                    if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                    u1 = ( spb2_(u1,PMom2) + Quarks(2)%Mass*u1 )*PropFac2
                    ubar0 = vqg(u1,eps2)       ! re-checked
                   rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                   rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                   PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom   + Boson%Mom! this PMom2 will be re-used below
                   if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(2)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(2)%Mass2).lt.PropCut ) cycle
                      u1 = ( spb2_(u1,PMom2) + Quarks(2)%Mass*u1 )*PropFac2
                   endif
                endif
                 TmpQuark(1)%Mom  => PMom2(:)
                 TmpQuark(1)%Pol  => ubar0(:)
                 TmpQuark(1)%Mass => Quarks(2)%Mass
                 TmpQuark(1)%Mass2=> Quarks(2)%Mass2
            elseif( Quark1PartType.eq.-Quarks(4)%PartType) then
                if( Quarks(4)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vbqg(u1,eps2)
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(4)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vqg(u1,eps2)
                  rIn = NumGlu(1)+NumGlu(2)+NumGlu(3)+n4a+1
                  rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(5)%Mom + Quarks(6)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(4)%Mass
                TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            endif
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpPartType = -Quark1PartType
            TmpQuark(1)%PartType => TmpPartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),Quark1PartType,(/counter-1,n1a,n6b/) )

            Res(:) = Res(:) + tmp(:)
         enddo
      enddo
   enddo


endif



!   (D)
    if( Quarks(2)%PartType.eq.-Quarks(3)%PartType .AND. ( &
        (Quark1PartType.eq.-Quarks(4)%PartType .and. (Quarks(4)%ExtRef.ne.-1.or.tag_f.ne.1)) &
   .OR. (Quark1PartType.eq.-Quarks(6)%PartType .and. (Quarks(6)%ExtRef.ne.-1.or.tag_f.ne.1))  )  &
   .AND. ( BosonVertex.eq.1  )  &
     ) then
!       print *, 'D-1'
      do n1a=0,NumGlu(1)
      do n3a=0,NumGlu(3)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n3b = NumGlu(3)-n3a
         n6b = NumGlu(6)-n6a

            counter=1
            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+n3a
            Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(2:3),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
            PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
            PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
            if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1

            rIn = NumGlu(1)+NumGlu(2)+n3a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(4:6),Quark1PartType,(/n3b+NumGlu(4)+NumGlu(5)+n6a,n3b,NumGlu(4),NumGlu(5),n6a/),tag_f,0)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom

            if( Quark1PartType.eq.-Quarks(4)%PartType ) then
                if( Quarks(4)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgbq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom   !this PMom2 will be re-used below  ! CAN be written as PMom2=PMom2+PMom1
                  if( n1a.ge.1 .or. n6b.ge.1 .or. BosonVertex.eq.1 .or. BosonVertex.eq.6 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(4)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1  .or. BosonVertex.eq.1 .or. BosonVertex.eq.6 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(4)%Mass
                TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            elseif(Quark1PartType.eq.-Quarks(6)%PartType) then
               if( Quarks(6)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(6)%Mass*u1 )*PropFac2
                  ubar0 = vgbq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 .or. BosonVertex.eq.1 .or. BosonVertex.eq.6 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(6)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(6)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(6)%Mass*u1 )*PropFac2
                  ubar0 = vgq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1  .or. BosonVertex.eq.1 .or. BosonVertex.eq.6 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(6)%Mass*ubar0 )*PropFac2
                  endif
                endif
            endif
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpPartType = -Quark1PartType
            TmpQuark(1)%PartType => TmpPartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Tmp(:) = cur_f_2fV(TmpGluons(1:counter-1),TmpQuark(1:1),Quark1PartType,Boson,(/counter-1,n1a,n6b/))

            Res(:)  = Res(:) + Tmp(:)

      enddo
      enddo
      enddo
   endif



!   (D)
    if( Quarks(2)%PartType.eq.-Quarks(3)%PartType .AND. ( &
        (Quark1PartType.eq.-Quarks(4)%PartType .and. (Quarks(4)%ExtRef.ne.-1.or.tag_f.ne.1)) &
   .OR. (Quark1PartType.eq.-Quarks(6)%PartType .and. (Quarks(6)%ExtRef.ne.-1.or.tag_f.ne.1))  ) &
   .AND. ( BosonVertex.eq.2)  &
     ) then
!       print *, 'D-2'
      do n1a=0,NumGlu(1)
      do n3a=0,NumGlu(3)
      do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n3b = NumGlu(3)-n3a
         n6b = NumGlu(6)-n6a

            counter=1
            rIn = n1a+1
            rOut= NumGlu(1)+NumGlu(2)+n3a
            Eps2 = cur_g_2fV(Gluons(rIn:rOut),Quarks(2:3),Boson,(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
            PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom +  Boson%Mom
            PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
            if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
            Eps2 = Eps2*PropFac1

            rIn = NumGlu(1)+NumGlu(2)+n3a+1
            rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
            BosonVertex_mod = BosonVertex
            u1 = cur_f_4f(Gluons(rIn:rOut),Quarks(4:6),Quark1PartType,(/n3b+NumGlu(4)+NumGlu(5)+n6a,n3b,NumGlu(4),NumGlu(5),n6a/),tag_f,0)
            PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom

            if( Quark1PartType.eq.-Quarks(4)%PartType ) then
                if( Quarks(4)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgbq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom +  Boson%Mom   !this PMom2 will be re-used below  ! CAN be written as PMom2=PMom2+PMom1
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(4)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Boson%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                endif
            elseif(Quark1PartType.eq.-Quarks(6)%PartType) then
               if( Quarks(6)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(6)%Mass*u1 )*PropFac2
                  ubar0 = vgbq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom+  Boson%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(6)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(6)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(6)%Mass*u1 )*PropFac2
                  ubar0 = vgq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom + Boson%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(6)%Mass*ubar0 )*PropFac2
                  endif
                endif
            endif
            TmpQuark(1)%Mom  => PMom2(:)
            TmpQuark(1)%Pol  => ubar0(:)
            TmpQuark(1)%Mass => Quarks(4)%Mass
            TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpPartType = -Quark1PartType
            TmpQuark(1)%PartType => TmpPartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),Quark1PartType,(/counter-1,n1a,n6b/))

            Res(:)  = Res(:) + Tmp(:)
      enddo
      enddo
      enddo
    endif


!   (D)
    if( Quarks(2)%PartType.eq.-Quarks(3)%PartType .AND. ( &
        (Quark1PartType.eq.-Quarks(4)%PartType .and. (Quarks(4)%ExtRef.ne.-1.or.tag_f.ne.1)) &
   .OR. (Quark1PartType.eq.-Quarks(6)%PartType .and. (Quarks(6)%ExtRef.ne.-1.or.tag_f.ne.1))  ) &
   .AND. ( BosonVertex .ge. 3 .and. BosonVertex .le. 6  )  ) then
!       print *, 'D-3'
       do n1a=0,NumGlu(1)
       do n3a=0,NumGlu(3)
       do n6a=0,NumGlu(6)
         n1b = NumGlu(1)-n1a
         n3b = NumGlu(3)-n3a
         n6b = NumGlu(6)-n6a

         counter=1
         rIn = n1a+1
         rOut= NumGlu(1)+NumGlu(2)+n3a
         Eps2 = cur_g_2f(Gluons(rIn:rOut),Quarks(2:3),(/1+n1b+NumGlu(2)+n3a,n1b,NumGlu(2),n3a/))
         PMom1(:) = SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
         PropFac1 = (0d0,-1d0)/sc_(PMom1,PMom1)
         if( abs(sc_(PMom1,PMom1)).lt.PropCut ) cycle
         Eps2 = Eps2*PropFac1

         rIn = NumGlu(1)+NumGlu(2)+n3a+1
         rOut= NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a
         BosonVertex_mod=BosonVertex-2
         u1 = cur_f_4fV(Gluons(rIn:rOut),Quarks(4:6),Quark1PartType,Boson,BosonVertex_mod,(/n3b+NumGlu(4)+NumGlu(5)+n6a,n3b,NumGlu(4),NumGlu(5),n6a/),tag_f,0)
         PMom2  = SumMom(Gluons,rIn,rOut)  + Quarks(4)%Mom + Quarks(5)%Mom + Quarks(6)%Mom+Boson%Mom

            if( Quark1PartType.eq.-Quarks(4)%PartType ) then
                if( Quarks(4)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgbq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom   !this PMom2 will be re-used below  ! CAN be written as PMom2=PMom2+PMom1
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(4)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(4)%Mass*u1 )*PropFac2
                  ubar0 = vgq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(4)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(4)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(4)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(4)%Mass
                TmpQuark(1)%Mass2=> Quarks(4)%Mass2
            elseif(Quark1PartType.eq.-Quarks(6)%PartType) then
                if( Quarks(6)%PartType.lt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                  u1 = (-spi2_(PMom2,u1) + Quarks(6)%Mass*u1 )*PropFac2
                  ubar0 = vgbq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                      ubar0 = (-spi2_(PMom2,ubar0) + Quarks(6)%Mass*ubar0 )*PropFac2
                  endif
                elseif( Quarks(6)%PartType.gt.0 ) then
                  PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                  if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                  u1 = ( spb2_(u1,PMom2) + Quarks(6)%Mass*u1 )*PropFac2
                  ubar0 = vgq(eps2,u1)       ! re-checked
                  rIn = n1a+1
                  rOut= NumGlu(1)+NumGlu(2)+n3a
                  PMom2 = PMom2 + SumMom(Gluons,rIn,rOut) + Quarks(2)%Mom + Quarks(3)%Mom
                  if( n1a.ge.1 .or. n6b.ge.1 ) then
                      PropFac2 = (0d0,1d0)/(sc_(PMom2,PMom2) - Quarks(6)%Mass2)
                      if( abs(sc_(PMom2,PMom2) - Quarks(6)%Mass2).lt.PropCut ) cycle
                      ubar0 = ( spb2_(ubar0,PMom2) + Quarks(6)%Mass*ubar0 )*PropFac2
                  endif
                endif
                TmpQuark(1)%Mom  => PMom2(:)
                TmpQuark(1)%Pol  => ubar0(:)
                TmpQuark(1)%Mass => Quarks(6)%Mass
                TmpQuark(1)%Mass2=> Quarks(6)%Mass2
            endif
            TmpExtRef = -1
            TmpQuark(1)%ExtRef => TmpExtRef
            TmpPartType = -Quark1PartType
            TmpQuark(1)%PartType => TmpPartType
            counter=1
            rIn =1
            rOut=n1a
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            rIn =NumGlu(1)+NumGlu(2)+NumGlu(3)+NumGlu(4)+NumGlu(5)+n6a+1
            rOut=NumGlu(0)
            do i=rIn,rOut
              call CopyParticlePtr(Gluons(i),TmpGluons(counter))
              counter=counter+1
            enddo
            Tmp(:) = cur_f_2f(TmpGluons(1:counter-1),TmpQuark(1:1),Quark1PartType,(/counter-1,n1a,n6b/))

            Res(:)  = Res(:) + Tmp(:)
      enddo
      enddo
      enddo
    endif


return
END FUNCTION






!---------------------------------------


FUNCTION SumMom(Particles,i1,i2)
implicit none
complex(8) :: SumMom(1:Dv)
type(PtrToParticle) :: Particles(:)
integer :: i1,i2,j

   SumMom(1:Dv)= (0d0,0d0)
   if (i2.ge.i1) then
      do j=i1,i2
         SumMom(1:Dv) = SumMom(1:Dv) + Particles(j)%Mom(1:Dv)
      enddo
   endif
END FUNCTION



FUNCTION FourVecDot(p1,p2)
implicit none
complex(8), intent(in) :: p1(1:Dv),p2(1:Dv)
complex(8)  :: FourVecDot
integer :: mu

   FourVecDot = p1(1)*p2(1)
!DEC$ UNROLL
   do mu=2,Dv
      FourVecDot = FourVecDot - p1(mu)*p2(mu)
   enddo
return
END FUNCTION FourVecDot

FUNCTION eval_TripVert(k1,k2,v1,v2)
implicit none
complex(8) :: eval_TripVert(1:Dv)
complex(8) :: k1(1:Dv),k2(1:Dv),v1(1:Dv),v2(1:Dv)
complex(8), parameter :: IOverSqrt2=(0d0,1d0)/dsqrt(2d0)

   eval_TripVert(1:Dv) = IOverSqrt2 * ( (k1(1:Dv)-k2(1:Dv))*(v1.Ndot.v2)  &
                           - 2d0*v1(1:Dv)*(k1.Ndot.v2)  &
                           + 2d0*v2(1:Dv)*(k2.Ndot.v1) )
return
END FUNCTION eval_TripVert

FUNCTION eval_QuadVert(k1,k2,k3)
implicit none
complex(8) :: eval_QuadVert(1:Dv)
complex(8) :: k1(1:Dv),k2(1:Dv),k3(1:Dv)
complex(8), parameter :: I=(0d0,1d0)

   eval_QuadVert(1:Dv) = I * (-k1(1:Dv)*(k2.Ndot.k3)*0.5d0  &
                           + k2(1:Dv)*(k1.Ndot.k3)  &
                           - k3(1:Dv)*(k1.Ndot.k2)*0.5d0 )
return
END FUNCTION eval_QuadVert




!DEC$ ATTRIBUTES INLINE :: linear_map
FUNCTION linear_map(i1,i2,Ngluons)
implicit none
integer :: linear_map,i1,i2,Ngluons

   linear_map = i2+Ngluons*(i2-i1)-((i2-i1)*(i2-i1+1))/2
return
END FUNCTION linear_map






SUBROUTINE CopyParticlePtr(InPointer,OutPointer)
implicit none
type(PtrToParticle), intent(in) :: InPointer
type(PtrToParticle), intent(out):: OutPointer

   OutPointer%PartType => InPointer%PartType
   OutPointer%ExtRef => InPointer%ExtRef
   OutPointer%Mass => InPointer%Mass
   OutPointer%Mass2 => InPointer%Mass2
   OutPointer%Helicity => InPointer%Helicity
   OutPointer%Mom => InPointer%Mom
   OutPointer%Pol => InPointer%Pol
return
END SUBROUTINE







!---------------------------------------




!-----------modified procedure
!---------- Dv is the dimensionality of the vector space
!---------- Ds is the dimensionality of the spinorial representation
         subroutine spi2(Dv,Ds,v,sp,f)
         implicit none
         integer i,i1,i2,i3,imax,Dv,Ds
         double complex sp(Ds),v(Dv),f(Ds)
         double complex x0(4,4),xx(4,4),xy(4,4)
         double complex xz(4,4),x5(4,4)
         double complex y1,y2,y3,y4,bp,bm,cp,cm
         double complex test(Ds)


         imax = Ds/4

           do i=1,imax
           i1= 1+4*(i-1)
           i2=i1+3

           y1=sp(i1)
           y2=sp(i1+1)
           y3=sp(i1+2)
           y4=sp(i1+3)

           x0(1,i)=y1
           x0(2,i)=y2
           x0(3,i)=-y3
           x0(4,i)=-y4


           xx(1,i) = y4
           xx(2,i) = y3
           xx(3,i) = -y2
           xx(4,i) = -y1


           xy(1,i)=dcmplx(0d0,-1d0)*y4
           xy(2,i)=dcmplx(0d0,1d0)*y3
           xy(3,i)=dcmplx(0d0,1d0)*y2
           xy(4,i)=dcmplx(0d0,-1d0)*y1

           xz(1,i)=y3
           xz(2,i)=-y4
           xz(3,i)=-y1
           xz(4,i)=y2

           x5(1,i)=y3
           x5(2,i)=y4
           x5(3,i)=y1
           x5(4,i)=y2

           enddo

           if(Dv.eq.4) then

           do i=1,4

           f(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)
           enddo

           endif


          if (Dv.eq.6) then
           bp = (v(5)+dcmplx(0d0,1d0)*v(6))
           bm=(v(5)-dcmplx(0d0,1d0)*v(6))


           do i=1,4

           f(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)-bp*x5(i,2)

            i1=i+4

            f(i1)=v(1)*x0(i,2)-v(2)*xx(i,2)-v(3)*xy(i,2)-v(4)*xz(i,2)+bm*x5(i,1)

            enddo

          endif

          if (Dv.eq.8) then

           bp = (v(5)+dcmplx(0d0,1d0)*v(6))
           bm=(v(5)-dcmplx(0d0,1d0)*v(6))
           cp=(v(7)+dcmplx(0d0,1d0)*v(8))
           cm=(v(7)-dcmplx(0d0,1d0)*v(8))



           do i=1,4

           f(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)-bp*x5(i,2)+ cp*x5(i,3)

            i1=i+4

            f(i1)=v(1)*x0(i,2)-v(2)*xx(i,2)-v(3)*xy(i,2)-v(4)*xz(i,2)+bm*x5(i,1)-cp*x5(i,4)

             i2=i1+4

            f(i2)=v(1)*x0(i,3)-v(2)*xx(i,3)-v(3)*xy(i,3)-v(4)*xz(i,3)-bp*x5(i,4)-cm*x5(i,1)

             i3=i2+4

             f(i3)=v(1)*x0(i,4)-v(2)*xx(i,4)-v(3)*xy(i,4)-v(4)*xz(i,4)+bm*x5(i,3)+cm*x5(i,2)


            enddo

            endif

! if(Ds.eq.16) then
!   test(:) = VSpiL(v,sp)
!   print *, ""
!   print *, "spi2",f(1:Ds)
!   print *, "test",test(1:Ds)
!   print *, "diff",test(1:Ds)-f(1:Ds)
!   if( any( abs(test(1:Ds)-f(1:Ds) ).gt.1d-10  ) ) pause
! endif

           return
           end subroutine



         function spi2_(v,sp)
         implicit none
         double complex, intent(in) :: sp(:),v(:)
         double complex :: spi2_(size(sp)) ,tmp(size(sp))
         integer :: Dv,Ds

          Ds = size(sp)
          if (Ds == 4) Dv = 4
          if (Ds == 8) Dv = 6
          if (Ds == 16) Dv = 8
          call spi2(Dv,Ds,v,sp,spi2_)

!  tmp = VSpiL(v,sp)
!  print *, "diff1",tmp-spi2_
! pause

        return
        end function



      function SpiVL(sp,v)   ! SpiVL=sp.(v_mu*gamma^mu) =spb2(4,4,...)
      implicit none
      double complex :: sp(:),v(:)
      double complex :: SpiVL(size(sp))

            SpiVL(1) = sp(1)*v(1) + sp(4)*(v(2) + (0d0,1d0)*v(3)) + sp(3)*v(4)
            SpiVL(2) = sp(2)*v(1) + sp(3)*(v(2) - (0d0,1d0)*v(3)) - sp(4)*v(4)
            SpiVL(3) =-sp(3)*v(1) - sp(2)*(v(2) + (0d0,1d0)*v(3)) - sp(1)*v(4)
            SpiVL(4) =-sp(4)*v(1) - sp(1)*(v(2) - (0d0,1d0)*v(3)) + sp(2)*v(4)

      return
      end function



subroutine spb2(Dv,Ds,sp,v,f)
         implicit none
         integer i,i1,i2,i3,Dv,Ds,imax
         double complex sp(Ds),v(Dv),f(Ds)
         double complex x0(4,4),xx(4,4),xy(4,4)
         double complex xz(4,4),x5(4,4)
         double complex y1,y2,y3,y4,bp,bm,cp,cm
         double complex test(Ds)

           imax = Ds/4

           do i=1,imax
           i1= 1+4*(i-1)
           i2=i1+3

           y1=sp(i1)
           y2=sp(i1+1)
           y3=sp(i1+2)
           y4=sp(i1+3)

           x0(1,i)=y1
           x0(2,i)=y2
           x0(3,i)=-y3
           x0(4,i)=-y4

           xx(1,i) = -y4
           xx(2,i) = -y3
           xx(3,i) = y2
           xx(4,i) = y1

           xy(1,i)=dcmplx(0d0,-1d0)*y4
           xy(2,i)=dcmplx(0d0,1d0)*y3
           xy(3,i)=dcmplx(0d0,1d0)*y2
           xy(4,i)=dcmplx(0d0,-1d0)*y1

           xz(1,i)=-y3
           xz(2,i)=y4
           xz(3,i)=y1
           xz(4,i)=-y2

           x5(1,i)=y3
           x5(2,i)=y4
           x5(3,i)=y1
           x5(4,i)=y2

           enddo

           if (Dv.eq.4) then

           do i=1,4

           f(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)

           enddo

           endif

           if (Dv.eq.6) then
           bp = (v(5)+dcmplx(0d0,1d0)*v(6))
           bm=(v(5)-dcmplx(0d0,1d0)*v(6))

           do i=1,4

           f(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)+bm*x5(i,2)

            i1 = i+4

            f(i1)= v(1)*x0(i,2)-v(2)*xx(i,2)-v(3)*xy(i,2)-v(4)*xz(i,2)-bp*x5(i,1)


            enddo

           endif

           if (Dv.eq.8) then
           bp=(v(5)+dcmplx(0d0,1d0)*v(6))
           bm=(v(5)-dcmplx(0d0,1d0)*v(6))
           cp=(v(7)+dcmplx(0d0,1d0)*v(8))
           cm=(v(7)-dcmplx(0d0,1d0)*v(8))

           do i=1,4

           f(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)+bm*x5(i,2)-cm*x5(i,3)

            i1 = i+4

            f(i1)= v(1)*x0(i,2)-v(2)*xx(i,2)-v(3)*xy(i,2)-v(4)*xz(i,2)-bp*x5(i,1)+cm*x5(i,4)

             i2 = i1+4

             f(i2)=v(1)*x0(i,3)-v(2)*xx(i,3)-v(3)*xy(i,3)-v(4)*xz(i,3)+bm*x5(i,4)+cp*x5(i,1)

              i3=i2+4

              f(i3)=v(1)*x0(i,4)-v(2)*xx(i,4)-v(3)*xy(i,4)-v(4)*xz(i,4)-bp*x5(i,3)-cp*x5(i,2)

              enddo

              endif



! if(Ds.eq.16) then
!   test(:) = SpiVL(sp,v)
!   print *, ""
!   print *, "spb2",f(1:Ds)
!   print *, "test",test(1:Ds)
!   print *, "diff",test(1:Ds)-f(1:Ds)
!   if( any( abs(test(1:Ds)-f(1:Ds) ).gt.1d-10  ) ) pause
! endif




               return
end subroutine


        function spb2_(sp,v)
        implicit none
        double complex, intent(in) :: sp(:),v(:)
        double complex :: spb2_(size(sp)) ,tmp(size(sp))
        integer :: Dv,Ds

          Ds = size(sp)
          if (Ds == 4) Dv = 4
          if (Ds == 8) Dv = 6
          if (Ds == 16) Dv = 8
          call spb2(Dv,Ds,sp,v,spb2_)

!  tmp = SpiVL(sp(1:Ds),v(1:Dv))
!  print *, "diff2",tmp-spb2_
!  pause
        return
        end function




          function psp1_(sp1,sp2) result(res)
          implicit none
          complex(8), intent(in) :: sp1(:)
          complex(8), intent(in) :: sp2(:)
          complex(8) :: res

            res = sum(sp1(1:)*sp2(1:))

           end function




! RR -- new for D-dim chirality
             recursive function Chir(sign,sp) result(res)
               implicit none
               logical :: sign
               double complex :: sp(:)
               double complex :: res(size(sp))
               integer        :: D

               D = size(sp)
               if ( D .eq. 4) then
                  if(sign) then !omega_+
                     res(1) = 0.5d0*(sp(1)+sp(3))
                     res(2) = 0.5d0*(sp(2)+sp(4))
                     res(3) = res(1)
                     res(4) = res(2)
                  else !omega_-
                     res(1) = 0.5d0*(sp(1)-sp(3))
                     res(2) = 0.5d0*(sp(2)-sp(4))
                     res(3) =-res(1)
                     res(4) =-res(2)
                  endif
               else
                  res(1:D/2)     =  Chir(sign,sp(1:D/2))
                  res((D/2+1):D) =  Chir(sign,sp( (D/2+1):D ))
               endif

             end function Chir




             recursive function iChir(sign,sp) result(res)
! RR -- this function is needed for the electric and dipole moment like couplings,
!!      which have a 1 +\- i*gamma5 (see 0811.3842)
               implicit none
               logical :: sign
               double complex :: sp(:)
               double complex :: res(size(sp))
               double complex :: ci
               integer        :: D

               ci=(0d0,1d0)
               D = size(sp)
               if ( D .eq. 4) then
                  if(sign) then !omega_+
                     res(1) = 0.5d0*(sp(1)+ci*sp(3))
                     res(2) = 0.5d0*(sp(2)+ci*sp(4))
                     res(3) = 0.5d0*(ci*sp(1)+sp(3))
                     res(4) = 0.5d0*(ci*sp(2)+sp(4))
                  else !omega_-
                     res(1) = 0.5d0*(sp(1)-ci*sp(3))
                     res(2) = 0.5d0*(sp(2)-ci*sp(4))
                     res(3) = 0.5d0*(-ci*sp(1)+sp(3))
                     res(4) = 0.5d0*(-ci*sp(2)+sp(4))
                  endif
               else
                  res(1:D/2)     =  iChir(sign,sp(1:D/2))
                  res((D/2+1):D) =  iChir(sign,sp( (D/2+1):D ))
               endif

             end function iChir




      function vbqV(sp,e1,coupl_left,coupl_right)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8), intent(in) :: coupl_left,coupl_right
      complex(8) :: vbqV(size(sp))

!            vbqV = -(0d0,1d0)*( coupl_left*Chir(.false.,spb2_(sp,e1)) + coupl_right*Chir(.true.,spb2_(sp,e1)) )
           vbqV = -(0d0,1d0)*( coupl_left*Chir(.false.,sp) + coupl_right*Chir(.true.,sp) )

      return
      end function



      function vVq(e1,sp,coupl_left,coupl_right)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8), intent(in) :: coupl_left,coupl_right
      complex(8) :: vVq(size(sp))

!            vVq = -(0d0,1d0)*( coupl_left*Chir(.true., spi2_(e1,sp)) + coupl_right*Chir(.false., spi2_(e1,sp)) )
           vVq = -(0d0,1d0)*( coupl_left*Chir(.false.,sp) + coupl_right*Chir(.true.,sp) )


      return
      end function




          function sc_(p1,p2)
          implicit none
          complex(8) :: p1(:),p2(:)
          complex(8) :: sc_
          integer :: sizemin

              sizemin=min(size(p1),size(p2))
              call rsc_(sizemin,p1,p2,sc_)

          return
          end function



         subroutine rsc_(n,x,y,r)
         implicit none
         integer i,n
         complex(8) x(*),y(*)
         complex(8) r

            r = x(1)*y(1)
            do i=2, n
              r = r - x(i)*y(i)
            enddo

         return
         end subroutine






      function vggg(e1,k1,e2,k2)
      implicit none
      complex(8), intent(in) :: e1(:), e2(:)
      complex(8), intent(in) :: k1(:), k2(:)
      complex(8)             :: vggg(size(e1))
      complex(8):: sk1e2,se1e2,sk2e1,xx
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0


          sk1e2=sc_(k1,e2)
          sk2e1=sc_(k2,e1)
          se1e2=sc_(e1,e2)
          xx=(0.0d0,1.0d0)*sqrt2
          vggg = xx*(-sk1e2*e1+sk2e1*e2+se1e2/2d0*(k1-k2))

       end function vggg


       function  vgggg(e1,e2,e3)
       implicit none
       complex(8), intent(in) :: e1(:),e2(:),e3(:)
       complex(8)             :: vgggg(size(e1))
       complex(8):: se1e3,se2e3,se1e2

          se1e3=sc_(e1,e3)
          se2e3=sc_(e2,e3)
          se1e2=sc_(e1,e2)
          vgggg = (0.0d0,1.0d0)*(e2*se1e3-0.5d0*(e1*se2e3+ e3*se1e2))

       end function vgggg


      function vqg(sp,e1)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vqg(size(sp))
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          vqg = (0d0,1d0)/sqrt2*spb2_(sp,e1)

      end function vqg




      function vgq(e1,sp)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vgq(size(sp))
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          vgq = (0d0,-1d0)/sqrt2*spb2_(sp,e1)

      end function vgq





      function vbqg(sp,e1)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vbqg(size(sp))
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          vbqg = (0d0,-1d0)/sqrt2*spi2_(e1,sp)

      end function vbqg



      function vgbq(e1,sp)
      implicit none
      complex(8), intent(in) :: e1(:)
      complex(8), intent(in) :: sp(:)
      complex(8) :: vgbq(size(sp))
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

          vgbq = (0d0,1d0)/sqrt2*spi2_(e1,sp)

      end function vgbq



!       function vbqq(Dv,sp1,sp2)
!       implicit none
!       complex(8), intent(in) :: sp1(:), sp2(:)
!       integer, intent(in) ::  Dv
!       integer :: i
!       complex(8) :: vbqq(Dv)
!       complex(8) :: rr, va(Dv),sp1a(size(sp1))
!       real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0
!
!           va=(0d0,0d0)
!           vbqq=(0d0,0d0)
!
!           do i=1,Dv
!              if (i.eq.1) then
!                va(1)=(1d0,0d0)
!              else
!                va(i)=(-1d0,0d0)
!              endif
!              sp1a=spb2_(sp1,va)
!
!              rr=(0d0,-1d0)/sqrt2*psp1_(sp1a,sp2)
!              if (i.eq.1) then
!                   vbqq = vbqq + rr*va
!               else
!                   vbqq = vbqq - rr*va
!              endif
!              va(i)=(0d0,0d0)
!           enddo
!
!       end function vbqq





      function vbqq(Dummy,sp1,sp2)! this is my own simpler 4-dim version
      implicit none
      complex(8), intent(in) :: sp1(1:4), sp2(1:4)
      integer :: i,Dummy
      complex(8) :: vbqq(4)
      complex(8) :: sp1a(4)
      real(8) :: va(1:4,1:4)
      real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

         va(1,1:4)=(/+1d0,0d0,0d0,0d0/)
         va(2,1:4)=(/0d0,-1d0,0d0,0d0/)
         va(3,1:4)=(/0d0,0d0,-1d0,0d0/)
         va(4,1:4)=(/0d0,0d0,0d0,-1d0/)

          do i=1,4
             sp1a=SpiVL(sp1,dcmplx(va(i,1:4)))
             vbqq(i) = (sp1a(1)*sp2(1)+sp1a(2)*sp2(2)+sp1a(3)*sp2(3)+sp1a(4)*sp2(4)) * (0d0,-1d0)/sqrt2
          enddo


      end function vbqq




      function vvbqq(sp1,sp2)
      implicit none
      complex(8), intent(in) :: sp1(:), sp2(:)
      integer :: i,j
      complex(8) :: vvbqq(4,4)
      complex(8) :: sp1a(4)
      real(8) :: va(1:4,1:4)

         va(1,1:4)=(/+1d0,0d0,0d0,0d0/)
         va(2,1:4)=(/0d0,-1d0,0d0,0d0/)
         va(3,1:4)=(/0d0,0d0,-1d0,0d0/)
         va(4,1:4)=(/0d0,0d0,0d0,-1d0/)

         do i=1,4
          do j=1,4
             sp1a=spb2_(sp1, dcmplx(va(i,1:4)))
             sp1a=spb2_(sp1a,dcmplx(va(j,1:4)))
             vvbqq(i,j) = psp1_(sp1a,sp2)
          enddo
         enddo

      end function vvbqq







END MODULE
