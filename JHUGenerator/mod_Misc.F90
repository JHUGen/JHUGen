MODULE ModMisc
implicit none


INTERFACE OPERATOR (.dot.)
   MODULE PROCEDURE MinkowskyProduct
   MODULE PROCEDURE MinkowskyProductC
END INTERFACE OPERATOR (.dot.)

INTERFACE OPERATOR (.cross.)
   MODULE PROCEDURE VectorCross
END INTERFACE OPERATOR (.cross.)


contains



FUNCTION VectorCross(p1,p2)
implicit none
real(8), intent(in) :: p1(1:3),p2(1:3)
real(8)             :: VectorCross(1:3)

    VectorCross(1)  = p1(2)*p2(3) - p1(3)*p2(2)
    VectorCross(2)  = p1(3)*p2(1) - p1(1)*p2(3)
    VectorCross(3)  = p1(1)*p2(2) - p1(2)*p2(1)
END FUNCTION VectorCross



FUNCTION MinkowskyProduct(p1,p2)
implicit none
real(8), intent(in) :: p1(1:4),p2(1:4)
real(8)             :: MinkowskyProduct

     MinkowskyProduct = p1(1)*p2(1)  &
                      - p1(2)*p2(2)  &
                      - p1(3)*p2(3)  &
                       -p1(4)*p2(4)
END FUNCTION MinkowskyProduct


FUNCTION MinkowskyProductC(p1,p2)
implicit none
complex(8), intent(in) :: p1(1:4),p2(1:4)
complex(8)             :: MinkowskyProductC

     MinkowskyProductC = p1(1)*p2(1)  &
                      - p1(2)*p2(2)  &
                      - p1(3)*p2(3)  &
                       -p1(4)*p2(4)
END FUNCTION MinkowskyProductC


FUNCTION Get_PT(Mom)
implicit none
real(8) ::Mom(1:4),Get_PT

   Get_PT = dsqrt( Mom(2)**2 + Mom(3)**2 )

RETURN
END FUNCTION


FUNCTION Get_PT2(Mom)
implicit none
real(8) ::Mom(1:4),Get_PT2

   Get_PT2 = Mom(2)**2 + Mom(3)**2

RETURN
END FUNCTION


FUNCTION Get_MInv(Mom)
implicit none
real(8) ::Mom(1:4),Get_MInv

   Get_MInv = dsqrt( Mom(1:4).dot.Mom(1:4) )

RETURN
END FUNCTION

FUNCTION Get_ETA(Mom)
implicit none
real(8) ::Mom(1:4),Get_ETA

   Get_ETA = 0.5d0*dlog( (Mom(1)+Mom(4))/(Mom(1)-Mom(4)) )

RETURN
END FUNCTION


FUNCTION Get_PHI(Mom)
implicit none
real(8) ::Mom(1:4),Get_PHI

   Get_PHI = datan2( Mom(3),Mom(2) )

RETURN
END FUNCTION



FUNCTION Get_R(Mom1,Mom2)
implicit none
real(8),parameter :: Pi=3.141592653589793d0
real(8) :: Mom1(1:4),Mom2(1:4),Get_R
real(8) :: eta1,eta2,phi1,phi2,DeltaPhi,r2,delphi

   eta1 = Get_ETA(Mom1(1:4))
   eta2 = Get_ETA(Mom2(1:4))
   phi1 = Get_PHI(Mom1(1:4))
   phi2 = Get_PHI(Mom2(1:4))
   DeltaPhi = dabs(phi1-phi2)
   if( DeltaPhi.gt.Pi ) DeltaPhi=2d0*Pi-DeltaPhi
   Get_R = dsqrt((eta1-eta2)**2 + DeltaPhi**2)

RETURN
END FUNCTION



SUBROUTINE Error(Message,ErrNum)
implicit none
character(*) :: Message
integer,optional :: ErrNum

   if( present(ErrNum) ) then
      print *, "ERROR: ",Message,ErrNum
   else
      print *, "ERROR: ",Message
   endif
   stop
END SUBROUTINE


FUNCTION IsNan(x)
implicit none
logical IsNan
real(8) :: x

   if( .not.x.le.0d0 .and. .not.x.gt.0d0 ) then
       IsNaN=.true.
   else
      IsNaN=.false.
   endif
END FUNCTION





END MODULE

