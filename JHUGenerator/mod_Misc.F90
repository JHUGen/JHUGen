MODULE ModMisc
implicit none


INTERFACE swap
   MODULE PROCEDURE swapi
   MODULE PROCEDURE swapr
   MODULE PROCEDURE swapc
   MODULE PROCEDURE swap_mom
   MODULE PROCEDURE swap_cmom
END INTERFACE swap

INTERFACE OPERATOR (.dot.)
   MODULE PROCEDURE MinkowskyProduct
   MODULE PROCEDURE MinkowskyProductC
   MODULE PROCEDURE MinkowskyProductRC
   MODULE PROCEDURE MinkowskyProductCR
END INTERFACE OPERATOR (.dot.)

INTERFACE OPERATOR (.cross.)
   MODULE PROCEDURE VectorCross
END INTERFACE OPERATOR (.cross.)

interface BubleSort
   module procedure BubleSort_integerinteger,BubleSort_realinteger, BubleSort_stringlogical, BubleSort_stringinteger, BubleSort_stringreal8, BubleSort_stringcomplex8, BubleSort_stringstring
end interface BubleSort

interface ReOrder
   module procedure ReOrder_integerinteger,ReOrder_integerreal
end interface ReOrder


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
                     - p1(4)*p2(4)
END FUNCTION MinkowskyProductC

FUNCTION MinkowskyProductRC(p1,p2)
implicit none
real(8),    intent(in) :: p1(1:4)
complex(8), intent(in) :: p2(1:4)
complex(8)             :: MinkowskyProductRC

   MinkowskyProductRC = p1(1)*p2(1)  &
                      - p1(2)*p2(2)  &
                      - p1(3)*p2(3)  &
                      - p1(4)*p2(4)
END FUNCTION MinkowskyProductRC

FUNCTION MinkowskyProductCR(p1,p2)
implicit none
real(8),    intent(in) :: p2(1:4)
complex(8), intent(in) :: p1(1:4)
complex(8)             :: MinkowskyProductCR

   MinkowskyProductCR = p1(1)*p2(1)  &
                      - p1(2)*p2(2)  &
                      - p1(3)*p2(3)  &
                      - p1(4)*p2(4)
END FUNCTION MinkowskyProductCR


function LogicalToInteger(var)
implicit none
logical :: var
integer :: LogicalToInteger
   if (var) then
      LogicalToInteger = 1
   else
      LogicalToInteger = 0
   endif
end function


double complex function et1(e1,e2,e3,e4)
implicit none
complex(8), intent(in) :: e1(4), e2(4), e3(4), e4(4)
   et1 =  e1(1)*e2(2)*e3(3)*e4(4)-e1(1)*e2(2)*e3(4)*e4(3) &
   -e1(1)*e2(3)*e3(2)*e4(4)+e1(1)*e2(3)*e3(4)*e4(2) &
   +e1(1)*e2(4)*e3(2)*e4(3)-e1(1)*e2(4)*e3(3)*e4(2) &
   -e1(2)*e2(1)*e3(3)*e4(4)+e1(2)*e2(1)*e3(4)*e4(3) &
   +e1(2)*e2(3)*e3(1)*e4(4)-e1(2)*e2(3)*e3(4)*e4(1) &
   -e1(2)*e2(4)*e3(1)*e4(3)+e1(2)*e2(4)*e3(3)*e4(1) &
   +e1(3)*e2(1)*e3(2)*e4(4)-e1(3)*e2(1)*e3(4)*e4(2) &
   -e1(3)*e2(2)*e3(1)*e4(4)+e1(3)*e2(2)*e3(4)*e4(1) &
   +e1(3)*e2(4)*e3(1)*e4(2)-e1(3)*e2(4)*e3(2)*e4(1) &
   -e1(4)*e2(1)*e3(2)*e4(3)+e1(4)*e2(1)*e3(3)*e4(2) &
   +e1(4)*e2(2)*e3(1)*e4(3)-e1(4)*e2(2)*e3(3)*e4(1) &
   -e1(4)*e2(3)*e3(1)*e4(2)+e1(4)*e2(3)*e3(2)*e4(1)
   return
end function et1

double complex function sc(q1,q2)
implicit none
complex(8), intent(in) :: q1(4)
complex(8), intent(in) :: q2(4)
   sc = q1.dot.q2
   return
end function sc

double precision function scr(p1,p2)
implicit none
real(8), intent(in) :: p1(4),p2(4)
   scr = p1.dot.p2
   return
end function scr

double complex function scrc(p1,p2)
implicit none
real(8), intent(in) :: p1(4)
complex(8), intent(in) :: p2(4)
   scrc = p1.dot.p2
   return
end function scrc



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

   Get_MInv = dsqrt( dabs(Mom(1:4).dot.Mom(1:4)) )

RETURN
END FUNCTION

FUNCTION Get_MInv2(Mom)
implicit none
real(8) ::Mom(1:4),Get_MInv2

   Get_MInv2 = Mom(1:4).dot.Mom(1:4)

RETURN
END FUNCTION

FUNCTION Get_ETA(Mom)
implicit none
real(8) ::Mom(1:4),Get_ETA

   Get_ETA = 0.5d0*dlog( dabs((Mom(1)+Mom(4)+1d-16)/(Mom(1)-Mom(4) +1d-16) ))

RETURN
END FUNCTION


FUNCTION Get_PHI(Mom)
implicit none
real(8) ::Mom(1:4),Get_PHI

   Get_PHI = datan2( Mom(3),Mom(2) )

RETURN
END FUNCTION




FUNCTION Get_CosAlpha(Mom1,Mom2)
implicit none
real(8) ::Mom1(1:4),Mom2(1:4),Get_CosAlpha

    Get_CosAlpha = (Mom1(2)*Mom2(2)+Mom1(3)*Mom2(3)+Mom1(4)*Mom2(4))/dsqrt(Mom1(2)**2+Mom1(3)**2+Mom1(4)**2)/dsqrt(Mom2(2)**2+Mom2(3)**2+Mom2(4)**2)

RETURN
END FUNCTION


FUNCTION Get_CosTheta(Mom)! = Mom.nz/abs(Mom)
implicit none
real(8) ::Mom(1:4), Get_CosTheta

    Get_CosTheta = Mom(4)/dsqrt( Mom(2)**2+Mom(3)**2+Mom(4)**2 +1d-16)

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




SUBROUTINE pT_order(N,Mom)
implicit none
integer :: N
real(8) :: Mom(1:4,1:N),Mom_Tmp(1:4,1:N),pTList(1:N)
integer :: i,MomOrder(1:N)


    if(N.lt.1) return
    do i=1,N
      pTList(i) = get_PT(Mom(1:4,i))
      MomOrder(i) = i
    enddo

    call BubleSort(N,pTList(1:N),MomOrder(1:N))

    Mom_Tmp(1:4,1:N) = Mom(1:4,1:N)
    do i=1,N
        Mom(1:4,i) = Mom_Tmp(1:4,MomOrder(i))
    enddo

END SUBROUTINE



#define MakeBubleSort(typex, typey, name) \
SUBROUTINE name(N,X, IY)                 ;\
IMPLICIT NONE                            ;\
integer n                                ;\
typex :: x(1:n)                          ;\
typey :: iy(1:n)                         ;\
typex :: temp                            ;\
integer :: i, j, jmax                    ;\
typey :: itemp                           ;\
logical :: keepgoing                     ;\
                                         ;\
      keepgoing = .true.                 ;\
      jmax=n-1                           ;\
      do i=1,n-1                         ;\
         keepgoing = .false.             ;\
         do j=1,jmax                     ;\
            if(x(j).gt.x(j+1)) cycle     ;\
              keepgoing = .true.         ;\
              temp=x(j)                  ;\
              x(j)=x(j+1)                ;\
              x(j+1)=temp                ;\
              itemp=iy(j)                ;\
              iy(j)=iy(j+1)              ;\
              iy(j+1)=itemp              ;\
         enddo                           ;\
         if(.not.keepgoing) return       ;\
         jmax=jmax-1                     ;\
       enddo                             ;\
                                         ;\
RETURN                                   ;\
END SUBROUTINE

MakeBubleSort(integer(8), integer, BubleSort_integerinteger)
MakeBubleSort(real(8), integer, BubleSort_realinteger)
MakeBubleSort(character(len=100), logical, BubleSort_stringlogical)
MakeBubleSort(character(len=100), integer, BubleSort_stringinteger)
MakeBubleSort(character(len=100), real(8), BubleSort_stringreal8)
MakeBubleSort(character(len=100), complex(8), BubleSort_stringcomplex8)
MakeBubleSort(character(len=100), character(len=100), BubleSort_stringstring)

#undef MakeBubleSort

! check the routine
! real(8) :: x(1:10)
! integer :: iy(1:10)
!     x(1:10) = (/3d0,62d0,2d0,78d0,32d0,87d0,1d0,199d0,4d0,73d0/)
!     iy(1:10) = (/1,2,3,4,5,6,7,8,9,10/)
!     print *, x(:)
!     call  BubleSort(10,X(1:10), IY)
!     print *, x(:)
!     print *, iy(:)
!     stop


#define MakeReOrder(typey, typex, name)      \
SUBROUTINE name(N, iy,x)                    ;\
integer :: n,i                              ;\
typex :: x(1:n)                             ;\
typey :: iy(1:n)                            ;\
typex :: temp(1:n)                          ;\
                                            ;\
     do i=1,N                               ;\
         temp(i) = x( iy(i) )               ;\
     enddo                                  ;\
     x(1:N) = temp(1:N)                     ;\
                                            ;\
RETURN                                      ;\
END SUBROUTINE                              ;\

MakeReOrder(integer,integer,  ReOrder_integerinteger)
MakeReOrder(integer, real(8), ReOrder_integerreal)
#undef MakeReOrder


SUBROUTINE printMom(Mom)
implicit none
real(8) :: Mom(:,:)
integer :: n

    do n=1,size(Mom,2)
      write (*,"(A,I2,A,1F20.16,A,1F20.16,A,1F20.16,A,1F20.16,A)") "Mom(1:4,",n,")= (/",Mom(1,n),"d0,",Mom(2,n),"d0,",Mom(3,n),"d0,",Mom(4,n),"d0   /)"
    enddo

RETURN
END SUBROUTINE




SUBROUTINE Error(Message,ErrNum)
implicit none
character(*) :: Message
integer,optional :: ErrNum

   if( present(ErrNum) ) then
      print *, "ERROR: ",Message,ErrNum
   else
      print *, "ERROR: ",Message
   endif
   stop 1
END SUBROUTINE


FUNCTION IsNan(x)
implicit none
logical IsNan
real(8) :: x

   IsNaN = .not.(x.eq.x)

RETURN
END FUNCTION




SUBROUTINE swapi(i,j)
implicit none
integer :: i,j,temp

    temp=j
    j=i
    i=temp

END SUBROUTINE

SUBROUTINE swapr(i,j)
implicit none
real(8) :: i,j,temp

    temp=j
    j=i
    i=temp

END SUBROUTINE

SUBROUTINE swapc(i,j)
implicit none
complex(8) :: i,j,temp

    temp=j
    j=i
    i=temp

END SUBROUTINE

SUBROUTINE swap_mom(Mom1,Mom2)
implicit none
real(8) :: Mom1(1:4),Mom2(1:4),tmp(1:4)

    tmp(1:4) = Mom2(1:4)
    Mom2(1:4) = Mom1(1:4)
    Mom1(1:4) = tmp(1:4)

RETURN
END SUBROUTINE

SUBROUTINE swap_cmom(Mom1,Mom2)
implicit none
complex(8) :: Mom1(1:4),Mom2(1:4),tmp(1:4)

    tmp(1:4) = Mom2(1:4)
    Mom2(1:4) = Mom1(1:4)
    Mom1(1:4) = tmp(1:4)

RETURN
END SUBROUTINE




FUNCTION LeviCiv(e1,e2,e3,e4)
implicit none
complex(8), intent(in)  :: e1(4), e2(4), e3(4), e4(4)
complex(8)  ::LeviCiv

LeviCiv =  e1(1)*e2(2)*e3(3)*e4(4)-e1(1)*e2(2)*e3(4)*e4(3) &
          -e1(1)*e2(3)*e3(2)*e4(4)+e1(1)*e2(3)*e3(4)*e4(2) &
          +e1(1)*e2(4)*e3(2)*e4(3)-e1(1)*e2(4)*e3(3)*e4(2) &
          -e1(2)*e2(1)*e3(3)*e4(4)+e1(2)*e2(1)*e3(4)*e4(3) &
          +e1(2)*e2(3)*e3(1)*e4(4)-e1(2)*e2(3)*e3(4)*e4(1) &
          -e1(2)*e2(4)*e3(1)*e4(3)+e1(2)*e2(4)*e3(3)*e4(1) &
          +e1(3)*e2(1)*e3(2)*e4(4)-e1(3)*e2(1)*e3(4)*e4(2) &
          -e1(3)*e2(2)*e3(1)*e4(4)+e1(3)*e2(2)*e3(4)*e4(1) &
          +e1(3)*e2(4)*e3(1)*e4(2)-e1(3)*e2(4)*e3(2)*e4(1) &
          -e1(4)*e2(1)*e3(2)*e4(3)+e1(4)*e2(1)*e3(3)*e4(2) &
          +e1(4)*e2(2)*e3(1)*e4(3)-e1(4)*e2(2)*e3(3)*e4(1) &
          -e1(4)*e2(3)*e3(1)*e4(2)+e1(4)*e2(3)*e3(2)*e4(1)

END FUNCTION


FUNCTION pol_gluon_incoming(p,hel)
implicit none
complex(8) :: pol_gluon_incoming(1:4)
real(8) :: p(1:4)
integer :: hel
real(8) :: pv,ct,st,cphi,sphi


    pv=dsqrt(dabs(p(1)**2))
    ct=p(4)/pv
    st=dsqrt(abs(1.0d0-ct**2))

    if( st.le.1d-9 ) then
       cphi=1.0d0
       sphi=0.0d0
    else
       cphi= p(2)/pv/st
       sphi= p(3)/pv/st
    endif
    ! -- distinguish between positive and negative energies
!     if ( p(1) .gt. 0.0_dp) then
!        hel=hel
!     else
!        hel=-hel
!     endif

!     if (present(outgoing)) then
!        if (outgoing) pol = -pol
!     endif

    pol_gluon_incoming(1) = 0d0
    pol_gluon_incoming(2) = ct*cphi - (0d0,1d0)*hel*sphi
    pol_gluon_incoming(3) = ct*sphi + (0d0,1d0)*hel*cphi
    pol_gluon_incoming(4) = -st

    pol_gluon_incoming(2:4) = pol_gluon_incoming(2:4)/dsqrt(2.0d0)
END FUNCTION

FUNCTION pol_Zff_outgoing(pf,pfbar,hel) !  does not contain the division by 1/q^2 as in pol_dk2mom
implicit none
complex(8) :: pol_Zff_outgoing(4)
integer :: hel
real(8):: pf(1:4),pfbar(1:4)
complex(8) :: Ub(1:4),V(1:4)
    Ub(1:4) = ubar_spinor(pf,hel)
    V(1:4)  = v_spinor(pfbar,-hel)

    !   This is an expression for (-i)* (-i) Ub(+/-)) Gamma^\mu V(-/+)
    pol_Zff_outgoing(1)=-(Ub(2)*V(4)+V(2)*Ub(4)+Ub(1)*V(3)+V(1)*Ub(3))
    pol_Zff_outgoing(2)=-(-Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)+V(2)*Ub(3))
    pol_Zff_outgoing(3)=-(0d0,1d0)*(Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)-V(2)*Ub(3))
    pol_Zff_outgoing(4)=-(Ub(2)*V(4)-V(2)*Ub(4)-Ub(1)*V(3)+V(1)*Ub(3))
RETURN
END FUNCTION


! ubar spinor, massless
FUNCTION ubar_spinor(p,hel)
implicit none
real(8) :: p(1:4)
integer :: hel
complex(8) :: ubar_spinor(1:4)
complex(8) :: fc, fc2


    fc2 = p(1) + p(4)
    fc=sqrt(fc2)

    if( abs(fc2).gt.1d-9 ) then
       if(hel.eq.1) then
          ubar_spinor(1)=0d0
          ubar_spinor(2)=0d0
          ubar_spinor(3)=fc
          ubar_spinor(4)=(p(2)-(0d0,1d0)*p(3))/fc
       elseif (hel.eq.-1) then
          ubar_spinor(1)=(p(2)+(0d0,1d0)*p(3))/fc
          ubar_spinor(2)=-fc
          ubar_spinor(3)=0d0
          ubar_spinor(4)=0d0
       endif
    else
       if(hel.eq.1) then
          ubar_spinor(1) = 0d0
          ubar_spinor(2) = 0d0
          ubar_spinor(3) = 0d0
          ubar_spinor(4) = dsqrt(2*p(1))
       elseif (hel.eq.-1) then
          ubar_spinor(1) = dsqrt((2*p(1)))
          ubar_spinor(2) = 0d0
          ubar_spinor(3) = 0d0
          ubar_spinor(4) = 0d0
       endif
    endif

RETURN
END FUNCTION

! -- v0  spinor, massless
FUNCTION v_spinor(p,hel)
implicit none
real(8) :: p(1:4)
integer :: hel
complex(8) :: v_spinor(1:4)
complex(8) :: fc2, fc


    fc2 = p(1) + p(4)
    fc=sqrt(fc2)

    if( abs(fc2).gt.1d-9 ) then
       if(hel.eq.1) then
          v_spinor(1)=0d0
          v_spinor(2)=0d0
          v_spinor(3)=(p(2)-(0d0,1d0)*p(3))/fc
          v_spinor(4)=-fc
       elseif (hel.eq.-1) then
          v_spinor(1)=fc
          v_spinor(2)=(p(2)+(0d0,1d0)*p(3))/fc
          v_spinor(3)=0d0
          v_spinor(4)=0d0
       endif
    else
       if(hel.eq.1) then
          v_spinor(1)=0d0
          v_spinor(2)=0d0
          v_spinor(3)=dsqrt(2*p(1))
          v_spinor(4)=0d0
       elseif (hel.eq.-1) then
          v_spinor(1)=0d0
          v_spinor(2)=dsqrt(2*p(1))
          v_spinor(3)=0d0
          v_spinor(4)=0d0
       endif
    endif

RETURN
END FUNCTION

function Chir_Weyl(sign,sp)   ! Chir = sp.omega_sign = omega_sign.sp
implicit none
integer :: sign
double complex :: sp(1:4)
double complex :: Chir_Weyl(1:4)
   if(sign.eq.+1) then !omega_+
      Chir_Weyl(1) = sp(1)
      Chir_Weyl(2) = sp(2)
      Chir_Weyl(3) = 0d0
      Chir_Weyl(4) = 0d0
   else !omega_-
      Chir_Weyl(1) = 0d0
      Chir_Weyl(2) = 0d0
      Chir_Weyl(3) = sp(3)
      Chir_Weyl(4) = sp(4)
   endif
   return
end function

FUNCTION vbqq_Weyl(sp1,sp2)
implicit none
complex(8), intent(in) :: sp1(:), sp2(:)
integer, parameter ::  Dv=4
integer :: i
complex(8) :: vbqq_Weyl(Dv)
complex(8) :: rr, va(Dv),sp1a(4)

   va=(0d0,0d0)
   vbqq_Weyl=(0d0,0d0)

   do i=1,Dv
      if (i.eq.1) then
         va(1)=(1d0,0d0)
      else
         va(i)=(-1d0,0d0)
      endif
      sp1a=spb2_Weyl(sp1,va)


      rr=sum(sp1a(1:4)*sp2(1:4))
      if (i.eq.1) then
         vbqq_Weyl = vbqq_Weyl + rr*va
      else
         vbqq_Weyl = vbqq_Weyl - rr*va
      endif
      va(i)=(0d0,0d0)
   enddo

END FUNCTION

function spb2_Weyl(sp,v)
implicit none
complex(8), intent(in) :: sp(:),v(:)
complex(8) :: spb2_Weyl(4)
complex(8) :: x0(4,4),xx(4,4),xy(4,4)
complex(8) :: xz(4,4),x5(4,4)
complex(8) :: y1,y2,y3,y4,bp,bm,cp,cm
integer :: i,i1,i2,i3,Dv,Ds,imax

   Ds = 4
   Dv = 4
   imax = Ds/4

   do i=1,imax
      i1= 1+4*(i-1)
      i2=i1+3

      y1=sp(i1)
      y2=sp(i1+1)
      y3=sp(i1+2)
      y4=sp(i1+3)

      x0(1,i)=y3
      x0(2,i)=y4
      x0(3,i)=y1
      x0(4,i)=y2

      xx(1,i) = y4
      xx(2,i) = y3
      xx(3,i) = -y2
      xx(4,i) = -y1

      xy(1,i)=(0d0,1d0)*y4
      xy(2,i)=-(0d0,1d0)*y3
      xy(3,i)=-(0d0,1d0)*y2
      xy(4,i)=(0d0,1d0)*y1

      xz(1,i)=y3
      xz(2,i)=-y4
      xz(3,i)=-y1
      xz(4,i)=y2

      x5(1,i)=y1
      x5(2,i)=y2
      x5(3,i)=-y3
      x5(4,i)=-y4
   enddo


   do i=1,4
      spb2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)
   enddo
end function

function spi2_Weyl(v,sp)
implicit none
complex(8), intent(in) :: sp(:),v(:)
complex(8) :: spi2_Weyl(4)
complex(8) :: x0(4,4),xx(4,4),xy(4,4)
complex(8) :: xz(4,4),x5(4,4)
complex(8) ::  y1,y2,y3,y4,bp,bm,cp,cm
integer :: i,i1,i2,i3,imax,Dv,Ds

   Ds = 4
   Dv = 4

   imax = Ds/4

   do i=1,imax
      i1= 1+4*(i-1)
      i2=i1+3

      y1=sp(i1)
      y2=sp(i1+1)
      y3=sp(i1+2)
      y4=sp(i1+3)

      x0(1,i)=y3
      x0(2,i)=y4
      x0(3,i)=y1
      x0(4,i)=y2


      xx(1,i) = -y4
      xx(2,i) = -y3
      xx(3,i) = y2
      xx(4,i) = y1


      xy(1,i)=(0d0,1d0)*y4
      xy(2,i)=-(0d0,1d0)*y3
      xy(3,i)=-(0d0,1d0)*y2
      xy(4,i)=(0d0,1d0)*y1

      xz(1,i)=-y3
      xz(2,i)=y4
      xz(3,i)=y1
      xz(4,i)=-y2

      x5(1,i)=y1
      x5(2,i)=y2
      x5(3,i)=-y3
      x5(4,i)=-y4
   enddo
   do i=1,4
      spi2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1) -v(3)*xy(i,1)-v(4)*xz(i,1)
   enddo
end function




function FindInputFmt0(EventInfoLine)
implicit none
character(len=*) :: EventInfoLine
character(len=150) FindInputFmt0
integer :: i, j, fieldwidth, spaces(1:6)
integer :: ProcessIdCharacters, WeightScaleAqedAqcdCharacters(1:4), WeightScaleAqedAqcdAfterDecimal(1:4)
character(len=40) :: FormatParts(1:6)

!find the number of spaces at the beginning
spaces(1) = 0
i = 1
do while (EventInfoLine(i:i) .eq. " ")
    i = i+1
    spaces(1) = spaces(1)+1
end do
!now find the width of the number of particles, assume it's right aligned with a max of 2 digits
fieldwidth = 0
do while (EventInfoLine(i:i) .ne. " ")
    i = i+1
    fieldwidth = fieldwidth+1
end do
spaces(1) = spaces(1) + fieldwidth - 2

!spaces and process id.
!"The process ID’s are not intended to be generic" [arXiv:0109068]
!I will assume that however many digits they are, they're right aligned

ProcessIdCharacters = -1  !so that it's 0 after the first space
spaces(2) = 1             !this is not increased, it's exactly 1
do while (EventInfoLine(i:i) .eq. " ")
    i = i+1
    ProcessIdCharacters = ProcessIdCharacters+1
end do
do while (EventInfoLine(i:i) .ne. " ")
    i = i+1
    ProcessIdCharacters = ProcessIdCharacters+1
end do

!the rest of the fields (scale, alpha_QED, alpha_QCD) are decimals
do j=1,4
    spaces(j+2) = 0
    do while (EventInfoLine(i:i) .eq. " ")
        i = i+1
        spaces(j+2) = spaces(j+2)+1
    end do
    if (EventInfoLine(i:i) .eq. "-") then
        i = i+1
        WeightScaleAqedAqcdCharacters(j) = 1  !we are already past the -
    else if (spaces(j+2).eq.1) then
        WeightScaleAqedAqcdCharacters(j) = 0
    else
        spaces(j+2) = spaces(j+2)-1           !because the place where the - sign is supposed to go is not always a space
        WeightScaleAqedAqcdCharacters(j) = 1  !we are already past the -
    endif
    do while (EventInfoLine(i:i) .ne. " ")
        i = i+1
        WeightScaleAqedAqcdCharacters(j) = WeightScaleAqedAqcdCharacters(j)+1
    end do
    WeightScaleAqedAqcdAfterDecimal(j) = WeightScaleAqedAqcdCharacters(j)-7
end do

!now we construct the format string
if (spaces(1).eq.0) then
    FormatParts(1) = "(I2,"
else
    write(FormatParts(1), "(A,I1,A)") "(", spaces(1), "X,I2,"
endif
if (ProcessIdCharacters.lt.10) then
    write(FormatParts(2), "(I1,A,I1,A)") spaces(2), "X,I", ProcessIdCharacters !the comma is at the beginning of the next part
else
    write(FormatParts(2), "(I1,A,I2)") spaces(2), "X,I", ProcessIdCharacters !the comma is at the beginning of the next part
endif
do j=1,4
    if (WeightScaleAqedAqcdCharacters(j).lt.10) then
        write(FormatParts(j+2), "(A,I1,A,I1,A,I1)") ",", spaces(j+2), "X,1PE", WeightScaleAqedAqcdCharacters(j), ".", WeightScaleAqedAqcdAfterDecimal(j)
    elseif (WeightScaleAqedAqcdCharacters(j).lt.17) then
        write(FormatParts(j+2), "(A,I1,A,I2,A,I1)") ",", spaces(j+2), "X,1PE", WeightScaleAqedAqcdCharacters(j), ".", WeightScaleAqedAqcdAfterDecimal(j)
    else
        write(FormatParts(j+2), "(A,I1,A,I2,A,I2)") ",", spaces(j+2), "X,1PE", WeightScaleAqedAqcdCharacters(j), ".", WeightScaleAqedAqcdAfterDecimal(j)
    endif
end do

FindInputFmt0 = (trim(FormatParts(1))  &
              // trim(FormatParts(2))  &
              // trim(FormatParts(3))  &
              // trim(FormatParts(4))  &
              // trim(FormatParts(5))  &
              // trim(FormatParts(6))  &
              // ")")

return
end function FindInputFmt0



function FindInputFmt1(ParticleLine)
implicit none
character(len=*) :: ParticleLine
integer :: i, j, fieldwidth, spaces(1:13)
character(len=150) :: FindInputFmt1
integer :: MomentumCharacters(1:5), LifetimeCharacters, LifetimeDigitsAfterDecimal, SpinCharacters, SpinDigitsAfterDecimal
logical :: LifetimeIsExponential, SpinIsExponential
character(len=40) :: FormatParts(11)
!first find the number of spaces at the beginning
i = 1
spaces(1) = 0
do while (ParticleLine(i:i) .eq. " ")
    i = i+1
    spaces(1) = spaces(1)+1
end do
!now find the width of the id; assume it's right aligned.
fieldwidth = 0
do while (ParticleLine(i:i) .ne. " ")
    i = i+1
    fieldwidth = fieldwidth+1
end do
spaces(1) = spaces(1) + fieldwidth - 3

!number of spaces between the id and the status
spaces(2) = 0
do while (ParticleLine(i:i) .eq. " ")
    i = i+1
    spaces(2) = spaces(2)+1
end do
!width of the status (should be -1, so width 2, but just in case)
fieldwidth = 0
do while (ParticleLine(i:i) .ne. " ")
    i = i+1
    fieldwidth = fieldwidth+1
end do
spaces(2) = spaces(2) + fieldwidth - 2

!number of spaces between the status and the first mother
spaces(3) = 0
do while (ParticleLine(i:i) .eq. " ")
    i = i+1
    spaces(3) = spaces(3)+1
end do
fieldwidth = 0
do while (ParticleLine(i:i) .ne. " ")
    i = i+1
    fieldwidth = fieldwidth+1
end do
spaces(3) = spaces(3) + fieldwidth - 2

!number of spaces between the mothers
spaces(4) = 0
do while (ParticleLine(i:i) .eq. " ")
    i = i+1
    spaces(4) = spaces(4)+1
end do
fieldwidth = 0
do while (ParticleLine(i:i) .ne. " ")
    i = i+1
    fieldwidth = fieldwidth+1
end do
spaces(4) = spaces(4) + fieldwidth - 2

!number of spaces between the second mother and the color
spaces(5) = 0
do while (ParticleLine(i:i) .eq. " ")
    i = i+1
    spaces(5) = spaces(5)+1
end do
fieldwidth = 0
do while (ParticleLine(i:i) .ne. " ")
    i = i+1
    fieldwidth = fieldwidth+1
end do
spaces(5) = spaces(5) + fieldwidth - 3

!number of spaces between the color and the anticolor
spaces(6) = 0
do while (ParticleLine(i:i) .eq. " ")
    i = i+1
    spaces(6) = spaces(6)+1
end do
fieldwidth = 0
do while (ParticleLine(i:i) .ne. " ")
    i = i+1
    fieldwidth = fieldwidth+1
end do
spaces(6) = spaces(6) + fieldwidth - 3

!From now on the alignment is simpler, every row will have the same width
!   except for possibly a - sign at the beginning
!Unfortunately now there's number formatting to worry about
do j=7,11   !px, py, pz, E, m
    !spaces before the momentum component
    spaces(j) = 0
    do while (ParticleLine(i:i) .eq. " ")
        i = i+1
        spaces(j) = spaces(j)+1
    end do
    if (ParticleLine(i:i) .eq. "-") then
        i = i+1
    else
        spaces(j) = spaces(j)-1  !because the place where the - sign is supposed to go is not always a space
    endif
    !i is now on the first actual digit (not -) of the component

    !number of characters used for the component
    MomentumCharacters(j-6) = 1           !we are already past the -
    do while (ParticleLine(i:i) .ne. " ")
        i = i+1
        MomentumCharacters(j-6) = MomentumCharacters(j-6)+1
    end do
enddo

!number of spaces between m and lifetime
!lifetime is nonnegative, so no - sign
!but some generators (JHUGen) write 0.00000000000E+00
!while some (old JHUGen, some versions of MadGraph) write 0.
spaces(12) = 0
do while (ParticleLine(i:i) .eq. " ")
    i = i+1
    spaces(12) = spaces(12)+1
end do

!lifetime formatting
LifetimeCharacters = 0
LifetimeDigitsAfterDecimal = -1
LifetimeIsExponential = .false.
do while (ParticleLine(i:i) .ne. " ")
    if ((LifetimeDigitsAfterDecimal .ge. 0 .and. .not.LifetimeIsExponential) &
            .or. ParticleLine(i:i) .eq. ".") then
        LifetimeDigitsAfterDecimal = LifetimeDigitsAfterDecimal+1
    endif
    if (ParticleLine(i:i) .eq. "E" .or. ParticleLine(i:i) .eq. "e") then
        LifetimeIsExponential = .true.
    endif
    i = i+1
    LifetimeCharacters = LifetimeCharacters+1
enddo

!number of spaces between lifetime and spin
!spin has all the complications of lifetime, plus it might be negative
spaces(13) = 0
do while (ParticleLine(i:i) .eq. " ")
    i = i+1
    spaces(13) = spaces(13)+1
end do
if (ParticleLine(i:i) .eq. "-") then
    i = i+1
else
    spaces(13) = spaces(13)-1  !because the place where the - sign is supposed to go is not always a space
endif
!i is now on the first digit of the spin

!spin formatting
SpinCharacters = 1             !we are already past the -
SpinDigitsAfterDecimal = -1
SpinIsExponential = .false.
do while (ParticleLine(i:i) .ne. " ")
    if (SpinDigitsAfterDecimal .ge. 0 .or. ParticleLine(i:i) .eq. ".") then
        SpinDigitsAfterDecimal = SpinDigitsAfterDecimal+1
    endif
    if (ParticleLine(i:i) .eq. "E" .or. ParticleLine(i:i) .eq. "e") then
        SpinIsExponential = .true.
    endif
    i = i+1
    SpinCharacters = SpinCharacters+1
enddo

!check that spaces(i) > 0
!(besides for spaces(1), before the id, which can be 0)
!the only way this can happen, consistent with the definition of lhe format,
! is for a column not to leave extra space for a -
do i=2,13
    if (spaces(i).eq.0) then
        spaces(i) = 1
        if (i.ge.7 .and. i.le.11) then   !we counted the nonexistant space for - in MomentumCharacters(i-6)
            MomentumCharacters(i-6) = MomentumCharacters(i-6)-1
        endif
        if (i.eq.13) then  !same
            SpinCharacters = SpinCharacters-1
        endif
    endif
enddo

!Ok, now we construct the format string
!Initial spaces and id
if (spaces(1).eq.0) then
    FormatParts(1) = "(I3,"
else
    write(FormatParts(1),"(A,I1,A)") "(", spaces(1), "X,I3,"
endif
!spaces and status
write(FormatParts(2),"(I1,A)") spaces(2), "X,I2,"
!spaces and mothers
write(FormatParts(3),"(I1,A,I1,A)") spaces(3), "X,I2,", spaces(4), "X,I2,"
!spaces and colors
write(FormatParts(4),"(I1,A,I1,A)") spaces(5), "X,I3,", spaces(6), "X,I3,"

do i=5,9
    !spaces and momentum components (and mass)
    if (MomentumCharacters(i-4) .lt. 10) then
        write(FormatParts(i),"(I1,A,I1,A,I1,A)") spaces(i+2), "X,1PE", MomentumCharacters(i-4), ".", MomentumCharacters(i-4)-7, ","
    elseif (MomentumCharacters(i-4) .lt. 17) then
        write(FormatParts(i),"(I1,A,I2,A,I1,A)") spaces(i+2), "X,1PE", MomentumCharacters(i-4), ".", MomentumCharacters(i-4)-7, ","
    else
        write(FormatParts(i),"(I1,A,I2,A,I2,A)") spaces(i+2), "X,1PE", MomentumCharacters(i-4), ".", MomentumCharacters(i-4)-7, ","
    endif
enddo

!spaces and lifetime
if (LifetimeIsExponential) then
    if (LifetimeCharacters .lt. 10) then
        write(FormatParts(10),"(I1,A,I1,A,I1,A)") &
            spaces(12), "X,1PE",LifetimeCharacters,".",LifetimeDigitsAfterDecimal,","
    elseif (LifetimeDigitsAfterDecimal.lt.10) then
        write(FormatParts(10),"(I1,A,I2,A,I1,A)") &
            spaces(12), "X,1PE",LifetimeCharacters,".",LifetimeDigitsAfterDecimal,","
    else
        write(FormatParts(10),"(I1,A,I2,A,I2,A)") &
            spaces(12), "X,1PE",LifetimeCharacters,".",LifetimeDigitsAfterDecimal,","
    endif
else
    if (LifetimeCharacters .lt. 10) then
        write(FormatParts(10),"(I1,A,I1,A,I1,A)") &
            spaces(12), "X,1F",LifetimeCharacters,".",LifetimeDigitsAfterDecimal,","
    elseif (LifetimeDigitsAfterDecimal.lt.10) then
        write(FormatParts(10),"(I1,A,I2,A,I1,A)") &
            spaces(12), "X,1F",LifetimeCharacters,".",LifetimeDigitsAfterDecimal,","
    else
        write(FormatParts(10),"(I1,A,I2,A,I2,A)") &
            spaces(12), "X,1F",LifetimeCharacters,".",LifetimeDigitsAfterDecimal,","
    endif
endif

!spaces and spin
if (SpinIsExponential) then
    if (SpinCharacters .lt. 10) then
        write(FormatParts(11),"(I1,A,I1,A,I1,A)") &
            spaces(13), "X,1PE",SpinCharacters,".",SpinDigitsAfterDecimal,")"
    elseif (SpinDigitsAfterDecimal.lt.10) then
        write(FormatParts(11),"(I1,A,I2,A,I1,A)") &
            spaces(13), "X,1PE",SpinCharacters,".",SpinDigitsAfterDecimal,")"
    else
        write(FormatParts(11),"(I1,A,I2,A,I2,A)") &
            spaces(13), "X,1PE",SpinCharacters,".",SpinDigitsAfterDecimal,")"
    endif
else
    if (SpinCharacters .lt. 10) then
        write(FormatParts(11),"(I1,A,I1,A,I1,A)") &
            spaces(13), "X,1F",SpinCharacters,".",SpinDigitsAfterDecimal,")"
    elseif (SpinDigitsAfterDecimal.lt.10) then
        write(FormatParts(11),"(I1,A,I2,A,I1,A)") &
            spaces(13), "X,1F",SpinCharacters,".",SpinDigitsAfterDecimal,")"
    else
        write(FormatParts(11),"(I1,A,I2,A,I2,A)") &
            spaces(13), "X,1F",SpinCharacters,".",SpinDigitsAfterDecimal,")"
    endif
endif

!do i=1,11
!    print *, FormatParts(i)
!enddo

FindInputFmt1 = (trim(FormatParts(1))  &
              // trim(FormatParts(2))  &
              // trim(FormatParts(3))  &
              // trim(FormatParts(4))  &
              // trim(FormatParts(5))  &
              // trim(FormatParts(6))  &
              // trim(FormatParts(7))  &
              // trim(FormatParts(8))  &
              // trim(FormatParts(9))  &
              // trim(FormatParts(10)) &
              // trim(FormatParts(11)))
return
end function FindInputFmt1



function Capitalize(InputString)
implicit none
character(len=150) :: InputString
character(len=150) :: Capitalize
character(len=26), parameter :: lower = "abcdefghijklmnopqrstuvwxyz"
character(len=26), parameter :: upper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
integer :: i, j

Capitalize = ""

do i=1,len(trim(InputString))
    j = Index(lower,InputString(i:i))
    if (j.ne.0) then
        Capitalize(i:i) = upper(j:j)
    else
        Capitalize(i:i) = InputString(i:i)
    endif
enddo
return
end function Capitalize



subroutine HouseOfRepresentatives(XSecArray, NumberOfSeats, TotalNumberOfSeats)
!Want to generate TotalNumberOfSeats events among the partonic channels, proportional to the xsec.
!But can't generate half an event, so need to round in some way.
!Apparently this is a major political issue.  https://en.wikipedia.org/wiki/Apportionment_paradox
!The first US presidential veto was when George Washington vetoed a change in rounding scheme
! because Virginia would have lost seats in the House.
!Also, if the US had used a different system, Rutherford B. Hayes wouldn't have been president.

!For our case, it's much easier because we can and should distribute randomly, not deterministically
! so that if you run 1000 times with 50 events each the fluctuations average out.
!Also states with small xsec like Wyoming don't get a seat.
implicit none
real(8) :: totalxsec, CrossSecNormalized(-5:5,-5:5), yRnd
integer :: n, i, j
integer, intent(in) :: TotalNumberOfSeats
real(8), intent(in) :: XSecArray(-5:5,-5:5)
integer(8), intent(out) :: NumberOfSeats(-5:5,-5:5)
  totalxsec = 0
  NumberOfSeats(:,:) = 0
  if (TotalNumberOfSeats.lt.1) return

  do i=-5,5
    do j=-5,5
      totalxsec = totalxsec+abs(XSecArray(i,j))
    enddo
  enddo
  if (totalxsec.eq.0d0) then
    print *,"HouseOfRepresentatives: Total xsec is 0. Aborting..."
    stop 1
  endif
  do i=-5,5
    do j=-5,5
      CrossSecNormalized(i,j) = abs(XSecArray(i,j)) / totalxsec
    enddo
  enddo

  do n=1,TotalNumberOfSeats
    do while(.true.)
      call random_number(yRnd)
      do i=-5,5
        do j=-5,5
          if (yRnd .lt. CrossSecNormalized(i,j)) then
            NumberOfSeats(i,j) = NumberOfSeats(i,j)+1
            goto 99
          endif
          yRnd = yRnd - CrossSecNormalized(i,j)
        enddo
      enddo
      !in case a rounding error causes sum(CrossSecNormalized) != 1
      print *, "HouseOfRepresentatives: Warning, rounding error, try again"
    enddo
99  continue
  enddo

  if (sum(NumberOfSeats(:,:)).ne.TotalNumberOfSeats) then
    print *, "HouseOfRepresentatives: Wrong total number of events, shouldn't be able to happen"
    stop 1
  endif
end subroutine




subroutine HouseOfRepresentatives2(XSecArray, NumberOfSeats, TotalNumberOfSeats)
!same as above with one-dim. arrays
implicit none
real(8), intent(in) :: XSecArray(:)
real(8) :: totalxsec, CrossSecNormalized(size(XSecArray)), yRnd
integer :: n,i,nchannels
integer, intent(in) :: TotalNumberOfSeats
integer(8), intent(out) :: NumberOfSeats(:)
  totalxsec = 0
  NumberOfSeats(:) = 0
  if (TotalNumberOfSeats.lt.1) return

  nchannels = size(XSecArray)

  do i=1,nchannels
    totalxsec = totalxsec+abs(XSecArray(i))
  enddo
  if (totalxsec.eq.0d0) then
    print *,"HouseOfRepresentatives2: Total xsec is 0. Aborting..."
    stop 1
  endif
  do i=1,nchannels
    CrossSecNormalized(i) = abs(XSecArray(i)) / totalxsec
  enddo

  do n=1,TotalNumberOfSeats
    do while(.true.)
      call random_number(yRnd)
      do i=1,nchannels
          if (yRnd .lt. CrossSecNormalized(i)) then
            NumberOfSeats(i) = NumberOfSeats(i)+1
            goto 99
          endif
          yRnd = yRnd - CrossSecNormalized(i)
      enddo
      print *, "HouseOfRepresentatives2: Warning, rounding error, try again"
    enddo
99  continue
  enddo

  if (sum(NumberOfSeats(:)).ne.TotalNumberOfSeats) then
    print *, "HouseOfRepresentatives2: Wrong total number of events, shouldn't be able to happen"
    stop 1
  endif
end subroutine




!========================================================================

    subroutine convert_to_MCFM(p,pout)
      implicit none
! converts from (E,px,py,pz) to (px,py,pz,E)
      real(8) :: p(1:4),tmp(1:4)
      real(8), optional :: pout(1:4)

      if( present(pout) ) then
          pout(1)=p(2)
          pout(2)=p(3)
          pout(3)=p(4)
          pout(4)=p(1)
      else
          tmp(1)=p(1)
          tmp(2)=p(2)
          tmp(3)=p(3)
          tmp(4)=p(4)

          p(1)=tmp(2)
          p(2)=tmp(3)
          p(3)=tmp(4)
          p(4)=tmp(1)
      endif

    end subroutine convert_to_MCFM

!========================================================================
!borrowed from Passarino's file, some internal variable names or comments
!might not be correct in this context
!========================================================================

SUBROUTINE EvaluateSpline(EvalPoint, SplineData, SplineDataLength, TheResult)
!SplineData: SplineDataLength by 2 array
!    SplineData(1:SplineDataLength,1) are the x values
!    SplineData(1:SplineDataLength,2) are the corresponding y values

IMPLICIT NONE

INTEGER i,top,gdim,SplineDataLength
REAL(8) u,value,EvalPoint
REAL(8), intent(out) :: TheResult
REAL(8), dimension(SplineDataLength) :: bc,cc,dc
REAL(8) :: SplineData(1:SplineDataLength, 1:2)

! u value of M_H at which the spline is to be evaluated

gdim= SplineDataLength

CALL HTO_FMMsplineSingleHt(bc,cc,dc,top,gdim,SplineData(1:SplineDataLength,1),SplineData(1:SplineDataLength,2))

u= EvalPoint
CALL HTO_Seval3SingleHt(u,bc,cc,dc,top,gdim,value,xc=SplineData(1:SplineDataLength,1),yc=SplineData(1:SplineDataLength,2))

TheResult= value

RETURN

!-----------------------------------------------------------------------

CONTAINS

SUBROUTINE HTO_FMMsplineSingleHt(b,c,d,top,gdim,xc,yc)

!---------------------------------------------------------------------------

INTEGER k,n,i,top,gdim,l

REAL(8), dimension(gdim) :: xc,yc
REAL(8), dimension(gdim) :: x,y

REAL(8), DIMENSION(gdim) :: b
! linear coeff

REAL(8), DIMENSION(gdim) :: c
! quadratic coeff.

REAL(8), DIMENSION(gdim) :: d
! cubic coeff.

REAL(8) :: t
REAL(8),PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0

! The grid


n= gdim
FORALL(l=1:gdim)
x(l)= xc(l)
y(l)= yc(l)
ENDFORALL

!.....Set up tridiagonal system.........................................
!     b=diagonal, d=offdiagonal, c=right-hand side

d(1)= x(2)-x(1)
c(2)= (y(2)-y(1))/d(1)
DO k= 2,n-1
d(k)= x(k+1)-x(k)
b(k)= TWO*(d(k-1)+d(k))
c(k+1)= (y(k+1)-y(k))/d(k)
c(k)= c(k+1)-c(k)
END DO

!.....End conditions.  third derivatives at x(1) and x(n) obtained
!     from divided differences.......................................

b(1)= -d(1)
b(n)= -d(n-1)
c(1)= ZERO
c(n)= ZERO
IF (n > 3) THEN
c(1)= c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
c(n)= c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
c(1)= c(1)*d(1)*d(1)/(x(4)-x(1))
c(n)= -c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))
END IF

DO k=2,n    ! forward elimination
t= d(k-1)/b(k-1)
b(k)= b(k)-t*d(k-1)
c(k)= c(k)-t*c(k-1)
END DO

c(n)= c(n)/b(n)

! back substitution ( makes c the sigma of text)

DO k=n-1,1,-1
c(k)= (c(k)-d(k)*c(k+1))/b(k)
END DO

!.....Compute polynomial coefficients...................................

b(n)= (y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
DO k=1,n-1
b(k)= (y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
d(k)= (c(k+1)-c(k))/d(k)
c(k)= THREE*c(k)
END DO
c(n)= THREE*c(n)
d(n)= d(n-1)

RETURN

END SUBROUTINE HTO_FMMsplineSingleHt

!------------------------------------------------------------------------

SUBROUTINE HTO_Seval3SingleHt(u,b,c,d,top,gdim,f,fp,fpp,fppp,xc,yc)

! ---------------------------------------------------------------------------

REAL(8),INTENT(IN) :: u
! abscissa at which the spline is to be evaluated

INTEGER j,k,n,l,top,gdim

REAL(8), dimension(gdim) :: xc,yc
REAL(8), dimension(gdim) :: x,y
REAL(8), DIMENSION(gdim) :: b,c,d
! linear,quadratic,cubic coeff

REAL(8),INTENT(OUT),OPTIONAL:: f,fp,fpp,fppp
! function, 1st,2nd,3rd deriv

INTEGER, SAVE :: i=1
REAL(8)    :: dx
REAL(8),PARAMETER:: TWO=2.0, THREE=3.0, SIX=6.0

! The grid

n= gdim
FORALL(l=1:gdim)
x(l)= xc(l)
y(l)= yc(l)
ENDFORALL

!.....First check if u is in the same interval found on the
!     last call to Seval.............................................

IF (  (i<1) .OR. (i >= n) ) i=1
IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
i=1

! binary search

j= n+1
DO
k= (i+j)/2
IF (u < x(k)) THEN
j= k
ELSE
i= k
ENDIF
IF (j <= i+1) EXIT
ENDDO
ENDIF

dx= u-x(i)

! evaluate the spline

IF (Present(f))    f= y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
IF (Present(fp))   fp= b(i)+dx*(TWO*c(i) + dx*THREE*d(i))
IF (Present(fpp))  fpp= TWO*c(i) + dx*SIX*d(i)
IF (Present(fppp)) fppp= SIX*d(i)

RETURN

END SUBROUTINE HTO_Seval3SingleHt

END SUBROUTINE EvaluateSpline

function CenterWithStars(string, totallength, align, padleft, padright)
implicit none
character(len=*) :: string
integer :: totallength, nspaces, nleftspaces, nrightspaces
integer, optional :: align, padleft, padright
integer :: align2, padleft2, padright2
character(len=totallength) CenterWithStars

    align2 = 2
    padleft2 = 0
    padright2 = 0
    if (present(align)) align2 = align
    if (present(padleft)) padleft2 = padleft
    if (present(padright)) padright2 = padright

    if (len(trim(string))+padleft2+padright2 .gt. totallength-2) then
        call Error("len(trim(string))+padleft+padright > totallength-2!")
    endif

    nspaces = totallength - len(trim(string)) - 2

    if (align2.eq.1) then !left
      nleftspaces = padleft2
      nrightspaces = nspaces-nleftspaces
    elseif (align2.eq.2) then !center
      nleftspaces = (nspaces+padleft2-padright2)/2
      nrightspaces = nspaces-nleftspaces
    elseif (align2.eq.3) then !right
      nrightspaces = padright2
      nleftspaces = nspaces-nrightspaces
    else
      print *, "Unknown align value", align2
    endif

    CenterWithStars = "*" // repeat(" ", nleftspaces) // trim(string) // repeat(" ", nrightspaces) // "*"

end function

SUBROUTINE PrintLogo(TheUnit, title)
use modParameters
implicit none
integer :: TheUnit
integer, parameter :: linelength = 89
character(len=*) :: title

    write(TheUnit, *) " "
    write(TheUnit, *) " ", repeat("*", linelength)
    write(TheUnit, *) " ", CenterWithStars(trim(title)//" "//trim(JHUGen_Version), linelength)
    write(TheUnit, *) " ", repeat("*", linelength)
    write(TheUnit, *) " ", CenterWithStars("", linelength)
    write(TheUnit, *) " ", CenterWithStars("Spin and parity determination of single-produced resonances at hadron colliders", linelength)
    write(TheUnit, *) " ", CenterWithStars("", linelength)
    write(TheUnit, *) " ", CenterWithStars("I. Anderson, S. Bolognesi, F. Caola, J. Davis, Y. Gao, A. V. Gritsan,", linelength)
    write(TheUnit, *) " ", CenterWithStars("L. S. Mandacaru Guerra, Z. Guo, C. B. Martin, T. Martini, K. Melnikov, R. Pan,", linelength)
    write(TheUnit, *) " ", CenterWithStars("R. Rontsch, J. Roskes, U. Sarica, M. Schulze, N. V. Tran, A. Whitbeck, M. Xiao, Y. Zhou", linelength)
    write(TheUnit, *) " ", CenterWithStars("Phys.Rev. D81 (2010) 075022;  arXiv:1001.3396  [hep-ph],", linelength)
    write(TheUnit, *) " ", CenterWithStars("Phys.Rev. D86 (2012) 095031;  arXiv:1208.4018  [hep-ph],", linelength)
    write(TheUnit, *) " ", CenterWithStars("Phys.Rev. D89 (2014) 035007;  arXiv:1309.4819  [hep-ph],", linelength)
    write(TheUnit, *) " ", CenterWithStars("Phys.Rev. D94 (2016) 055023;  arXiv:1606.03107 [hep-ph],", linelength)
    write(TheUnit, *) " ", CenterWithStars("Phys.Rev. D102 (2020) 056022; arXiv:2002.09888 [hep-ph],", linelength)
    write(TheUnit, *) " ", CenterWithStars("Phys.Rev. D102 (2021) 055045; arXiv:2104.04277 [hep-ph].", linelength)
    write(TheUnit, *) " ", CenterWithStars("                              arXiv:2109.13363 [hep-ph].", linelength)
    write(TheUnit, *) " ", CenterWithStars("", linelength)
    write(TheUnit, *) " ", repeat("*", linelength)
    write(TheUnit, *) " "
return
END SUBROUTINE


real function infinity()
  implicit none
  real :: x
  x = 0
  infinity=-log(x)
end function infinity

function CalculatesXsec(Process)
  implicit none
  integer :: Process
  logical :: CalculatesXsec
  CalculatesXsec=.false.
  if (Process.le.2) then
    CalculatesXsec=.true.
  elseif (Process.eq.50) then
    CalculatesXsec=.false.
  elseif (Process.ge.51 .and. Process.le.52) then
    CalculatesXsec=.false.
  elseif (Process.eq.60) then
    CalculatesXsec=.true.
  elseif (Process.eq.61) then
    CalculatesXsec=.true.
  elseif (Process.eq.62) then
    CalculatesXsec=.false.
  elseif (Process.ge.66 .and. Process.le.75) then
    CalculatesXsec=.true.
  elseif (Process.eq.80 .or. Process.eq.90) then
    CalculatesXsec=.false.
  elseif (Process.ge.110 .and. Process.le.114) then
    CalculatesXsec=.true.
  elseif (Process.ge.115 .and. Process.le.117) then
    CalculatesXsec=.true.
  else
    print *, "Unknown process in CalculatesXsec", process, "; setting CalculatesXsec=false."
  endif
end function CalculatesXsec

END MODULE

