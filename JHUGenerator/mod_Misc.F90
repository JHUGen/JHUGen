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

   Get_ETA = 0.5d0*dlog( (Mom(1)+Mom(4))/(Mom(1)-Mom(4)) )

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

    Get_CosTheta = Mom(4)/dsqrt( Mom(2)**2+Mom(3)**2+Mom(4)**2 )

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



SUBROUTINE swap_mom(Mom1,Mom2)
implicit none
real(8) :: Mom1(1:4),Mom2(1:4),tmp(1:4)

    tmp(1:4) = Mom2(1:4)
    Mom2(1:4) = Mom1(1:4)
    Mom1(1:4) = tmp(1:4)

RETURN
END SUBROUTINE

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




SUBROUTINE BubleSort(N,X, IY)
IMPLICIT NONE
integer n
real(8) x(1:n)
integer iy(1:n)
real(8) temp
integer i, j, jmax, itemp

      jmax=n-1
      do i=1,n-1
         temp=1d38
         do j=1,jmax
            if(x(j).gt.x(j+1)) cycle
              temp=x(j)
              x(j)=x(j+1)
              x(j+1)=temp
              itemp=iy(j)
              iy(j)=iy(j+1)
              iy(j+1)=itemp
         enddo
         if(temp.eq.1d38) return
         jmax=jmax-1
       enddo

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

return
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
RETURN
END FUNCTION




SUBROUTINE swapi(i,j)
implicit none
integer :: i,j,temp

    temp=j
    j=i
    i=temp

END SUBROUTINE


function FindInputFmt0(EventInfoLine)
implicit none
character(len=*) :: EventInfoLine
character(len=150) FindInputFmt0
integer :: i
i = 1
do while (EventInfoLine(i+1:i+1) .eq. " ")
    i = i+1
end do
if (i.eq.1) then
    FindInputFmt0 = "(I2,A160)"
else
    write(FindInputFmt0, "(A,I2,A)") "(", i, "X,I2,A160)"
endif
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
    !i is now on the first actual digit (not -) of px

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







subroutine spinoru(N,p,za,zb,s)
!---Calculate spinor products      
!---taken from MCFM & modified by R. Rontsch, May 2015
!---extended to deal with negative energies ie with all momenta outgoing                                                                
!---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl,                                                                                  
!---za(i,j)*zb(j,i)=s(i,j)                      
      implicit none
      real(8) :: p(:,:),two
      integer, parameter :: mxpart=14
      complex(8):: c23(N),f(N),rt(N),za(:,:),zb(:,:),czero,cone,ci
      real(8)   :: s(:,:)
      integer i,j,N
      
      if (size(p,1) .ne. N) then
         call Error("spinorz: momentum mismatch")
      endif
      two=2d0
      czero=dcmplx(0d0,0d0)
      cone=dcmplx(1d0,0d0)
      ci=dcmplx(0d0,1d0)
      

!---if one of the vectors happens to be zero this routine fails.                                                                                                                
      do j=1,N
         za(j,j)=czero
         zb(j,j)=za(j,j)

!-----positive energy case                                                                                                                                                      
         if (p(j,4) .gt. 0d0) then
            rt(j)=dsqrt(p(j,4)+p(j,1))
            c23(j)=dcmplx(p(j,3),-p(j,2))
            f(j)=cone
         else
!-----negative energy case                                                                                                                                                      
            rt(j)=dsqrt(-p(j,4)-p(j,1))
            c23(j)=dcmplx(-p(j,3),p(j,2))
            f(j)=ci
         endif
      enddo
      do i=2,N
         do j=1,i-1
         s(i,j)=two*(p(i,4)*p(j,4)-p(i,1)*p(j,1)-p(i,2)*p(j,2)-p(i,3)*p(j,3))
         za(i,j)=f(i)*f(j)*(c23(i)*dcmplx(rt(j)/rt(i))-c23(j)*dcmplx(rt(i)/rt(j)))

         if (abs(s(i,j)).lt.1d-5) then
         zb(i,j)=-(f(i)*f(j))**2*dconjg(za(i,j))
         else
         zb(i,j)=-dcmplx(s(i,j))/za(i,j)
         endif
         za(j,i)=-za(i,j)
         zb(j,i)=-zb(i,j)
         s(j,i)=s(i,j)
         enddo
      enddo

    end subroutine spinoru

    
    
    
    
    
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


    

END MODULE

