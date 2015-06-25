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
          pout(1)=p(1)
          pout(2)=p(2)
          pout(3)=p(3)
          pout(4)=p(4)          
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

