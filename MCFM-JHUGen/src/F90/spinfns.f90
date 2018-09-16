  module spinfns
    use consts_dp
    implicit none

    public :: sc, spb2, spi2, spi5, spb5, psp1, ubar0, v0, pol_dk2mom

  contains

  function ubar0(p,i)
    double complex, intent(in) :: p(4)
    integer, intent(in)       :: i
    ! -------------------------------      
    double complex :: ubar0(4)
    double complex :: fc, fc2
    double precision    :: p0,px,py,pz
    
!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=dreal(p(1))
    px=dreal(p(2))
    py=dreal(p(3))
    pz=dreal(p(4))
!^^^END

    fc2 = p0 + pz 
    fc=sqrt(fc2)

    if (abs(fc2).gt. tol) then 
       
       if (i.eq.1) then 
          ubar0(1)=czero
          ubar0(2)=czero
          ubar0(3)=fc
          ubar0(4)=(px-ci*py)/fc
       elseif (i.eq.-1) then 
          ubar0(1)=(px+ci*py)/fc
          ubar0(2)=-fc
          ubar0(3)=czero
          ubar0(4)=czero
       else
          stop 'ubar0: i out of range' 
       endif
       
    else
!	if (verbose) write(*,*) 'here',i       
       if (i.eq.1) then 
          ubar0(1) = czero
          ubar0(2) = czero
          ubar0(3) = czero
          ubar0(4) = sqrt(cone*two*p0)
       elseif (i.eq.-1) then 
          ubar0(1) = sqrt(cone*(two*p0))
          ubar0(2) = czero
          ubar0(3) = czero
          ubar0(4) = czero
       else
          stop 'ubar0: i out of range' 
       endif
!       if (verbose) write(*,*) 'and here' 
    endif
    
    
  end function ubar0
  
  
  
  ! -- v0  spinor, massless
  function v0(p,i)
    double complex, intent(in) :: p(4)
    integer, intent(in)       :: i
    ! -------------------------------      
    double complex :: v0(4)
    double complex :: fc2, fc
    double precision    :: p0,px,py,pz
    
!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=dreal(p(1))
    px=dreal(p(2))
    py=dreal(p(3))
    pz=dreal(p(4))
!^^^END

    fc2 = p0 + pz 
    fc=sqrt(fc2)

    if (abs(fc2).gt. tol) then 
       
       if (i.eq.1) then 
          v0(1)=czero
          v0(2)=czero
          v0(3)=(px-ci*py)/fc
          v0(4)=-fc
       elseif (i.eq.-1) then
          v0(1)=fc
          v0(2)=(px+ci*py)/fc
          v0(3)=czero
          v0(4)=czero
       else
          stop 'v0: i out of range' 
       endif

    else
       
       if (i.eq.1) then 
          v0(1)=czero
          v0(2)=czero
          v0(3)=sqrt(cone*two*p0)
          v0(4)=czero
       elseif (i.eq.-1) then
          v0(1)=czero
          v0(2)=sqrt(cone*two*p0)
          v0(3)=czero
          v0(4)=czero
       else
          stop 'v0: i out of range' 
       endif
       
    endif
    
  end function v0

  ! -- massless vector polarization subroutine
  function pol_mless(p,i,outgoing)
    double complex, intent(in)    :: p(4)
    integer, intent(in)          :: i
    logical, intent(in),optional :: outgoing
    ! -------------------------------      
    integer :: pol
    double precision :: p0,px,py,pz
    double precision :: pv,ct,st,cphi,sphi
    double complex :: pol_mless(4)
    
!^^^IFmp
!    p0=(p(1)+conjg(p(1)))/two
!    px=(p(2)+conjg(p(2)))/two
!    py=(p(3)+conjg(p(3)))/two
!    pz=(p(4)+conjg(p(4)))/two
!^^^ELSE
    p0=dreal(p(1))
    px=dreal(p(2))
    py=dreal(p(3))
    pz=dreal(p(4))
!^^^END   


    pv=sqrt(abs(p0**2))
    ct=pz/pv
    st=sqrt(abs(1.0d0-ct**2))
    
    if (st < tol) then 
       cphi=1.0d0
       sphi=0.0d0
    else
       cphi= px/pv/st
       sphi= py/pv/st
    endif
    

    ! -- distinguish between positive and negative energies
    if ( p0 > 0.0d0) then  
       pol=i
    else
       pol=-i
    endif
    
    ! -- take complex conjugate for outgoing 
    if (present(outgoing)) then 
       if (outgoing) pol = -pol 
    endif
    
    pol_mless(1)=czero 
    pol_mless(2)=ct*cphi/sqrt2 - ci*pol*sphi/sqrt2
    pol_mless(3)=ct*sphi/sqrt2 + ci*pol*cphi/sqrt2
    pol_mless(4)=-st/sqrt2

  end function pol_mless


  !--------massive vector boson with decay polarization
  ! send by Keith with email on 12Set08
  function pol_dk2mom(plepton,antilepton,i,outgoing)
    integer, intent(in) :: i
    integer :: j
    double complex, intent(in) :: plepton(:),antilepton(:)
    logical, intent(in),optional :: outgoing
    double complex :: pol_dk2mom(4),Ub(4),V(4),q(4),qsq
    
    
    q=plepton+antilepton
    qsq=q(1)**2-q(2)**2-q(3)**2-q(4)**2
    
    Ub(:)=ubar0(plepton,i)
    V(:)=v0(antilepton,-i)
    
    !---Now return in Kirill's notation  1=E,2=px,3=py,4=pz
    !   This is an expression for (-i)/qsq* (-i) Ub(+/-)) Gamma^\mu V(-/+)
    pol_dk2mom(1)=-(Ub(2)*V(4)+V(2)*Ub(4)+Ub(1)*V(3)+V(1)*Ub(3))
    pol_dk2mom(2)=-(-Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)+V(2)*Ub(3))
    pol_dk2mom(3)=-ci*(Ub(1)*V(4)+V(1)*Ub(4)-Ub(2)*V(3)-V(2)*Ub(3))
    pol_dk2mom(4)=-(Ub(2)*V(4)-V(2)*Ub(4)-Ub(1)*V(3)+V(1)*Ub(3))
    
    do j=1,4
       pol_dk2mom(j)=pol_dk2mom(j)/qsq
    enddo

    ! -- do nothing in this case
    if (present(outgoing)) then 
       !if (outgoing) pol_dk2mom = conjg(pol_dk2mom) 
    endif
  end function pol_dk2mom
  

    
    function sc(p1,p2)
      double complex, intent(in) :: p1(:)
      double complex, intent(in) :: p2(:)
      double complex             :: sc
      integer :: sizemin 
      
      sizemin=min(size(p1),size(p2))

      sc = p1(1)*p2(1)
      sc = sc - sum(p1(2:sizemin)*p2(2:sizemin))
      
    end function sc
      

    function spb2(sp,v) 
      double complex, intent(in) :: sp(:),v(:)
      double complex :: spb2(size(sp))
      double complex :: x0(4,4),xx(4,4),xy(4,4)
      double complex :: xz(4,4),x5(4,4)
      double complex :: y1,y2,y3,y4,bp,bm,cp,cm
      integer :: i,i1,i2,i3,Dv,Ds,imax

      Ds = size(sp)

      if (Ds == 4) then
         Dv = 4
      elseif (Ds == 8) then
         Dv = 6
      elseif (Ds == 16) then
         Dv = 8
      else
         stop 'spb2:Dv not allowed'
      endif

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

         xy(1,i)=ci*y4
         xy(2,i)=-ci*y3           
         xy(3,i)=-ci*y2
         xy(4,i)=ci*y1

         xz(1,i)=y3
         xz(2,i)=-y4
         xz(3,i)=-y1
         xz(4,i)=y2

         x5(1,i)=y1
         x5(2,i)=y2
         x5(3,i)=-y3
         x5(4,i)=-y4
      enddo

      if (Dv.eq.4) then 

         do i=1,4
            spb2(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)
         enddo

      elseif (Dv.eq.6) then 
         bp = (v(5)+ci*v(6))
         bm=(v(5)-ci*v(6))

         do i=1,4

            spb2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
                 &          -v(3)*xy(i,1)-v(4)*xz(i,1)+bm*x5(i,2)
            
            i1 = i+4
            spb2(i1)= v(1)*x0(i,2)-v(2)*xx(i,2) &
                 &            -v(3)*xy(i,2)-v(4)*xz(i,2)-bp*x5(i,1)
         enddo
      elseif (Dv.eq.8) then 
         bp=(v(5)+ci*v(6))
         bm=(v(5)-ci*v(6))
         cp=(v(7)+ci*v(8))
         cm=(v(7)-ci*v(8))
         
         do i=1,4

            spb2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
            &          -v(3)*xy(i,1)-v(4)*xz(i,1) &
            &         +bm*x5(i,2)-cm*x5(i,3)

            i1 = i+4
 
            spb2(i1) = v(1)*x0(i,2)-v(2)*xx(i,2) &
            &            -v(3)*xy(i,2)-v(4)*xz(i,2) &
            &            -bp*x5(i,1)+cm*x5(i,4)

            i2 = i1+4

            spb2(i2)=v(1)*x0(i,3)-v(2)*xx(i,3) &
            &             -v(3)*xy(i,3)-v(4)*xz(i,3) & 
            &             +bm*x5(i,4)+cp*x5(i,1)

            i3=i2+4

            spb2(i3)=v(1)*x0(i,4)-v(2)*xx(i,4) &
            &             -v(3)*xy(i,4)-v(4)*xz(i,4) &
            &             -bp*x5(i,3)-cp*x5(i,2)

         enddo
      else
         stop 'spb2: Dv out of bound' 
      endif

    end function spb2



    function spi2(v,sp)
      double complex, intent(in) :: sp(:),v(:)
      double complex :: spi2(size(sp))
      double complex :: x0(4,4),xx(4,4),xy(4,4)
      double complex :: xz(4,4),x5(4,4)
      double complex ::  y1,y2,y3,y4,bp,bm,cp,cm
      integer :: i,i1,i2,i3,imax,Dv,Ds

      Ds = size(sp)

      if (Ds == 4) then
         Dv = 4
      elseif (Ds == 8) then
         Dv = 6
      elseif (Ds == 16) then
         Dv = 8
      else
         stop 'spi2:Dv not allowed'
      endif

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


         xy(1,i)=ci*y4
         xy(2,i)=-ci*y3           
         xy(3,i)=-ci*y2
         xy(4,i)=ci*y1

         xz(1,i)=-y3
         xz(2,i)=y4
         xz(3,i)=y1
         xz(4,i)=-y2

         x5(1,i)=y1
         x5(2,i)=y2
         x5(3,i)=-y3
         x5(4,i)=-y4

      enddo

      if(Dv.eq.4) then 

         do i=1,4

            spi2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
            &           -v(3)*xy(i,1)-v(4)*xz(i,1)
         enddo

      elseif (Dv.eq.6) then 
         bp = (v(5)+ci*v(6))
         bm=(v(5)-ci*v(6))


         do i=1,4

            spi2(i)=v(1)*x0(i,1)-v(2)*xx(i,1) &
            &           -v(3)*xy(i,1)-v(4)*xz(i,1) &
            &           -bp*x5(i,2)

            i1=i+4

            spi2(i1)=v(1)*x0(i,2)-v(2)*xx(i,2) &
            &             -v(3)*xy(i,2)-v(4)*xz(i,2) &
            &             +bm*x5(i,1)

         enddo

      elseif (Dv.eq.8) then 
         
         bp = (v(5)+ci*v(6))
         bm=(v(5)-ci*v(6))
         cp=(v(7)+ci*v(8))
         cm=(v(7)-ci*v(8))

         do i=1,4

            spi2(i)=v(1)*x0(i,1)-v(2)*xx(i,1)&
            &           -v(3)*xy(i,1)-v(4)*xz(i,1)&
            &           -bp*x5(i,2)+ cp*x5(i,3)

            i1=i+4

            spi2(i1)=v(1)*x0(i,2)-v(2)*xx(i,2)&
            &             -v(3)*xy(i,2)-v(4)*xz(i,2)&
            &             +bm*x5(i,1)-cp*x5(i,4)

            i2=i1+4

            spi2(i2)=v(1)*x0(i,3)-v(2)*xx(i,3)&
            &           -v(3)*xy(i,3)-v(4)*xz(i,3)&
            &          -bp*x5(i,4)-cm*x5(i,1)

            i3=i2+4

            spi2(i3)=v(1)*x0(i,4)-v(2)*xx(i,4)&
            &           -v(3)*xy(i,4)-v(4)*xz(i,4)&
            &           +bm*x5(i,3)+cm*x5(i,2)


         enddo

      else
         stop 'spi2: Dv out of bounds' 
      end if

    end function spi2


    function  psp1(sp1,sp2)
      double complex, intent(in) :: sp1(:)
      double complex, intent(in) :: sp2(:)
      double complex :: psp1
      
      psp1 = sum(sp1(1:)*sp2(1:))
      
    end function psp1


    ! -- multiplication of spinor with gamma_5 on the left  
    function spi5(sp)
      double complex, intent(in) :: sp(:)
      double complex :: spi5(size(sp))
      integer :: i,j,imax, Ds  

      Ds = size(sp)
      imax = Ds/4
      do i=1,imax
         spi5(4*(i-1)+1) = sp(4*(i-1)+1) 
         spi5(4*(i-1)+2) = sp(4*(i-1)+2) 
         spi5(4*(i-1)+3) = -sp(4*(i-1)+3) 
         spi5(4*(i-1)+4) = -sp(4*(i-1)+4) 
      enddo
!      do i=1,imax
!         do j=1,4
!            if (j <2) then 
!               spi5(4*(i-1)+j) = sp(4*(i-1)+j) 
!            else
!               spi5(4*(i-1)+j) = -sp(4*(i-1)+j) 
!            endif
!         enddo
!      enddo

    end function spi5


    ! -- multiplication of bspinor with gamma_5 on the right  
    function spb5(sp)
      double complex, intent(in) :: sp(:)
      double complex :: spb5(size(sp))
      integer :: i,j,imax, Ds  
      
      spb5 = spi5(sp) 
      !Ds = size(sp)
      !imax = Ds/4
      !do i=1,imax
      !   do j=1,4
      !      if (j <2) then 
      !         spb5(4*(i-1)+j) = sp(4*(i-1)+j) 
      !      else
      !         spb5(4*(i-1)+j) = -sp(4*(i-1)+j) 
      !      endif
      !   enddo
      !enddo

    end function spb5



  end module spinfns