
  !--- spinor products
  subroutine spinoru(j,p,zza,zzb,ss1) 
    implicit none 
    integer, intent(in) :: j
    complex(dp), intent(out)  :: zza(j,j), zzb(j,j)
    real(dp), intent(out) :: ss1(j,j)
    real(dp), intent(in) :: p(4,j)
    integer :: i1,i2
    include 'includeVars.F90'

    zza = czero
    zzb = czero
    ss1 = zero

    do i1=1,j
       do i2=i1+1,j
          
          zza(i1,i2) = aa(p(:,i1),p(:,i2))
          zzb(i1,i2) = bb(p(:,i1),p(:,i2)) 
          zza(i2,i1) = -zza(i1,i2)
          zzb(i2,i1) = -zzb(i1,i2)

          ss1(i1,i2) = zza(i1,i2)*zzb(i2,i1)
          ss1(i2,i1) = ss1(i1,i2)

       enddo
    enddo

    return 
  end subroutine spinoru


  function aa(p2,p1)
    implicit none
    real(dp), intent(in) :: p1(4), p2(4) 
    complex(dp) :: aa

    aa = sum(bspa(p2)*spa(p1))

  end function aa


  ! [1 2] 

  function bb(p2,p1)
    implicit none
    real(dp), intent(in) :: p1(4), p2(4) 
    complex(dp) :: bb
    
    bb = sum(bspb(p2)*spb(p1))
    
  end function bb


  ! | p > spinor

  function spa(p)  
    implicit none 
    real(dp), intent(in) :: p(4)
    complex(dp) :: spa(4)
      
    spa = u0(p,1)

  end function spa


  ! | p ] spinor

  function spb(p) 
    implicit none  
    real(dp), intent(in) :: p(4)
    complex(dp) :: spb(4)
      
    spb = u0(p,-1)

  end function spb


  ! < p | spinor

  function bspa(p)  
    implicit none 
    real(dp), intent(in) :: p(4)
    complex(dp) :: bspa(4)
  
    bspa = ubar0(p,-1)

  end function bspa


  ! [ p | spinor

  function bspb(p)  
    implicit none 
    real(dp), intent(in) :: p(4)
    complex(dp) :: bspb(4)
    
    bspb = ubar0(p,1)

  end function bspb

  ! -- ubar spinor, massless
  function ubar0(p,i) 
    implicit none 
    real(dp), intent(in) :: p(4)
    integer, intent(in) :: i
    complex(dp) :: ubar0(4)
    real(dp)    :: p0,px,py,pz, theta, phi, rrr
    complex(dp) :: cp0
    integer :: sgn
    logical :: flag_nan
    include 'includeVars.F90'

    flag_nan = .false.

    p0=p(1)
    px=p(2)
    py=p(3)
    pz=p(4)
    

    cp0 = cmplx(p0,kind=dp)

    if (p0.eq.zero) then
       write(6,*) 'error in ubar -> p0=0'
    elseif(px.eq.zero.and.py.eq.zero) then
       if ((pz/p0).gt.0.0_dp) theta = zero
       if ((pz/p0).lt.0.0_dp) theta = pi
       phi = zero
    elseif(px.eq.zero.and.py.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = py/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = asin(rrr)
    elseif(py.eq.zero.and.px.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    else
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    endif
    

    call get_NaN(theta,flag_nan)
    if (flag_nan.eqv. .true.) then
       write(6,*) 'ubar-th', p0,px,py,pz
       stop
    endif


    call get_NaN(phi,flag_nan)
    if (flag_nan.eqv. .true.) then
       write(6,*) 'ubar-phi', p0,px,py,pz
       stop
    endif


    if (i.eq.1) then
       ubar0(1)=czero
       ubar0(2)=czero
       ubar0(3)=sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
       ubar0(4)=sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),-sin(phi),kind=dp)
    elseif(i.eq.-1) then
       ubar0(1)= sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),sin(phi),kind=dp)
       ubar0(2)=-sqrt2*sqrt(cp0)*cmplx(abs(cos(theta/two)),0.0_dp,kind=dp)
       ubar0(3)=czero
       ubar0(4)=czero
    else
       stop 'ubar0: i out of range'
    endif

  end function ubar0


  ! -- u0  spinor, massless
  function u0(p,i) 
    implicit none 
    real(dp), intent(in) :: p(4)
    integer, intent(in) :: i
    complex(dp) :: u0(4)
    real(dp)    :: p0,px,py,pz, theta, phi, rrr
    complex(dp) :: cp0
    integer :: sgn
    logical :: flag_nan
    include 'includeVars.F90'
    
    flag_nan = .false.
    
    p0=p(1)
    px=p(2)
    py=p(3)
    pz=p(4)
    
    cp0 = cmplx(p0,kind=dp)

    if (p0.eq.zero) then
       write(6,*) 'error in v0 -> p0=0'
    elseif(px.eq.zero.and.py.eq.zero) then
       if ((pz/p0).gt.0.0_dp) theta = zero
       if ((pz/p0).lt.0.0_dp) theta = pi
       phi = zero
    elseif(px.eq.zero.and.py.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       if (py/p0.gt.zero)  phi = pi/two
       if (py/p0.lt.zero)  phi = 3.0_dp*pi/two
    elseif(py.eq.zero.and.px.ne.zero) then
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       if (px/p0.gt.zero)  phi = zero
       if (px/p0.lt.zero)  phi = pi
    else
       rrr = pz/p0
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       theta=acos(rrr)
       rrr = px/p0/sin(theta)
       if (rrr.lt.-1.0_dp) rrr = -1.0_dp
       if (rrr.gt.1.0_dp)  rrr = 1.0_dp
       phi = acos(rrr)
       if (py/p0.gt.zero) phi = phi
       if (py/p0.lt.zero) phi = -phi
    endif

    call get_NaN(theta,flag_nan)
    if (flag_nan.eqv. .true.) then
       write(6,*) 'th-v0', p0,px,py,pz
       write(6,*) 'ratio', pz/p0, acos(pz/p0)
       stop
    endif
    

    call get_NaN(phi,flag_nan)
    if (flag_nan.eqv. .true.) then
       write(6,*) 'ph-v0', p0,px,py,pz, px/p0/sin(theta)
       stop
    endif
    
    if (i.eq.-1) then
       u0(1)=czero
       u0(2)=czero
       u0(3)=sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),-sin(phi),kind=dp)
       u0(4)=-sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
    elseif (i.eq.1) then
       u0(1)= sqrt2*sqrt(cp0)*cmplx(cos(theta/two),0.0_dp,kind=dp)
       u0(2)= sqrt2*sqrt(cp0)*sin(theta/two)*cmplx(cos(phi),sin(phi),kind=dp)
       u0(3)=czero
       u0(4)=czero
    else
       stop 'u0: i out of range'
    endif

  end function u0

  subroutine get_NaN(value,flag_nan)
    implicit none
    real(dp), intent(in)  :: value
    logical, intent(out) :: flag_nan
    
    flag_nan = .false.
    
    if (.not.value.le.0.0_dp .and. .not.value.gt.0.0_dp ) then
       flag_nan =.true.
    endif
    
  end subroutine get_NaN

  function scr(p1,p2)   !scalar product of real vectors 
    implicit none 
    real(dp), intent(in) :: p1(4), p2(4)
    real(dp) :: scr
    scr = p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)-p1(4)*p2(4)
  end function scr

 !--- scalar products, 2*pi.pj for massless particles
  subroutine sprodu(j,p,sprod)
    implicit none
    integer, intent(in) :: j
    real(dp), intent(in) :: p(4,j)
    real(dp), intent(out) :: sprod(j,j)
    integer :: i1, i2
    include 'includeVars.F90'

    sprod = zero

    do i1 = 1, j
       do i2 = i1+1, j
          sprod(i1,i2) = two * scr(p(:,i1),p(:,i2))
          sprod(i2,i1) = sprod(i1,i2)
       enddo
    enddo

    return

  end subroutine sprodu
