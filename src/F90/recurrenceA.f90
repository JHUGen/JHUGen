  !------------------------------------------------------------------------!
  ! Authors: Giulia Zanderighi...	                                   !
  !------------------------------------------------------------------------!

!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECrecurrenceGbitsone.f90.

module recurrenceA
  use consts_dp
  use spinfns
  implicit none

  public :: vggg, vgggg,vqg,vgq,vbqg,vgbq,vbqq,vbqW,vWq!,g 
  public :: vgluon,f,bf

  logical :: verbose = .false. 
  
  private
  
  contains 
  
  function vggg(e1,k1,e2,k2)
    double complex, intent(in) :: e1(:), e2(:)
    double complex, intent(in) :: k1(:), k2(:)
    double complex             :: vggg(size(e1))
    ! -----------------------------------------
    double complex:: sk1e2,se1e2,sk2e1,xx

    sk1e2=sc(k1,e2)
    sk2e1=sc(k2,e1)
    se1e2=sc(e1,e2)
    
    xx=ci*sqrt2
    
    vggg = xx*(-sk1e2*e1+sk2e1*e2+se1e2/two*(k1-k2))

  end function vggg


  function  vgggg(e1,e2,e3)              
    double complex, intent(in) :: e1(:),e2(:),e3(:)
    double complex             :: vgggg(size(e1))
    ! -----------------------------------------
    double complex:: se1e3,se2e3,se1e2

    se1e3=sc(e1,e3)
    se2e3=sc(e2,e3)
    se1e2=sc(e1,e2)
    
    vgggg = ci*(e2*se1e3-chalf*(e1*se2e3+ e3*se1e2))

  end function vgggg

  function vqg(sp,e1)
    double complex, intent(in) :: e1(:)
    double complex, intent(in) :: sp(:)
    double complex :: vqg(size(sp))

    vqg = ci/sqrt2*spb2(sp,e1)

  end function vqg


  function vgq(e1,sp)
    double complex, intent(in) :: e1(:)
    double complex, intent(in) :: sp(:)
    double complex :: vgq(size(sp))

    vgq = -ci/sqrt2*spb2(sp,e1)

  end function vgq


  function vbqg(sp,e1)
    double complex, intent(in) :: e1(:)
    double complex, intent(in) :: sp(:)
    double complex :: vbqg(size(sp))

    vbqg = -ci/sqrt2*spi2(e1,sp)

  end function vbqg


  function vgbq(e1,sp)
    double complex, intent(in) :: e1(:)
    double complex, intent(in) :: sp(:)
    double complex :: vgbq(size(sp))

    vgbq = ci/sqrt2*spi2(e1,sp)

  end function vgbq


  function vbqq(Dv,sp1,sp2)
    double complex, intent(in) :: sp1(:), sp2(:)
    integer, intent(in) ::  Dv
    ! -------------------------------------------       
    double complex :: vbqq(Dv)
    double complex :: rr, va(Dv),sp1a(size(sp1))
    integer :: i


    va=czero
    vbqq=czero

    do i=1,Dv

       if (i.eq.1) then 
          va(1)=cone
       else
          va(i)=-cone
       endif

       sp1a=spb2(sp1,va)

       rr=psp1(sp1a,sp2)
       !       rr=-ci/sqrt2*psp1(sp1a,sp2)

       if (i.eq.1) then 
          vbqq = vbqq + rr*va
       else
          vbqq = vbqq - rr*va
       endif

       va(i)=czero 
    enddo
    vbqq = -ci/sqrt2*vbqq

  end function vbqq

  function vbqW(sp,e1)
    double complex, intent(in) :: e1(:)
    double complex, intent(in) :: sp(:)
    double complex :: vbqW(size(sp)),sp5(size(sp))

!    if (v_coupling == v_vector) then 
       vbqW = -ci*spb2(sp,e1)
!    else
!       sp5 = -ci*spb2(sp,e1)
!       vbqW = spb5(sp5) 
!    endif

  end function vbqW


  function vWq(e1,sp)
    double complex, intent(in) :: e1(:)
    double complex, intent(in) :: sp(:)
    double complex :: vWq(size(sp)),sp5(size(sp))

!    if (v_coupling == v_vector) then 
       vWq = -ci*spi2(e1,sp)
!    else
!       sp5 = spi5(sp)
!       vWq = -ci*spi2(e1,sp5)
!    endif

  end function vWq

  function g(e,k,giarray,pol_int) result(res)
    double complex, intent(in) :: e(:,:),k(:,:)
    integer, intent(in),optional       :: giarray(:),pol_int 
    double complex             :: res(size(e,dim=1))
    res = vgluon(e,k,giarray,pol_int)
  end function g

  !  ---  recurrence for gluon currents
  recursive function vgluon(e,k,giarray,pol_int) result(res)
    double complex, intent(in) :: e(:,:),k(:,:)
    integer, intent(in),optional       :: giarray(:),pol_int 
    double complex             :: res(size(e,dim=1))
    ! -----------------------------------------------------
    double complex :: k1(size(e,dim=1)),k2(size(e,dim=1)),k3(size(e,dim=1))
    double complex :: e1(size(e,dim=1)),e2(size(e,dim=1)),e3(size(e,dim=1))
    double complex :: tmp(size(e,dim=1))
    double complex :: k1sq, k2sq, k3sq
    integer       :: npart, m, m1
    logical :: done 

    npart = size(e,dim=2)

    if (npart == 0) then 
       res = czero 
       return 
    endif
    if (npart == 1) then 
       res = e(:,1)
       return 
    endif

    !if (verbose) write(*,*) 'entering vgluon:npart',npart
    !if (verbose .and. present(giarray)) write(*,*) 'entering vgluon:gi',giarray
    done = .false. 
!    if (present(giarray)) then 
!       if (npart /= size(giarray)) stop 'vgluon: npart/= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray)
!       if (done) return 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'vgluon: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif


    if (npart == 1) then 
       res = e(:,1)

    elseif (npart == 2) then 
       res = vggg(e(:,1),k(:,1),e(:,2),k(:,2))

    else

       res = czero

       do m=1,npart-1

          k1=sum(k(:,1:m),dim=2)
          k2=sum(k(:,m+1:npart),dim=2)

          e1 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)
          e2 = vgluon(e(:,m+1:npart),k(:,m+1:npart),&
               &giarray(m+1:npart),pol_int)

          tmp = vggg(e1,k1,e2,k2)

          if (m > 1) then 
             k1sq = sc(k1,k1)
             if (abs(k1sq) > propcut) then 
                tmp = -ci*tmp/k1sq
             else
                tmp = czero
             endif
          endif

          if (m + 1 < npart) then 
             k2sq = sc(k2,k2)
             if (abs(k2sq) > propcut) then 
                tmp = -ci*tmp/k2sq
             else
                tmp = czero 
             endif
          endif

          res = res + tmp

          if (m <= npart-2) then 

             do m1=m+1,npart-1
                e2=vgluon(e(:,m+1:m1),k(:,m+1:m1),giarray(m+1:m1),pol_int)
                e3=vgluon(e(:,m1+1:npart),k(:,m1+1:npart),&
                     &giarray(m1+1:npart),pol_int)
                k2 = sum(k(:,m+1:m1),dim=2)
                k3 = sum(k(:,m1+1:npart),dim=2)
                tmp = vgggg(e1,e2,e3)
                if (m > 1) then  
                   k1sq = sc(k1,k1)
                   if(abs(k1sq) > propcut) then 
                      tmp = -ci*tmp/k1sq
                   else 
                      tmp = czero
                   endif
                endif
                if (m+1 < m1) then 
                   k2sq = sc(k2,k2)
                   if(abs(k2sq) > propcut) then 
                      tmp = -ci*tmp/k2sq
                   else
                      tmp = czero
                   endif
                endif
                if (m1+1 < npart) then 
                   k3sq = sc(k3,k3)
                   if(abs(k3sq) > propcut) then 
                      tmp = -ci*tmp/k3sq
                   else
                      tmp = czero
                   endif
                endif
                res = res + tmp
             enddo

          endif

       enddo

    endif

    ! -- store current 
!    if (present(giarray))   call store_result(pol_int,res,giarray)


  end function vgluon



  ! ---- fermion current 
  recursive function f(e,k,sp,p,flout,flin,ms,giarray,qiarray,pol_int) &
       &result(res)
    double complex, intent(in) :: e(:,:), k(:,:), sp(:), p(:)
    character(len=3),  intent(in)  :: flin    ! flavor of off-shell f-line
    character(len=3),  intent(in)  :: flout   ! flavor of on-shell f-line
    integer, intent(in) ::  ms
    integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int 
    ! -----------------------------------------------------------------------
    double complex :: res(size(sp))
    double complex :: tmp(size(sp))
    double complex :: k1(size(p))
    double complex :: k2(size(p))
    double complex :: sp2(size(sp))
    double complex :: e1(size(e,dim=1))
    double complex :: k1sq,k2sq
!    real(dp)    :: mass,mass2
    integer       :: ms1,m,ng1, ng2, ngluon
    logical       :: done 

    ! NB: this is must not cached, because caching does not have info about 
    !     flavour 
    if (flout.ne.flin) then 
       res = czero 
       return
    endif

    if (size(e,dim=2) == 0) then 
       res = sp
       return
    endif

    !if (verbose) write(*,*) 'entering f:ngluon,ng1,ng2',&
    !     &size(e,dim=2),ms,size(e,dim=2)-ms
    !if (verbose .and. present(giarray)) write(*,*) 'entering f:giarray',giarray
    !if (verbose) write(*,*) 'entering f:qiarray',qiarray
    
    done = .false. 
!    if (present(giarray)) then 
!       !if (size(qiarray) /= 1) stop 'f: wrong size qiarray'
!       !if (size(e,dim=2) /= size(giarray)) stop 'f: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray)
!       if (done) return ! XXX 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'f: giarray missing'
!          i_warn = i_warn+1
!       endif
!!       stop 
!    endif

    !mass = mt
    !mass2 = mass**2

    ngluon = size(e,dim=2)
    ng1 = ms           ! #gluons to the left of the f-line 
    ng2 = ngluon - ms  ! #gluons to the right of the f-line


    if ((ng1 < 0) .or. (ng2 < 0)) write(*,*) 'WRONG DEFINITION OF CURRENT f'

    if (ngluon == 0) then 
       res = sp

    else

       res = czero

       do m=0,ng2-1
          if (ng1+1+m<=ngluon) then 
             k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
          else
             k1 = czero 
          endif

          e1=vgluon(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon), &
               &giarray(ng1+1+m:ngluon),pol_int)

          k1sq=sc(k1,k1)

          if (1<=ng1+m) then 
             k2 = sum(k(:,1:ng1+m),dim=2)
          else
             k2 = czero
          endif

          k2 = k2 + p
          k2sq = sc(k2,k2)!-mass2
          sp2 = f(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,flout,flin,ng1,&
               &giarray(1:ng1+m),qiarray,pol_int)

          if (ng1 >0.or.m>0) sp2 = spb2(sp2,k2)!+mass*sp2

          tmp = vqg(sp2,e1)

          if (m < ng2-1)  then      
             if(abs(k1sq) > propcut) then 
                tmp = -ci/k1sq*tmp
             else
                tmp = czero
             endif
          endif

          if (ng1>0.or.m>0) then 
             if (abs(k2sq) > propcut) then 
                tmp =  ci/k2sq*tmp
             else
                tmp = czero 
             endif
          endif

          res = res + tmp


       enddo


       do m=1,ng1

          k1 = sum(k(:,1:m),dim=2)
          e1=vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)

          k1sq = sc(k1,k1)

          if (m+1<=ngluon) then 
             k2 = sum(k(:,m+1:ngluon),dim=2)
          else
             k2 = czero 
          endif

          k2 = k2 + p
          k2sq = sc(k2,k2) !- mass2
          ms1 = ng1 - m
          sp2=f(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,flout,flin,ms1,&
               &giarray(m+1:ngluon),qiarray,pol_int)

          if (ng2 > 0.or.m < ng1) sp2 = spb2(sp2,k2)!+mass*sp2

          tmp = vgq(e1,sp2)

          if (m > 1) then  
             if (abs(k1sq) > propcut) then 
                tmp=-ci/k1sq*tmp
             else
                tmp = czero 
             endif
          endif

          if (ng2 > 0.or. m < ng1) then 
             if (abs(k2sq) > propcut) then 
                tmp=ci/k2sq*tmp
             else
                tmp = czero 
             endif
          endif

          res = res + tmp

       enddo

    endif

    ! -- store current 
!    if (present(giarray)) then 
!       call store_result(pol_int,res,giarray,qiarray)
!       !if (verbose) write(*,*) 'f: storing',ngluon,ng1,ng2,res
!    endif

  end function f



  ! ---- bar-fermion current
  recursive function bf(e,k,sp,p,flout,flin,ms,giarray,qiarray,pol_int) &
       &result(res)
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp(:), p(:)
    integer, intent(in) ::  ms
    character, intent(in) :: flout*3  ! flavor of on-shell f-line
    character, intent(in) :: flin*3   ! flavor of off-shell f-line
    integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int 
    ! -----------------------------------------------------------------------
    double complex             :: res(size(sp))
    double complex             :: tmp(size(sp))
    double complex             :: k1(size(p))
    double complex             :: k2(size(p))
    double complex             :: sp2(size(sp))
    double complex             :: e1(size(e,dim=1))
    double complex             :: k1sq,k2sq
!    real(dp)                :: mass,mass2
    integer                   :: ms1,m,ng1, ng2, ngluon
    logical                   :: done 

    !if (verbose) write(*,*) 'entering bf:ngluon,ng1,ng2',&
    !     &size(e,dim=2),ms,size(e,dim=2)-ms
    !if (verbose .and. present(giarray)) write(*,*) 'entering bf:giarray',giarray
    !if (verbose) write(*,*) 'entering bf:qiarray',qiarray

    if (flout.ne.flin) then 
       res = czero 
       return 
    endif

    if (size(e,dim=2) == 0) then 
       res = sp
       return
    endif


    done = .false. 
!    if (present(giarray)) then 
!       !if (size(qiarray) /= 1) stop 'bf: wrong size qiarray'
!       !if (size(e,dim=2) /= size(giarray)) then 
!       !   write(*,*) size(giarray), size(e,dim=2)
!       !   stop 'bf: ng= size(giarray)' 
!       !endif
!       call memory_check(pol_int,res,done,giarray,qiarray)
!       if (done) return 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'bf: giarray missing'
!          i_warn = i_warn+1
!       endif
!!       stop 
!    endif

!    mass = mt
!    mass2 = mt**2

    ngluon = size(e,dim=2)
    ng1 = ms   !#gluons to the left of a f-line 
    ng2 = ngluon - ms  !#gluons to the right of the f-line

    !if (verbose) write(*,*) 'in function bf', ng1, ngluon

    if ((ng1 < 0) .and. (ng2 < 0)) write(*,*) 'bf: WRONG DEFINITION OF CURRENT'

    if (ngluon == 0) then 
       res = sp

    else

       res = czero

       do m=0,ng2-1
          if (ng1+1+m<=ngluon) then 
             k1 = sum(k(:,ng1+1+m:ngluon),dim=2)
          else
             k1 = czero 
          endif
          e1=vgluon(e(:,ng1+1+m:ngluon),k(:,ng1+1+m:ngluon),&
               &giarray(ng1+1+m:ngluon),pol_int)

          k1sq=sc(k1,k1)

          if (1<=ng1+m) then 
             k2 = sum(k(:,1:ng1+m),dim=2)
          else
             k2 = czero 
          endif
          k2 = -k2 - p
          k2sq = sc(k2,k2)!-mass2
          sp2 = bf(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,flout,flin,ng1,&
               &giarray(1:ng1+m),qiarray,pol_int)

          if (ng1 >0.or.m>0) sp2 = spi2(k2,sp2)!+mass*sp2

          tmp = vbqg(sp2,e1)

          if (m < ng2-1) then       
             if (abs(k1sq) > propcut) then 
                tmp = -ci/k1sq*tmp
             else
                tmp = czero 
             endif
          endif
          if (ng1>0.or.m>0) then 
             if (abs(k2sq) > propcut) then 
                tmp =  ci/k2sq*tmp
             else
                tmp = czero
             endif
          endif

          res = res + tmp


       enddo


       do m=1,ng1

          k1 = sum(k(:,1:m),dim=2)
          e1=vgluon(e(:,1:m),k(:,1:m), &
               &giarray(1:m),pol_int)
          k1sq = sc(k1,k1)

          if (m+1<=ngluon) then 
             k2 = sum(k(:,m+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2 = -k2 - p
          k2sq = sc(k2,k2) !- mass2
          ms1 = ng1 - m
          sp2=bf(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,flout,flin,ms1,&
               &giarray(m+1:ngluon),qiarray,pol_int)

          if (ng2 > 0.or.m < ng1) sp2 = spi2(k2,sp2)!+mass*sp2

          tmp = vgbq(e1,sp2)

          if (m > 1) then 
             if (abs(k1sq) > propcut) then 
                tmp=-ci/k1sq*tmp
             else 
                tmp = czero 
             endif
          endif

          if (ng2 > 0.or. m < ng1) then 
             if (abs(k2sq) > propcut) then 
                tmp=ci/k2sq*tmp
             else
                tmp = czero 
             endif
          endif

          res = res + tmp

       enddo

    endif


!    ! -- store current 
!    if (present(giarray)) then
!       call store_result(pol_int,res,giarray,qiarray)
!       !if (verbose) write(*,*) 'bf: storing',ngluon,ng1,ng2,res
!    endif
  end function bf

end module recurrenceA
