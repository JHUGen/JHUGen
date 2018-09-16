  !------------------------------------------------------------------------!
  ! Authors: Giulia Zanderighi...	                                   !
  !------------------------------------------------------------------------!

!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECrecurrenceGbitsthree.f90.

module recurrenceC
  use consts_dp
  use spinfns
  use recurrenceA
  use recurrenceB
  implicit none

  public :: fW,bfW,gW_fbf,gW_bff,fW_bffbf
!$$$   public :: fW_bffbf_2W
!$$$   public :: fsW_fbfbf
!$$$   public :: gW_sbsfbf,bfW_fbff,fW_bffbf_1, gW_sbsfbf_1, &
!$$$        & gW_sbsfbf_2,gW_sbsfbf_3
!$$$   public ::  fW_bffbffbf, fW_bffbf_2
!$$$   public :: fW_fbfbf ! for fermion loops with Z 
  
  logical :: verbose = .false. 


  private

  contains 

  !----- currents with the W boson
  recursive function fW(e,k,sp,p,flon,floff,&
       &eW,kW,ms,giarray,qiarray,Wid,pol_int) result(res)
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp(:), p(:)
    double complex, intent(in) :: eW(:), kW(:)
    integer, intent(in) ::  ms
    character, intent(in) :: flon*3, floff*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
    ! -----------------------------------------------------------------------
    integer             :: ms1,m,ng1, ng2
    integer :: ngluon
    double complex             :: res(size(sp))
    double complex             :: res_stored(size(sp))
    double complex             :: tmp(size(sp))
    double complex             :: k1(size(p))
    double complex             :: k2(size(p))
    double complex             :: sp2(size(sp))
    double complex             :: e1(size(e,dim=1))
    double complex             :: k1sq,k2sq!,k3sq
!    real(dp) :: mass,mass2
    logical                   :: done 

    !if (verbose) write(*,*) 'entering fW:ng1,ng2',ms,size(e,dim=2)-ms
    !if (verbose .and. present(giarray)) write(*,*) 'entering fW:gi',giarray
    !if (verbose .and. present(giarray)) write(*,*) 'entering fW:qi',qiarray
    !if (verbose .and. present(Wid)) write(*,*) 'entering fW:Wid',Wid


    if ((flon.eq.'str'.or.floff.eq.'str') .or. (flon.eq.'dwn'.or.floff.eq.'dwn')) then 
       res = czero 
       return 
    endif

    done = .false. 
!    if (present(giarray)) then 
!       !if (size(qiarray) /= 1) stop 'fW: wrong size qiarray'
!       !if (size(e,dim=2) /= size(giarray)) stop 'fW: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!       if (done) return 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'fw: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif

!    mass = mt
!    mass2 = mass**2

    ngluon = size(e,dim=2)
    ng1 = ms   !#gluons to the left of a f-line 
    ng2 = ngluon - ms  !#gluons to the right of the f-line

    if ((ng1 < 0) .or. (ng2 < 0)) write(*,*) 'WRONG DEFINITION OF CURRENT fW'



    if (ngluon == 0) then 
       res = vbqW(sp,eW)
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
          k2 = k2 + p + kW
          k2sq = sc(k2,k2)!-mass2
          sp2 = fW(e(:,1:ng1+m),k(:,1:ng1+m),&
               &sp,p,flon,floff,eW,kW,ng1,&
               &giarray(1:ng1+m),qiarray,Wid,pol_int) 

          sp2 = spb2(sp2,k2)!+mass*sp2

          tmp = vqg(sp2,e1)

          if (m < ng2-1)  then      
             if(abs(k1sq) > propcut) then 
                tmp = -ci/k1sq*tmp
             else
                tmp = czero
             endif
          endif


          if (abs(k2sq) > propcut) then 
             tmp =  ci/k2sq*tmp
          else
             tmp = czero 
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
          k2 = k2 + p + kW
          k2sq = sc(k2,k2) !- mass2
          ms1 = ng1 - m
          sp2=fW(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,&
               &flon,floff,eW,kW,ms1,&
               &giarray(m+1:ngluon),qiarray,Wid,pol_int) 

          sp2 = spb2(sp2,k2)!+mass*sp2

          tmp = vgq(e1,sp2)

          if (m > 1) then  
             if (abs(k1sq) > propcut) then 
                tmp=-ci/k1sq*tmp
             else
                tmp = czero 
             endif
          endif


          if (abs(k2sq) > propcut) then 
             tmp=ci/k2sq*tmp
          else
             tmp = czero 
          endif


          res = res + tmp
       enddo

       ! -------- the new term, finally 

       sp2 = f(e,k,sp,p,flon,flon,ms,giarray,qiarray,pol_int)    
       if (1<=ngluon) then 
          k2 = sum(k(:,1:ngluon),dim=2)
       else
          k2 = czero 
       endif
       k2 = k2 + p 
       k2sq = sc(k2,k2) !- mass**2

       sp2 = spb2(sp2,k2)!+mass*sp2

       tmp = vbqW(sp2,eW)

       if (abs(k2sq) > propcut) then 
          tmp = ci/k2sq*tmp
       else 
          tmp = czero
       endif

       res = res + tmp
          
    endif

    ! -- store current 
!    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)

  end function fW


  recursive function bfW(e,k,sp,p,flon,floff,eW,kW,&
       &ms,giarray,qiarray,Wid,pol_int) result(res)
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp(:), p(:)
    double complex, intent(in) :: ew(:), kW(:)
    integer, intent(in) ::  ms
    character, intent(in) :: flon*3, floff*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
    ! -----------------------------------------------------------------------
    integer             :: ms1,m,ng1, ng2
    integer :: ngluon
    double complex             :: res(size(sp))
    double complex             :: tmp(size(sp))
    double complex             :: k1(size(p))
    double complex             :: k2(size(p))
    double complex             :: sp2(size(sp))
    double complex             :: e1(size(e,dim=1))
    double complex             :: k1sq,k2sq
    !real(dp)                :: mass,mass2
    logical                   :: done 

    !if (verbose) write(*,*) 'entering bfW'

    if ((flon.eq.'str'.or.floff.eq.'str') .or. (flon.eq.'dwn'.or.floff.eq.'dwn')) then 
       res = czero 
       return 
    endif

    done = .false. 
!    if (present(giarray)) then 
!       !if (size(qiarray) /= 1) stop 'bfW: wrong size qiarray'
!       !if (size(e,dim=2) /= size(giarray)) stop 'bfW: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!       if (done) return 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'bfW: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif

    !mass = mt
    !mass2 = mass**2

    ngluon = size(e,dim=2)
    ng1 = ms   !#gluons to the left of a f-line 
    ng2 = ngluon - ms  !#gluons to the right of the f-line


    if (ng2 < 0) write(*,*) 'WRONG DEFINITION OF CURRENT D'


    if (ngluon == 0) then 
       res = vWq(eW,sp)

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

          sp2 = bfW(e(:,1:ng1+m),k(:,1:ng1+m),sp,p,&
               &flon,floff,eW,kW,ng1,&
               &giarray(1:ng1+m),qiarray,Wid,pol_int)    
          if (1<=ng1+m) then 
             k2 = sum(k(:,1:ng1+m),dim=2)
          else
             k2 = czero 
          endif
          k2 = -k2 - p - kW
          k2sq = sc(k2,k2)!-mass**2


          sp2 = spi2(k2,sp2)!+mass*sp2

          tmp = vbqg(sp2,e1)

          if (m < ng2-1) then       
             if (abs(k1sq) > propcut) then 
                tmp = -ci/k1sq*tmp
             else
                tmp = czero 
             endif
          endif


          if (abs(k2sq) > propcut) then 
             tmp =  ci/k2sq*tmp
          else
             tmp = czero
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
          k2 = -k2 - p - kW
          k2sq = sc(k2,k2) !- mass**2
          ms1 = ng1 - m
          sp2=bfW(e(:,m+1:ngluon),k(:,m+1:ngluon),sp,p,&
               &flon,floff,eW,kW,ms1,&
               &giarray(m+1:ngluon),qiarray,Wid,pol_int)    

          sp2 = spi2(k2,sp2)!+mass*sp2

          tmp = vgbq(e1,sp2)

          if (m > 1) then 
             if (abs(k1sq) > propcut) then 
                tmp=-ci/k1sq*tmp
             else 
                tmp = czero 
             endif
          endif

          if (abs(k2sq) > propcut) then 
             tmp=ci/k2sq*tmp
          else
             tmp = czero 
          endif

          res = res + tmp

       enddo


       !----------finally the new term

       sp2 = bf(e,k,sp,p,flon,flon,ms,giarray,qiarray,pol_int)    
       if (1<=ngluon) then 
          k2 = sum(k(:,1:ngluon),dim=2)
       else
          k2 = czero 
       endif
       k2 = -k2 - p 
       k2sq = sc(k2,k2) !- mass**2

       sp2 = spi2(k2,sp2)!+mass*sp2

       tmp = vWq(eW,sp2)

       if (abs(k2sq) > propcut) then 
          tmp = ci/k2sq*tmp
       else 
          tmp = czero
       endif

       res = res + tmp

    endif

    ! -- store current 
!    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)

  end function bfW




  recursive function gW_fbf(e,k,sp1,p1,fl1,sp2,p2,fl2,&
       &eW,kW,ng1,ng2,giarray,qiarray,Wid,pol_int) result(res)
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp1(:), p1(:), sp2(:), p2(:)
    double complex, intent(in) :: eW(:),kW(:)
    integer, intent(in) ::  ng1,ng2
    character, intent(in) :: fl1*3,fl2*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
    ! -----------------------------------------------------------------------
    integer             :: m1,m2, ms1a, ms2a,Dv!,m,m3
    integer :: ngluon,  ng3, ngL 
    double complex             :: res(size(e,dim=1))
    double complex             :: tmp(size(e,dim=1))
    double complex             :: k1(size(p1))
    double complex             :: k2(size(p1))
    double complex             :: k3(size(p1))
    double complex             :: k4(size(p1))
    double complex             :: sp3(size(sp1))
    double complex             :: sp4(size(sp1))
    double complex             :: e1(size(e,dim=1))
    double complex             :: e2(size(e,dim=1))
    double complex             :: e3(size(e,dim=1))
    double complex  :: k1sq,k2sq,k3sq,k4sq
    !real(dp) :: mass,mass2
    logical                   :: done 

    !if (verbose) write(*,*) 'entering gW_fbf'

    if ((fl1.eq.'str'.or.fl2.eq.'str') .or. (fl1.eq.'dwn'.or.fl2.eq.'dwn')) then 
       res = czero 
       return 
    endif

    if (fl1.eq.fl2) stop 'gw_fbf: flav'

    done = .false. 
!    if (present(giarray)) then 
!       !if (size(qiarray) /= 2) stop 'gW_fbf: wrong size qiarray'
!       !if (size(e,dim=2) /= size(giarray)) stop 'gW_fbf: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!       if (done) return 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'gw_fbf: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif

    !mass = mt
    !mass2 = mass**2

    Dv = size(e,dim=1)

    ngluon = size(e,dim=2)
    ng3 = ngluon - ng1 - ng2


    if (ng3 < 0) write(*,*) 'ERROR IN CURRENT D'

    if (ngluon == 0) then 

       res = czero

       sp4 = fW(e,k,sp2,p2,fl2,fl1,eW,kW,0,giarray,qiarray(2:2),Wid,pol_int)    
       k2 = p2 + kW
       k2sq = sc(k2,k2) !-mass2
       sp4 = ci/k2sq*(spb2(sp4,k2))!+mass*sp4)
       res = -cone*vbqq(Dv,sp4,sp1)

       sp4 = bfW(e,k,sp1,p1,fl1,fl2,eW,kW,0,giarray,qiarray(1:1),Wid,pol_int)    
       k2 = -p1 - kW
       k2sq = sc(k2,k2) !- mass2
       sp4 = ci/k2sq*(spi2(k2,sp4) )!+ mass*sp4)
       tmp = -cone*vbqq(Dv,sp2,sp4)

       res = res + tmp

    else

       res = czero

       do m1=1,ng1

          e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)    
          k1=sum(k(:,1:m1),dim=2)
          k1sq = sc(k1,k1)

          ms1a = ng1-m1
          e2=gW_fbf(e(:,m1+1:ngluon),k(:,m1+1:ngluon),&
               &sp1,p1,fl1,sp2,p2,fl2,eW,kW,ms1a,ng2,&
               &giarray(m1+1:ngluon),qiarray(1:2),Wid,pol_int)    
          if (m1+1<=ngluon) then 
             k2 = sum(k(:,m1+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2 = k2 + p1 + p2+kW
          k2sq = sc(k2,k2)

          if (abs(k2sq) > propcut) then 
             tmp = -ci/k2sq*vggg(e1,k1,e2,k2)
          else 
             tmp = czero 
          endif

          if (m1 > 1) then  
             if (abs(k1sq) > propcut) then 
                tmp = -ci/k1sq*tmp
             else 
                tmp = czero 
             endif
          endif

          res = res + tmp

       enddo  !#1


       do m1=0,ng3-1 

          e1 = gW_fbf(e(:,1:ng1+ng2+m1),k(:,1:ng1+ng2+m1),&
               &sp1,p1,fl1,sp2,p2,fl2,eW,kW,ng1,ng2,&
               &giarray(1:ng1+ng2+m1),qiarray,Wid,pol_int)    
          if (1<=ng1+ng2+m1) then 
             k1=sum(k(:,1:ng1+ng2+m1),dim=2)
          else
             k1 = czero 
          endif
          k1 = k1 + p1 + p2+kW
          k1sq=sc(k1,k1)

          e2=vgluon(e(:,ng1+ng2+m1+1:ngluon),&
               &k(:,ng1+ng2+m1+1:ngluon),giarray(ng1+ng2+m1+1:ngluon),pol_int)    
          if (ng1+ng2+m1+1<=ngluon) then 
             k2=sum(k(:,ng1+ng2+m1+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2sq = sc(k2,k2)

          if (abs(k1sq) > propcut) then 
             tmp = -ci/k1sq*vggg(e1,k1,e2,k2)
          else
             tmp = czero 
          endif

          if (m1 < ng3-1) then 
             if (abs(k2sq) > propcut) then 
                tmp = -ci/k2sq*tmp
             else 
                tmp = czero 
             endif
          endif

          res = res + tmp

       enddo   !#2


       do m1 =1, ng1-1

          e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)    
          k1=sum(k(:,1:m1),dim=2)
          k1sq=sc(k1,k1)

          do m2 = m1+1,ng1   

             e2=vgluon(e(:,m1+1:m2),k(:,m1+1:m2),giarray(m1+1:m2),pol_int)    
             if (m1+1<=m2) then 
                k2=sum(k(:,m1+1:m2),dim=2)
             else
                k2 = czero 
             endif
             k2sq=sc(k2,k2)


             ms1a=ng1-m2
             e3=gW_fbf(e(:,m2+1:ngluon),k(:,m2+1:ngluon),&
                  &sp1,p1,fl1,sp2,p2,fl2,eW,kW,ms1a,ng2,&
                  &giarray(m2+1:ngluon),qiarray,Wid,pol_int)    
             if (m2+1<=ngluon) then 
                k3=sum(k(:,m2+1:ngluon),dim=2)
             else
                k3 = czero 
             endif
             k3 = k3 + p1 + p2+kW
             k3sq=sc(k3,k3)

             if (abs(k3sq) > propcut) then 
                tmp = -ci/k3sq*vgggg(e1,e2,e3)
             else 
                tmp = czero 
             endif

             if (m1>1) then 
                if (abs(k1sq) > propcut) then 
                   tmp = -ci/k1sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             if (m2 > m1+1) then 
                if (abs(k2sq) > propcut) then 
                   tmp = -ci/k2sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             res = res + tmp

          enddo

       enddo   !#3


       do m1=1,ng1

          e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)    
          k1=sum(k(:,1:m1),dim=2)
          k1sq=sc(k1,k1)

          do m2 = 0,ng3-1    

             ms1a=ng1-m1
             e2=gW_fbf(e(:,m1+1:ng1+ng2+m2),k(:,m1+1:ng1+ng2+m2),&
                  &sp1,p1,fl1,sp2,p2,fl2,eW,kW,ms1a,ng2,&
                  &giarray(m1+1:ng1+ng2+m2),qiarray,Wid,pol_int)    
             if (m1+1<=ng1+ng2+m2) then 
                k2=sum(k(:,m1+1:ng1+ng2+m2),dim=2)
             else
                k2 = czero 
             endif
             k2 = k2 + p1 + p2 + kW
             k2sq=sc(k2,k2)


             e3=vgluon(e(:,ng1+ng2+m2+1:ngluon)&
                  &,k(:,ng1+ng2+m2+1:ngluon),&
                  &giarray(ng1+ng2+m2+1:ngluon),pol_int)    
             if (ng1+ng2+m2+1<=ngluon) then 
                k3=sum(k(:,ng1+ng2+m2+1:ngluon),dim=2)
             else
                k3 = czero 
             endif
             k3sq=sc(k3,k3)

             if (abs(k2sq) > propcut) then 
                tmp = -ci/k2sq*vgggg(e1,e2,e3)
             else 
                tmp = czero 
             endif

             if (m1>1) then 
                if (abs(k1sq) > propcut) then 
                   tmp = -ci/k1sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             if (m2+1+ng1+ng2 < ngluon) then 
                if (abs(k3sq) > propcut) then 
                   tmp = -ci/k3sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             res = res + tmp

          enddo

       enddo  !#4



       do m1=0,ng3-2

          ngL=ng1+ng2+m1

          e1=gW_fbf(e(:,1:ngL),k(:,1:ngL),&
               &sp1,p1,fl1,sp2,p2,fl2,eW,kW,ng1,ng2,&
               &giarray(1:ngL),qiarray,Wid,pol_int)    
          if (1<=ngL) then 
             k1=sum(k(:,1:ngL),dim=2)
          else
             k1 = czero 
          endif
          k1 = k1 + p1 + p2+kW
          k1sq=sc(k1,k1)

          do m2=ngL+1,ngluon-1

             e2=vgluon(e(:,ngL+1:m2),&
                  &k(:,ngL+1:m2),giarray(ngL+1:m2),pol_int)    
             if (ngL+1<=m2) then 
                k2=sum(k(:,ngL+1:m2),dim=2)
             else
                k2 = czero 
             endif
             k2sq=sc(k2,k2)

             e3=vgluon(e(:,m2+1:ngluon),k(:,m2+1:ngluon),&
                  &giarray(m2+1:ngluon),pol_int)    
             if (m2+1<=ngluon) then 
                k3=sum(k(:,m2+1:ngluon),dim=2)
             else
                k3 = czero 
             endif
             k3sq=sc(k3,k3)


             if (abs(k1sq) > propcut) then 
                tmp = -ci/k1sq*vgggg(e1,e2,e3)
             else 
                tmp = czero 
             endif

             if (m2 > ngL+1) then 
                if (abs(k2sq)> propcut) then 
                   tmp=-ci/k2sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             if (m2 < ng3-1) then 
                if (abs(k3sq) > propcut) then  
                   tmp=-ci/k3sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             res = res + tmp

          enddo

       enddo  !#5


       do m1=0,ng2

          ms1a=m1+ng1   
          sp3=bf(e(:,1:ms1a),k(:,1:ms1a),sp1,p1,fl1,fl1,ng1,&
               &giarray(1:ms1a),qiarray(1:1),pol_int)    
          if (1<=ms1a) then 
             k3=sum(k(:,1:ms1a),dim=2)
          else
             k3 = czero 
          endif
          k3 = -k3 - p1
          k3sq = sc(k3,k3)!-mass2

          if (ng1 > 0.or.m1 > 0) sp3 = spi2(k3,sp3)!+mass*sp3

          ms2a=ng2-m1
          sp4=fW(e(:,ms1a+1:ngluon),k(:,ms1a+1:ngluon),&
               &sp2,p2,fl2,fl1,eW,kW,ms2a,&
               &giarray(ms1a+1:ngluon),qiarray(2:2),Wid,pol_int)    
          if (ms1a+1<=ngluon) then 
             k4=sum(k(:,ms1a+1:ngluon),dim=2)
          else
             k4 = czero 
          endif
          k4 =  k4 + p2+kW
          k4sq = sc(k4,k4)!-mass2
          sp4 = spb2(sp4,k4)!+mass*sp4

          tmp = -cone*vbqq(Dv,sp4,sp3)

          if (ng1 > 0.or.m1 > 0) then 
             if (abs(k3sq) > propcut) then 
                tmp= ci/k3sq*tmp
             else 
                tmp = czero 
             endif
          endif


          if (abs(k4sq) > propcut) then 
             tmp= ci/k4sq*tmp 
          else 
             tmp = czero 
          endif


          res = res + tmp

       enddo  !#6


       do m1=0,ng2

          ms1a=m1+ng1   
          sp3=bfW(e(:,1:ms1a),k(:,1:ms1a),sp1,p1,fl1,fl2,eW,kW,ng1,&
               &giarray(1:ms1a),qiarray(1:1),Wid,pol_int)    
          if (1<=ms1a) then 
             k3=sum(k(:,1:ms1a),dim=2)
          else
             k3 = czero 
          endif
          k3 = -k3 - p1 - kW
          k3sq = sc(k3,k3)!-mass2
          sp3 = spi2(k3,sp3)!+mass*sp3

          ms2a=ng2-m1
          sp4=f(e(:,ms1a+1:ngluon),k(:,ms1a+1:ngluon),&
               &sp2,p2,fl2,fl2,ms2a,&
               &giarray(ms1a+1:ngluon),qiarray(2:2),pol_int)    
          if (ms1a+1<=ngluon) then 
             k4=sum(k(:,ms1a+1:ngluon),dim=2)
          else
             k4 = czero 
          endif
          k4 =  k4 + p2
          k4sq = sc(k4,k4)!-mass2

          if (ng3 > 0.or.ng2-m1>0) sp4 = spb2(sp4,k4)!+mass*sp4

          tmp = -cone*vbqq(Dv,sp4,sp3)

          if (abs(k3sq) > propcut) then 
             tmp= ci/k3sq*tmp
          else 
             tmp = czero 
          endif


          if (ng3 > 0.or.ng2-m1>0) then 
             if (abs(k4sq) > propcut) then 
                tmp= ci/k4sq*tmp 
             else 
                tmp = czero 
             endif
          endif

          res = res + tmp

       enddo  !#7


    endif   !end condition for ngluon 



    ! -- store current 
!    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)



  end function gW_fbf



  recursive function gW_bff(e,k,sp1,p1,fl1,&
       &sp2,p2,fl2,eW,kW,ms1,ms2,giarray,qiarray,Wid,pol_int) result(res)
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp1(:), p1(:), sp2(:), p2(:)
    double complex, intent(in) :: eW(:), kW(:)
    integer, intent(in) ::  ms1,ms2
    character, intent(in) :: fl1*3, fl2*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
    ! -----------------------------------------------------------------------
    integer, parameter :: Ndumm=0
    integer             :: m1,m2, ms1a, ms2a,Dv!,m,m3
    integer :: ngluon, ng1, ng2, ng3, ngL
    double complex             :: res(size(e,dim=1))
    double complex             :: tmp(size(e,dim=1))
    double complex             :: k1(size(p1))
    double complex             :: k2(size(p1))
    double complex             :: k3(size(p1))
    double complex             :: k4(size(p1))
    double complex             :: sp3(size(sp1))
    double complex             :: sp4(size(sp1))
    double complex             :: e1(size(e,dim=1))
    double complex             :: e2(size(e,dim=1))
    double complex             :: e3(size(e,dim=1))
    double complex  :: k1sq,k2sq,k3sq,k4sq
    double complex             :: kdumm(size(k,dim=1),Ndumm)
    double complex             :: edumm(size(e,dim=1),Ndumm)
    !real(dp) :: mass,mass2
    logical                   :: done 

    !if (verbose) write(*,*) 'entering gW_bff'

    ! added 06/07/07
    if ((fl1.eq.'str'.or.fl2.eq.'str') .or. (fl1.eq.'dwn'.or.fl2.eq.'dwn')) then 
       res = czero 
       return 
    endif

    done = .false. 
!    if (present(giarray)) then 
!       !if (size(qiarray) /= 2) stop 'gW_bff: wrong size qiarray'
!       !if (size(e,dim=2) /= size(giarray)) stop 'gW_bff: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!       if (done) return 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'gW_bff: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif

    !mass = mt
    !mass2 = mass**2

    if (size(sp1) == 4) then 
       Dv=4
    elseif (size(sp1) == 8) then 
       Dv=6
    elseif (size(sp1) == 16) then
       Dv=8
    else
       stop 'Dv undefined' 
    endif

    ngluon = size(e,dim=2)
    ng1 = ms1
    ng2 = ms2
    ng3 = ngluon - ms1 - ms2

    if (ng3 < 0) write(*,*) 'ERROR IN CURRENT A'



    if (ngluon == 0) then 

       res = czero

       sp3=fW(edumm,kdumm,sp1,p1,fl1,fl2,eW,kW,0,giarray,qiarray(1:1),Wid,pol_int)    
       k3 = p1+kW
       k3sq = sc(k3,k3)!-mass**2

       sp3 = spb2(sp3,k3)!+mass*sp3


       sp4=bf(edumm,kdumm,&
            &sp2,p2,fl2,fl2,0,giarray,qiarray(2:2),pol_int)    


       tmp = vbqq(Dv,sp3,sp4)

       if (abs(k3sq) > propcut) then 
          tmp= ci/k3sq*tmp
       else
          tmp = czero 
       endif

       res = res + tmp



       sp3=f(edumm,kdumm,sp1,p1,fl1,fl1,0,giarray,qiarray(1:1),pol_int)    


       sp4=bfW(edumm,kdumm,&
            &sp2,p2,fl2,fl1,eW,kW,0,giarray,qiarray(2:2),Wid,pol_int)    
       k4 = - p2-kW
       k4sq = sc(k4,k4)!-mass**2

       sp4 = spi2(k4,sp4)!+mass*sp4

       tmp = vbqq(Dv,sp3,sp4)

       if (abs(k4sq) > propcut) then 
          tmp= ci/k4sq*tmp 
       else
          tmp = czero 
       endif

       res = res + tmp

    else
       res = czero


       do m1=1,ng1

          e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)    
          if (1<=m1) then 
             k1=sum(k(:,1:m1),dim=2)
          else
             k1 = czero 
          endif
          k1sq = sc(k1,k1)

          ms1a = ng1-m1
          e2=gW_bff(e(:,m1+1:ngluon),k(:,m1+1:ngluon),&
               &sp1,p1,fl1,sp2,p2,fl2,eW,kW,ms1a,ng2,&
               &giarray(m1+1:ngluon),qiarray,Wid,pol_int)    
          if (m1+1<=ngluon) then 
             k2 = sum(k(:,m1+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2 = k2 + p1 + p2+kW
          k2sq = sc(k2,k2)

          if (abs(k2sq) > propcut) then 
             tmp = -ci/k2sq*vggg(e1,k1,e2,k2)
          else
             tmp = czero 
          endif

          if (m1 > 1) then 
             if (abs(k1sq) > propcut) then 
                tmp = -ci/k1sq*tmp
             else
                tmp = czero 
             endif
          endif



          res = res + tmp

       enddo  !#1


       do m1=0,ng3-1 

          e1 = gW_bff(e(:,1:ng1+ng2+m1),k(:,1:ng1+ng2+m1),&
               &sp1,p1,fl1,sp2,p2,fl2,eW,kW,ng1,ng2,&
               &giarray(1:ng1+ng2+m1),qiarray,Wid,pol_int)    
          if (1<=ng1+ng2+m1) then 
             k1=sum(k(:,1:ng1+ng2+m1),dim=2)
          else
             k1 = czero 
          endif
          k1 = k1 + p1 + p2+kW
          k1sq=sc(k1,k1)

          e2=vgluon(e(:,ng1+ng2+m1+1:ngluon),&
               &k(:,ng1+ng2+m1+1:ngluon),giarray(ng1+ng2+m1+1:ngluon),pol_int)    
          if (ng1+ng2+m1+1<=ngluon) then 
             k2=sum(k(:,ng1+ng2+m1+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2sq = sc(k2,k2)

          if (abs(k1sq) > propcut) then 
             tmp = -ci/k1sq*vggg(e1,k1,e2,k2)
          else 
             tmp = czero 
          endif


          if (m1 < ng3-1) then 
             if (abs(k2sq) > propcut) then 
                tmp = -ci/k2sq*tmp
             else
                tmp = czero 
             endif
          endif

          res = res + tmp

       enddo   !#2


       do m1 =1, ng1-1

          e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)    
          k1=sum(k(:,1:m1),dim=2)
          k1sq=sc(k1,k1)

          do m2 = m1+1,ng1   

             e2=vgluon(e(:,m1+1:m2),k(:,m1+1:m2),giarray(m1+1:m2),pol_int)    
             if (m1+1<=m2) then 
                k2=sum(k(:,m1+1:m2),dim=2)
             else
                k2 = czero 
             endif
             k2sq=sc(k2,k2)


             ms1a=ng1-m2
             e3=gW_bff(e(:,m2+1:ngluon),k(:,m2+1:ngluon),&
                  &sp1,p1,fl1,sp2,p2,fl2,eW,kW,ms1a,ng2,&
                  &giarray(m2+1:ngluon),qiarray,Wid,pol_int)    
             if (m2+1<=ngluon) then 
                k3=sum(k(:,m2+1:ngluon),dim=2)
             else
                k2 = czero 
             endif
             k3 = k3 + p1 + p2 + kW
             k3sq=sc(k3,k3)

             if (abs(k3sq) > propcut) then 
                tmp = -ci/k3sq*vgggg(e1,e2,e3)
             else
                tmp = czero 
             endif


             if (m1>1)   then    
                if (abs(k1sq) > propcut) then 
                   tmp = -ci/k1sq*tmp
                else
                   tmp = czero 
                endif
             endif

             if (m2 > m1+1)  then 
                if (abs(k2sq) > propcut) then 
                   tmp = -ci/k2sq*tmp
                else
                   tmp = czero 
                endif
             endif

             res = res + tmp

          enddo

       enddo   !#3


       do m1=1,ng1

          e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)    
          k1=sum(k(:,1:m1),dim=2)
          k1sq=sc(k1,k1)

          do m2 = 0,ng3-1    

             ms1a=ng1-m1
             e2=gW_bff(e(:,m1+1:ng1+ng2+m2),k(:,m1+1:ng1+ng2+m2),&
                  &sp1,p1,fl1,sp2,p2,fl2,eW,kW,ms1a,ng2,&
                  &giarray(m1+1:ng1+ng2+m2),qiarray,Wid,pol_int)    
             if (m1+1<=ng1+ng2+m2) then 
                k2=sum(k(:,m1+1:ng1+ng2+m2),dim=2)
             else
                k2 = czero 
             endif 
             k2 = k2 + p1 + p2+kW
             k2sq=sc(k2,k2)


             e3=vgluon(e(:,ng1+ng2+m2+1:ngluon)&
                  &,k(:,ng1+ng2+m2+1:ngluon),&
                  &giarray(ng1+ng2+m2+1:ngluon),pol_int)    
             if (ng1+ng2+m2+1<=ngluon) then 
                k3=sum(k(:,ng1+ng2+m2+1:ngluon),dim=2)
             else
                k3 = czero 
             endif
             k3sq=sc(k3,k3)


             if (abs(k2sq) > propcut) then 
                tmp = -ci/k2sq*vgggg(e1,e2,e3)
             else 
                tmp = czero 
             endif

             if (m1>1) then  
                if (abs(k1sq) > propcut) then 
                   tmp = -ci/k1sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             if (m2+1+ng1+ng2 < ngluon) then  
                if (abs(k3sq) > propcut) then 
                   tmp = -ci/k3sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             res = res + tmp



          enddo

       enddo  !#4



       do m1=0,ng3-2

          ngL=ng1+ng2+m1

          e1=gW_bff(e(:,1:ngL),k(:,1:ngL),&
               &sp1,p1,fl1,sp2,p2,fl2,eW,kW,ng1,ng2,&
               &giarray(1:ngL),qiarray,Wid,pol_int)    
          if (1<=ngL) then 
             k1=sum(k(:,1:ngL),dim=2)
          else
             k1 = czero 
          endif
          k1 = k1 + p1 + p2
          k1sq=sc(k1,k1)

          do m2=ngL+1,ngluon-1

             e2=vgluon(e(:,ngL+1:m2),&
                  &k(:,ngL+1:m2),giarray(ngL+1:m2),pol_int)    
             if (ngL+1<=m2) then 
                k2=sum(k(:,ngL+1:m2),dim=2)
             else
                k2 = czero 
             endif
             k2sq=sc(k2,k2)

             e3=vgluon(e(:,m2+1:ngluon),k(:,m2+1:ngluon),giarray(m2+1:ngluon),pol_int)  
             if (m2+1<=ngluon) then 
                k3=sum(k(:,m2+1:ngluon),dim=2)
             else
                k3 = czero 
             endif
             k3sq=sc(k3,k3)

             if (abs(k1sq) > propcut) then 
                tmp = -ci/k1sq*vgggg(e1,e2,e3)
             else 
                tmp = czero 
             endif

             if (m2 > ngL+1) then
                if (abs(k2sq) > propcut) then 
                   tmp=-ci/k2sq*tmp
                else
                   tmp = czero 
                endif
             endif

             if (m2 < ng3-1) then 
                if (abs(k3sq) > propcut) then 
                   tmp=-ci/k3sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             res = res + tmp



          enddo

       enddo  !#5


       do m1=0,ng2

          ms1a=m1+ng1   
          sp3=fW(e(:,1:ms1a),k(:,1:ms1a),sp1,p1,fl1,fl2,eW,kW,ng1,&
               &giarray(1:ms1a),qiarray(1:1),Wid,pol_int)    
          if (1<=ms1a) then 
             k3=sum(k(:,1:ms1a),dim=2)
          else
             k3 = czero 
          endif
          k3 = k3 + p1+kW
          k3sq = sc(k3,k3)!-mass**2

          sp3 = spb2(sp3,k3)!+mass*sp3

          ms2a=ng2-m1
          sp4=bf(e(:,ms1a+1:ngluon),k(:,ms1a+1:ngluon),&
               &sp2,p2,fl2,fl2,ms2a,&
               &giarray(ms1a+1:ngluon),qiarray(2:2),pol_int)    
          if (ms1a+1<=ngluon) then 
             k4=sum(k(:,ms1a+1:ngluon),dim=2)
          else
             k4 = czero
          endif 
          k4 = - k4 - p2
          k4sq = sc(k4,k4)!-mass**2

          if (ng3 > 0.or.ng2-m1>0) sp4 = spi2(k4,sp4)!+mass*sp4

          tmp = vbqq(Dv,sp3,sp4)

          if (abs(k3sq) > propcut) then 
             tmp= ci/k3sq*tmp
          else
             tmp = czero 
          endif

          if (ng3 > 0.or.ng2-m1>0) then  
             if (abs(k4sq) > propcut) then 
                tmp= ci/k4sq*tmp 
             else
                tmp = czero 
             endif
          endif


          res = res + tmp

       enddo  !#6


       do m1=0,ng2

          ms1a=m1+ng1   
          sp3=f(e(:,1:ms1a),k(:,1:ms1a),sp1,p1,fl1,fl1,ng1,&
               &giarray(1:ms1a),qiarray(1:1),pol_int)    
          if (1<=ms1a) then 
             k3=sum(k(:,1:ms1a),dim=2)
          else
             k3 = czero 
          endif
          k3 = k3 + p1
          k3sq = sc(k3,k3)!-mass**2

          if (ng1 > 0.or.m1 > 0) sp3 = spb2(sp3,k3)!+mass*sp3

          ms2a=ng2-m1
          sp4=bfW(e(:,ms1a+1:ngluon),k(:,ms1a+1:ngluon),&
               &sp2,p2,fl2,fl1,eW,kW,ms2a,&
               &giarray(ms1a+1:ngluon),qiarray(2:2),Wid,pol_int)    
          if (ms1a+1<=ngluon) then 
             k4=sum(k(:,ms1a+1:ngluon),dim=2)
          else
             k4 = czero 
          endif
          k4 = - k4 - p2-kW
          k4sq = sc(k4,k4)!-mass**2

          sp4 = spi2(k4,sp4)!+mass*sp4

          tmp = vbqq(Dv,sp3,sp4)


          if (ng1 > 0.or.m1 > 0) then  
             if (abs(k3sq) > propcut) then 
                tmp= ci/k3sq*tmp
             else
                tmp = czero 
             endif
          endif


          if (abs(k4sq) > propcut) then 
             tmp= ci/k4sq*tmp 
          else
             tmp = czero 
          endif

          res = res + tmp

       enddo  !#7



    endif   !end condition for ngluon 

    ! -- store current 
!    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)


  end function gW_bff



  recursive function fW_bffbf(e,k,sp,p,fll,fl0,&
       &eW,kW,ng1,ng2,ng3,sw,giarray,qiarray,Wid,pol_int) result(res)
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp(:,:),p(:,:)
    double complex, intent(in) :: eW(:),kW(:)
    integer, intent(in) ::  ng1,ng2,ng3,sw
    character, intent(in) :: fll(:)*3
    character, intent(in) :: fl0*3   ! flavor off-shell f-line
    integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
    ! -----------------------------------------------------------------------
    character :: flaux*3
    character :: fl1*3,fl2*3,fl3*3
    integer :: ngluon, ng4, ngL,m!,m1
    integer, parameter :: Ndumm=0
    double complex             :: res(size(sp,dim=1))
    double complex             :: tmp(size(sp,dim=1))
    double complex             :: k1(size(k,dim=1))
    double complex             :: k2(size(k,dim=1))
    double complex             :: k4(size(k,dim=1))
    double complex             :: sp1(size(sp,dim=1))
    double complex             :: sp2(size(sp,dim=1))
    double complex             :: sp4(size(sp,dim=1))
    double complex             :: e1(size(e,dim=1))
    double complex             :: e2(size(e,dim=1))
    double complex             :: kdumm(size(k,dim=1),Ndumm)
    double complex             :: edumm(size(e,dim=1),Ndumm)
    double complex             :: k1sq,k2sq,k4sq
    !real(dp)                :: mass,mass2
    logical                   :: done 

    !if (verbose)    write(*,*) 'entering fW_bffbf',ng1,ng2,ng3

    done = .false. 
!    if (present(giarray)) then 
!       !if (size(qiarray) /= 3) stop 'fW_bffbf: wrong size qiarray'
!       !if (size(e,dim=2) /= size(giarray)) stop 'fW_bffbf: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!       if (done) return 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'fW_bffbf: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif

    !mass = mt
    !mass2 = mass**2

    ngluon = size(e,dim=2)
    ng4 = ngluon - ng1 - ng2-ng3

    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)


    if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C'

    if (ngluon == 0) then 

       res = czero

       if (sw.eq.3) then 

          e2 = gW_fbf(edumm,kdumm,sp(:,2),p(:,2),fl2,&
               &sp(:,3),p(:,3),fl3,eW,kW,0,0,&
               &giarray,qiarray(2:3),Wid,pol_int)    
          k1 = p(:,2)+p(:,3)+kW
          k1sq=sc(k1,k1)

          if (abs(k1sq) > propcut.and.fl0==fl1) then 
             tmp = -ci/k1sq*vqg(sp(:,1),e2)
          else 
             tmp = czero 
          endif

          res = res + tmp  ! #1

       endif

! XXXXX
       if (sw.eq.2 .and. ((case_b2 .eqv. .true.) .or. &
            &(qbq_WW_and_gluons .eqv. .false.))) then 

          e2 = gW_bff(edumm,kdumm,sp(:,1),p(:,1),fl1,&
               &sp(:,2),p(:,2),fl2,eW,kW,0,0,&
               &giarray,qiarray(1:2),Wid,pol_int)    
          k2 = p(:,1)+p(:,2)+kW 
          k2sq = sc(k2,k2)

          if (abs(k2sq) > propcut.and.fl0==fl3) then 
             tmp = -ci/k2sq*vgq(e2,sp(:,3))
          else 
             tmp = czero 
          endif

          !        used to have czero here (now)
          res = res + tmp   ! #2

       endif


       if ((sw.eq.1).or.(sw.eq.2)) then 

          e2 = g_fbf(edumm,kdumm,sp(:,2),p(:,2),fl2,&
               &sp(:,3),p(:,3),fl3,0,0,&
               &giarray,qiarray(2:3),pol_int)    

          k2 = p(:,2)+p(:,3)
          k2sq = sc(k2,k2)

          sp4 = fW(edumm,kdumm,sp(:,1),p(:,1),fl1,fl0,eW,kW,0,&
               &giarray,qiarray(1:1),Wid,pol_int)    
          k4  = p(:,1) + kW
          k4sq = sc(k4,k4)
          sp4 = spb2(sp4,k4) !+ mass*sp4

          tmp = vqg(sp4,e2)

          if (abs(k2sq) > propcut) then 
             tmp = -ci/k2sq*tmp
          else 
             tmp = czero 
          endif

          if (abs(k4sq) > propcut) then 
             tmp = ci/k4sq*tmp 
          else 
             tmp = czero 
          endif

          res = res + tmp  ! #3

       endif

!XXXXX
       if (((sw.eq.3).or.(sw.eq.4)) .and. ((qbq_WW_and_gluons .eqv. .false.) .or. (case_a3 .eqv. .true.))) then 
          
          e2 = g_bff(edumm,kdumm,sp(:,1),p(:,1),fl1,&
               & sp(:,2),p(:,2),fl2,0,0,&
               &giarray,qiarray(1:2),pol_int)    
          k2 = p(:,1)+p(:,2)
          k2sq = sc(k2,k2)


          sp4 = fW(edumm,kdumm,sp(:,3),p(:,3),fl3,fl0,eW,kW,0,&
               &giarray,qiarray(3:3),Wid,pol_int)    
          k4  = p(:,3) + kW
          k4sq = sc(k4,k4)
          sp4 = spb2(sp4,k4) !+ mass*sp4

          tmp = vgq(e2,sp4)   

          if (abs(k2sq) > propcut) then 
             tmp = -ci/k2sq*tmp
          else 
             tmp = czero 
          endif

          if (abs(k4sq) > propcut) then 
             tmp = ci/k4sq*tmp 
          else 
             tmp = czero 
          endif

          !---------used to have czero here 


          res = res + tmp  !# 4
!XXX       endif

       endif

       if ((sw.eq.1).or.(sw.eq.4)) then 

          if (fl0.eq.'top') flaux = 'bot'
          if (fl0.eq.'bot') flaux = 'top'
          if (fl0.eq.'chr') flaux = 'chr'

! XXXX
          if ((qbq_WW_and_gluons .eqv. .true.) .and. (case_a3 .eqv. .false.)) then
             sp4 = f_bffbf_2(edumm,kdumm,sp,p,fll,flaux,&
                  &ng1,ng2,ng3,&
                  &giarray,qiarray,pol_int)   
          else
             sp4 = f_bffbf(edumm,kdumm,sp,p,fll,flaux,&
                  &ng1,ng2,ng3,&
                  &giarray,qiarray,pol_int)    

          endif

          k4 = p(:,1)+p(:,2)+p(:,3)
          sp4 = spb2(sp4,k4) !+ mass*sp4
          k4sq = sc(k4,k4)

          tmp = vbqW(sp4,eW)

          if (abs(k4sq) > propcut) then 
             tmp = ci/k4sq*tmp
          else 
             tmp = czero
          endif



          res = res + tmp   ! #5

       endif

    else  ! -- this is for ngluon > 0

       res = czero

       if (sw.eq.1.or.sw.eq.2) then 

          do m=0,ng2

             sp1=fW(e(:,1:ng1+m),k(:,1:ng1+m),sp(:,1),p(:,1),&
                  &fl1,fl0,eW,kW,ng1,&
                  &giarray(1:ng1+m),qiarray(1:1),Wid,pol_int)    
             if (1<=ng1+m) then 
                k1=sum(k(:,1:ng1+m),dim=2)
             else
                k1 = czero 
             endif
             k1 = k1 + p(:,1)+kW
             k1sq = sc(k1,k1) !- mass2

             sp1 = spb2(sp1,k1)!+mass*sp1

             e2 = g_fbf(e(:,ng1+m+1:ngluon),k(:,ng1+m+1:ngluon),&
                  &sp(:,2),p(:,2),fl2,sp(:,3),p(:,3),fl3,ng2-m,ng3,&
                  &giarray(ng1+m+1:ngluon),qiarray(2:3),pol_int)    
             if (ng1+m+1<=ngluon) then 
                k2=sum(k(:,ng1+m+1:ngluon),dim=2)
             else
                k2 = czero 
             endif
             k2 = k2 + p(:,2)+p(:,3)
             k2sq=sc(k2,k2)


             if (abs(k2sq) > propcut) then 
                tmp = -ci/k2sq*vqg(sp1,e2)
             else 
                tmp = czero 
             endif


             if (abs(k1sq) > propcut) then 
                tmp = ci/k1sq*tmp
             else 
                tmp = czero 
             endif

             res = res + tmp

          enddo  !#1

       endif


       !--------- next step is valid for all sw

       do m=0,ng4-1

          ngL = ng1+ ng2+ng3+m      

          sp1=fW_bffbf(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0,&
               &eW,kW,ng1,ng2,ng3,sw,&
               &giarray(1:ngL),qiarray,Wid,pol_int)    
          if (1<=ngL) then 
             k1=sum(k(:,1:ngL),dim=2)
          else
             k1 = czero 
          endif
          k1 = k1 + p(:,1)+p(:,2)+p(:,3)+kW
          k1sq = sc(k1,k1) !- mass2

          sp1 = spb2(sp1,k1) !+ mass*sp1

          e2 = vgluon(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
               &giarray(ngL+1:ngluon),pol_int)    
          if (ngL+1<=ngluon) then 
             k2=sum(k(:,ngL+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2sq=sc(k2,k2)

          if (abs(k1sq) > propcut) then 
             tmp = ci/k1sq*vqg(sp1,e2)
          else 
             tmp = czero 
          endif

          if (m < ng4-1) then 
             if (abs(k2sq) > propcut) then 
                tmp = -ci/k2sq*tmp
             else 
                tmp = czero 
             endif
          endif

          res = res + tmp


       enddo  !#2


       !-------- next step is valid for any sw

       do m=1,ng1

          e1 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)    
          k1=sum(k(:,1:m),dim=2)
          k1sq=sc(k1,k1)

          sp2=fW_bffbf(e(:,m+1:ngluon),k(:,m+1:ngluon),&
               &sp,p,fll,fl0,eW,kW,ng1-m,ng2,ng3,sw,&
               &giarray(m+1:ngluon),qiarray,Wid,pol_int)    
          if (m+1<=ngluon) then 
             k2=sum(k(:,m+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2 = k2 + p(:,1)+p(:,2)+p(:,3)+kW
          k2sq = sc(k2,k2) !- mass2
          sp2 = spb2(sp2,k2)!+mass*sp2


          if (abs(k2sq) > propcut) then 
             tmp = ci/k2sq*vgq(e1,sp2)
          else 
             tmp = czero 
          endif


          if (m > 1) then 
             if (abs(k1sq) > propcut) then 
                tmp = -ci/k1sq*tmp
             else 
                tmp = czero 
             endif
          endif

          res = res + tmp


       enddo  !#3

! XXXXXX
       if ((sw.eq.3.or.sw.eq.4) .and. (WWqqqq .eqv. .false.)) then 

          do m=0,ng3

             ngL = ng1+ ng2+m      

             e1 = g_bff(e(:,1:ngL),k(:,1:ngL),sp(:,1),p(:,1),fl1,&
                  &sp(:,2),p(:,2),fl2,ng1,ng2,&
                  &giarray(1:ngL),qiarray(1:2),pol_int)    
             if (1<=ngL) then 
                k1=sum(k(:,1:ngL),dim=2)+p(:,1)+p(:,2)
             else
                k1 = p(:,1)+p(:,2)
             endif
             k1sq=sc(k1,k1)


             sp2=fW(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
                  &sp(:,3),p(:,3),fl3,fl0,eW,kW,ng3-m,&
                  &giarray(ngL+1:ngluon),qiarray(3:3),Wid,pol_int)    
             if (ngL+1<=ngluon) then 
                k2=sum(k(:,ngL+1:ngluon),dim=2)
             else
                k2 = czero 
             endif
             k2 = k2 + p(:,3)+kW
             k2sq = sc(k2,k2) !- mass2


             sp2 = spb2(sp2,k2)!+mass*sp2


             if (abs(k1sq) > propcut) then 
                tmp = -ci/k1sq*vgq(e1,sp2)
             else 
                tmp = czero 
             endif




             if (abs(k2sq) > propcut) then 
                tmp = ci/k2sq*tmp
             else 
                tmp = czero 
             endif

             !----------used to have czero here (now)
             res = res + tmp

          enddo  !#4

       endif

! XXXXX
       if ((sw.lt.4).and. ((case_b2 .eqv. .true.) .or. (WWqqqq .eqv. .false.))) then

          do m=0,ng3

             ngL = ng1+ ng2+m      

             e1 = gW_bff(e(:,1:ngL),k(:,1:ngL),sp(:,1),p(:,1),fl1,&
                  &sp(:,2),p(:,2),fl2,eW,kW,ng1,ng2,&
                  &giarray(1:ngL),qiarray(1:2),Wid,pol_int)    
             if (1<=ngL ) then 
                k1=sum(k(:,1:ngL),dim=2)+p(:,1)+p(:,2)+kW
             else
                k1=p(:,1)+p(:,2)+kW
             endif
             k1sq=sc(k1,k1)


             sp2=f(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
                  &sp(:,3),p(:,3),fl3,fl0,ng3-m,&
                  &giarray(ngL+1:ngluon),qiarray(3:3),pol_int)    
             if (ngL+1<=ngluon) then 
                k2=sum(k(:,ngL+1:ngluon),dim=2)
             else
                k2 = czero 
             endif
             k2 = k2 + p(:,3)
             k2sq = sc(k2,k2) !- mass2

             if (ng4 > 0.or. m < ng3) then 
                sp2 = spb2(sp2,k2)!+mass*sp2
             endif


             if (abs(k1sq) > propcut) then 
                tmp = -ci/k1sq*vgq(e1,sp2)
             else 
                tmp = czero 
             endif

             if (ng4 > 0.or.m < ng3) then 
                if (abs(k2sq) > propcut) then 
                   tmp = ci/k2sq*tmp
                else 
                   tmp = czero 
                endif
             endif


             !-------used to have czero here (now)

             if ((case_b2 .eqv. .false.) .and. (qbq_and_gluons .eqv. .true.))   tmp = czero
             if (extra_ferm_pair1) tmp = czero

             res = res + tmp

          enddo  !#5

       endif

! XXXXX
       if ((sw.gt.1) .and. (case_b2 .eqv. .false.)) then 

          do m=0,ng2

             sp1=f(e(:,1:ng1+m),k(:,1:ng1+m),sp(:,1),p(:,1),&
                  &fl1,fl0,ng1,&
                  &giarray(1:ng1+m),qiarray(1:1),pol_int)    
             if (1<=ng1+m) then 
                k1=sum(k(:,1:ng1+m),dim=2)
             else
                k1 = czero 
             endif 
             k1 = k1 + p(:,1)
             k1sq = sc(k1,k1) !- mass2



             if (ng1 > 0.or.m>0) then 
                sp1 = spb2(sp1,k1)!+mass*sp1
             endif



             e2 = gW_fbf(e(:,ng1+m+1:ngluon),k(:,ng1+m+1:ngluon),&
                  &sp(:,2),p(:,2),fl2,sp(:,3),p(:,3),fl3,eW,kW,ng2-m,ng3,&
                  &giarray(ng1+m+1:ngluon),qiarray(2:3),Wid,pol_int)    
             if (ng1+m+1<=ngluon) then 
                k2=sum(k(:,ng1+m+1:ngluon),dim=2)
             else
                k2 = czero 
             endif 
             k2 = k2 + p(:,2)+p(:,3)+kW
             k2sq=sc(k2,k2)


             if (abs(k2sq) > propcut) then 
                tmp = -ci/k2sq*vqg(sp1,e2)
             else 
                tmp = czero 
             endif

             if (ng1 > 0.or.m>0) then 
                if (abs(k1sq) > propcut) then 
                   tmp = ci/k1sq*tmp
                else 
                   tmp = czero 
                endif
             endif

             res = res + tmp



          enddo  !#6

       endif

       if ((sw.eq.1).or.(sw.eq.4)) then 

          if (fl0.eq.'top') flaux = 'bot'
          if (fl0.eq.'bot') flaux = 'top'
          if (fl0.eq.'chr') flaux = 'chr'

! XXXXXX
          if (WWqqqq) then 
             sp4 = f_bffbf_2(e,k,sp,p,fll,flaux,&
                  &ng1,ng2,ng3,&
                  &giarray,qiarray,pol_int)    
          else
             sp4 = f_bffbf(e,k,sp,p,fll,flaux,&
                  &ng1,ng2,ng3,&
                  &giarray,qiarray,pol_int)    
          endif


          if (1<=ngluon) then 
             k4 = sum(k(:,1:ngluon),dim=2)+p(:,1)+p(:,2)+p(:,3)
          else
             k4 = p(:,1)+p(:,2)+p(:,3)
          endif

          sp4 = spb2(sp4,k4) !+ mass*sp4
          k4sq = sc(k4,k4)

          tmp = vbqW(sp4,eW)

          if (abs(k4sq) > propcut) then 
             tmp = ci/k4sq*tmp
          else 
             tmp = czero
          endif

          res = res + tmp   ! #7

       endif

    endif

    ! -- store current 
!    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)


  end function fW_bffbf

!$$$    recursive function fW_bffbf_2W(e,k,sp,p,fll,fl0,&
!$$$         &eW,kW,ng1,ng2,ng3,sw,giarray,qiarray,Wid,pol_int) result(res)
!$$$      double complex, intent(in) :: e(:,:), k(:,:)
!$$$      double complex, intent(in) :: sp(:,:),p(:,:)
!$$$      double complex, intent(in) :: eW(:),kW(:)
!$$$      integer, intent(in) ::  ng1,ng2,ng3,sw
!$$$      character, intent(in) :: fll(:)*3
!$$$      character, intent(in) :: fl0*3   ! flavor off-shell f-line
!$$$      integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
!$$$      ! -----------------------------------------------------------------------
!$$$      character :: flaux*3
!$$$      character :: fl1*3,fl2*3,fl3*3
!$$$      integer :: ngluon, ng4, ngL,m!,m1
!$$$      integer, parameter :: Ndumm=0
!$$$      double complex             :: res(size(sp,dim=1))
!$$$      double complex             :: tmp(size(sp,dim=1))
!$$$      double complex             :: k1(size(k,dim=1))
!$$$      double complex             :: k2(size(k,dim=1))
!$$$      double complex             :: k4(size(k,dim=1))
!$$$      double complex             :: sp1(size(sp,dim=1))
!$$$      double complex             :: sp2(size(sp,dim=1))
!$$$      double complex             :: sp4(size(sp,dim=1))
!$$$      double complex             :: e1(size(e,dim=1))
!$$$      double complex             :: e2(size(e,dim=1))
!$$$  
!$$$      double complex             :: kdumm(size(k,dim=1),Ndumm)
!$$$      double complex             :: edumm(size(e,dim=1),Ndumm)
!$$$      double complex  :: k1sq,k2sq,k4sq
!$$$      !real(dp) :: mass,mass2
!$$$      logical                   :: done 
!$$$  
!$$$      !if (verbose) write(*,*) 'entering fW_bffbf_2W'
!$$$  
!$$$      done = .false. 
!$$$  
!$$$  !    if (present(giarray)) then 
!$$$  !       !if (size(qiarray) /= 3) stop 'fW_bffbf_2W: wrong size qiarray'
!$$$  !       !if (size(e,dim=2) /= size(giarray)) stop 'fW_bffbf: ng= size(giarray)' 
!$$$  !       call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!$$$  !       if (done) return 
!$$$  !    else
!$$$  !       if (i_warn < max_warn) then 
!$$$  !          write(*,*) 'fW_bffbf_2W: giarray missing', i_warn 
!$$$  !          i_warn = i_warn+1
!$$$  !       endif
!$$$  !    endif
!$$$  
!$$$      !mass = mt
!$$$      !mass2 = mass**2
!$$$  
!$$$      ngluon = size(e,dim=2)
!$$$      ng4 = ngluon - ng1 - ng2-ng3
!$$$  
!$$$      fl1 = fll(1)
!$$$      fl2 = fll(2)
!$$$      fl3 = fll(3)
!$$$  
!$$$      if ((sw.ne.3).and.(sw.ne.1)) stop 'sw in fW_bffbf_2W not implemented'
!$$$  
!$$$  
!$$$      if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C'
!$$$  
!$$$      if (ngluon == 0) then 
!$$$  
!$$$         res = czero
!$$$  
!$$$         if (sw.eq.3) then 
!$$$  
!$$$            e2 = gW_fbf(edumm,kdumm,sp(:,2),p(:,2),fl2,&
!$$$                 &sp(:,3),p(:,3),fl3,eW,kW,0,0,&
!$$$                 &giarray,qiarray(2:3),Wid,pol_int)    
!$$$            k1 = p(:,2)+p(:,3)+kW
!$$$            k1sq=sc(k1,k1)
!$$$  
!$$$            if (abs(k1sq) > propcut.and.fl0==fl1) then 
!$$$               tmp = -ci/k1sq*vqg(sp(:,1),e2)
!$$$            else 
!$$$               tmp = czero 
!$$$  	     write(*,*)'zerotmp'
!$$$            endif
!$$$  
!$$$            res = res + tmp  ! #1
!$$$         endif
!$$$  
!$$$         if ((sw.eq.1)) then 
!$$$  
!$$$            e2 = g_fbf(edumm,kdumm,sp(:,2),p(:,2),fl2,&
!$$$                 &sp(:,3),p(:,3),fl3,0,0,&
!$$$                 &giarray,qiarray(2:3),pol_int)    
!$$$  
!$$$            k2 = p(:,2)+p(:,3)
!$$$            k2sq = sc(k2,k2)
!$$$  
!$$$            sp4 = fW(edumm,kdumm,sp(:,1),p(:,1),fl1,fl0,eW,kW,0,&
!$$$                 &giarray,qiarray(1:1),Wid,pol_int)    
!$$$            k4  = p(:,1) + kW
!$$$            k4sq = sc(k4,k4)
!$$$            sp4 = spb2(sp4,k4) !+ mass*sp4
!$$$  
!$$$            tmp = vqg(sp4,e2)
!$$$  
!$$$            if (abs(k2sq) > propcut) then 
!$$$               tmp = -ci/k2sq*tmp
!$$$            else 
!$$$               tmp = czero 
!$$$  	     write(*,*)'czerotmp'
!$$$            endif
!$$$  
!$$$            if (abs(k4sq) > propcut) then 
!$$$               tmp = ci/k4sq*tmp 
!$$$            else 
!$$$               tmp = czero  
!$$$  	     write(*,*)'czerotmp'
!$$$            endif
!$$$  
!$$$            res = res + tmp  ! #3
!$$$  
!$$$            if (fl0.eq.'top') flaux = 'bot'
!$$$            if (fl0.eq.'bot') flaux = 'top'
!$$$  
!$$$            if (qbq_WW_and_gluons) then 
!$$$               sp4 = f_bffbf_2(edumm,kdumm,sp,p,fll,flaux,&
!$$$                    &ng1,ng2,ng3,&
!$$$                    &giarray,qiarray,pol_int)   
!$$$            else
!$$$               sp4 = f_bffbf(edumm,kdumm,sp,p,fll,flaux,&
!$$$                    &ng1,ng2,ng3,&
!$$$                    &giarray,qiarray,pol_int)    
!$$$  
!$$$            endif
!$$$  
!$$$            k4 = p(:,1)+p(:,2)+p(:,3)
!$$$            sp4 = spb2(sp4,k4) !+ mass*sp4
!$$$            k4sq = sc(k4,k4)
!$$$  
!$$$            tmp = vbqW(sp4,eW)
!$$$  
!$$$            if (abs(k4sq) > propcut) then 
!$$$               tmp = ci/k4sq*tmp
!$$$            else 
!$$$               tmp = czero 
!$$$  	     write(*,*)'czerotmp'
!$$$            endif
!$$$  
!$$$            res = res + tmp   ! #5
!$$$  
!$$$         endif
!$$$  
!$$$      else  ! -- this is for ngluon > 0
!$$$  
!$$$         res = czero
!$$$  
!$$$         if (sw.eq.1) then 
!$$$  
!$$$            do m=0,ng2
!$$$  
!$$$               sp1=fW(e(:,1:ng1+m),k(:,1:ng1+m),sp(:,1),p(:,1),&
!$$$                    &fl1,fl0,eW,kW,ng1,&
!$$$                    &giarray(1:ng1+m),qiarray(1:1),Wid,pol_int)    
!$$$               if (1<=ng1+m) then 
!$$$                  k1=sum(k(:,1:ng1+m),dim=2)
!$$$               else
!$$$                  k1 = czero 
!$$$               endif
!$$$               k1 = k1 + p(:,1)+kW
!$$$               k1sq = sc(k1,k1) !- mass2
!$$$  
!$$$               sp1 = spb2(sp1,k1)!+mass*sp1
!$$$  
!$$$               e2 = g_fbf(e(:,ng1+m+1:ngluon),k(:,ng1+m+1:ngluon),&
!$$$                    &sp(:,2),p(:,2),fl2,sp(:,3),p(:,3),fl3,ng2-m,ng3,&
!$$$                    &giarray(ng1+m+1:ngluon),qiarray(2:3),pol_int)    
!$$$               if (ng1+m+1<=ngluon) then 
!$$$                  k2=sum(k(:,ng1+m+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,2)+p(:,3)
!$$$               k2sq=sc(k2,k2)
!$$$  
!$$$  
!$$$               if (abs(k2sq) > propcut) then 
!$$$                  tmp = -ci/k2sq*vqg(sp1,e2)
!$$$               else 
!$$$                  tmp = czero  
!$$$  	     write(*,*)'czerotmp'
!$$$               endif
!$$$  
!$$$  
!$$$               if (abs(k1sq) > propcut) then 
!$$$                  tmp = ci/k1sq*tmp
!$$$               else 
!$$$                  tmp = czero  
!$$$  	     write(*,*)'czerotmp'
!$$$               endif
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$            enddo  !#1
!$$$  
!$$$         endif
!$$$  
!$$$  
!$$$         !--------- next step is valid for all sw
!$$$  
!$$$         do m=0,ng4-1
!$$$  
!$$$            ngL = ng1+ ng2+ng3+m      
!$$$  
!$$$            sp1=fW_bffbf_2W(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0,&
!$$$                 &eW,kW,ng1,ng2,ng3,sw,&
!$$$                 &giarray(1:ngL),qiarray,Wid,pol_int)    
!$$$            if (1<=ngL) then 
!$$$               k1=sum(k(:,1:ngL),dim=2)
!$$$            else
!$$$               k1 = czero 
!$$$            endif
!$$$            k1 = k1 + p(:,1)+p(:,2)+p(:,3)+kW
!$$$            k1sq = sc(k1,k1) !- mass2
!$$$  
!$$$            sp1 = spb2(sp1,k1) !+ mass*sp1
!$$$  
!$$$            e2 = vgluon(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
!$$$                 &giarray(ngL+1:ngluon),pol_int)    
!$$$            if (ngL+1<=ngluon) then 
!$$$               k2=sum(k(:,ngL+1:ngluon),dim=2)
!$$$            else
!$$$               k2 = czero 
!$$$            endif
!$$$            k2sq=sc(k2,k2)
!$$$  
!$$$            if (abs(k1sq) > propcut) then 
!$$$               tmp = ci/k1sq*vqg(sp1,e2)
!$$$            else 
!$$$               tmp = czero  
!$$$  	     write(*,*)'czerotmp'
!$$$            endif
!$$$  
!$$$            if (m < ng4-1) then 
!$$$               if (abs(k2sq) > propcut) then 
!$$$                  tmp = -ci/k2sq*tmp
!$$$               else 
!$$$                  tmp = czero  
!$$$  	     write(*,*)'czerotmp'
!$$$               endif
!$$$            endif
!$$$  
!$$$            res = res + tmp
!$$$  
!$$$  
!$$$         enddo  !#2
!$$$  
!$$$  
!$$$         !-valid for all sw
!$$$  
!$$$         do m=1,ng1
!$$$  
!$$$            e1 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)    
!$$$            k1=sum(k(:,1:m),dim=2)
!$$$            k1sq=sc(k1,k1)
!$$$  
!$$$            sp2=fW_bffbf_2W(e(:,m+1:ngluon),k(:,m+1:ngluon),&
!$$$                 &sp,p,fll,fl0,eW,kW,ng1-m,ng2,ng3,sw,&
!$$$                 &giarray(m+1:ngluon),qiarray,Wid,pol_int)    
!$$$            if (m+1<=ngluon) then 
!$$$               k2=sum(k(:,m+1:ngluon),dim=2)
!$$$            else
!$$$               k2 = czero 
!$$$            endif
!$$$            k2 = k2 + p(:,1)+p(:,2)+p(:,3)+kW
!$$$            k2sq = sc(k2,k2) !- mass2
!$$$            sp2 = spb2(sp2,k2)!+mass*sp2
!$$$  
!$$$  
!$$$            if (abs(k2sq) > propcut) then 
!$$$               tmp = ci/k2sq*vgq(e1,sp2)
!$$$            else 
!$$$               tmp = czero  
!$$$  	     write(*,*)'czerotmp'
!$$$            endif
!$$$  
!$$$  
!$$$            if (m > 1) then 
!$$$               if (abs(k1sq) > propcut) then 
!$$$                  tmp = -ci/k1sq*tmp
!$$$               else 
!$$$                  tmp = czero  
!$$$  	     write(*,*)'czerotmp'
!$$$               endif
!$$$            endif
!$$$  
!$$$            res = res + tmp
!$$$  
!$$$  
!$$$         enddo  !#3
!$$$  
!$$$  
!$$$         if (sw.eq.3) then 
!$$$  
!$$$            do m=0,ng2
!$$$  
!$$$               sp1=f(e(:,1:ng1+m),k(:,1:ng1+m),sp(:,1),p(:,1),&
!$$$                    &fl1,fl0,ng1,&
!$$$                    &giarray(1:ng1+m),qiarray(1:1),pol_int)    
!$$$               if (1<=ng1+m) then 
!$$$                  k1=sum(k(:,1:ng1+m),dim=2)
!$$$               else
!$$$                  k1 = czero 
!$$$               endif 
!$$$               k1 = k1 + p(:,1)
!$$$               k1sq = sc(k1,k1) !- mass2
!$$$  
!$$$  
!$$$  
!$$$               if (ng1 > 0.or.m>0) then 
!$$$                  sp1 = spb2(sp1,k1)!+mass*sp1
!$$$               endif
!$$$  
!$$$  
!$$$  
!$$$               e2 = gW_fbf(e(:,ng1+m+1:ngluon),k(:,ng1+m+1:ngluon),&
!$$$                    &sp(:,2),p(:,2),fl2,sp(:,3),p(:,3),fl3,eW,kW,ng2-m,ng3,&
!$$$                    &giarray(ng1+m+1:ngluon),qiarray(2:3),Wid,pol_int)    
!$$$               if (ng1+m+1<=ngluon) then 
!$$$                  k2=sum(k(:,ng1+m+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif 
!$$$               k2 = k2 + p(:,2)+p(:,3)+kW
!$$$               k2sq=sc(k2,k2)
!$$$  
!$$$  
!$$$               if (abs(k2sq) > propcut) then 
!$$$                  tmp = -ci/k2sq*vqg(sp1,e2)
!$$$               else 
!$$$                  tmp = czero  
!$$$  	     write(*,*)'czerotmp'
!$$$               endif
!$$$  
!$$$               if (ng1 > 0.or.m>0) then 
!$$$                  if (abs(k1sq) > propcut) then 
!$$$                     tmp = ci/k1sq*tmp
!$$$                  else 
!$$$                     tmp = czero  
!$$$  	     write(*,*)'czerotmp'
!$$$                  endif
!$$$               endif
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$  
!$$$  
!$$$            enddo  !#6
!$$$  
!$$$         endif
!$$$  
!$$$         if ((sw.eq.1)) then 
!$$$  
!$$$            if (fl0.eq.'top') flaux = 'bot'
!$$$            if (fl0.eq.'bot') flaux = 'top'
!$$$  
!$$$            sp4 = f_bffbf_2(e,k,sp,p,fll,flaux,&
!$$$                 &ng1,ng2,ng3,&
!$$$                 &giarray,qiarray,pol_int)    
!$$$  
!$$$  
!$$$            if (1<=ngluon) then 
!$$$               k4 = sum(k(:,1:ngluon),dim=2)+p(:,1)+p(:,2)+p(:,3)
!$$$            else
!$$$               k4 = p(:,1)+p(:,2)+p(:,3)
!$$$            endif
!$$$  
!$$$            sp4 = spb2(sp4,k4) !+ mass*sp4
!$$$            k4sq = sc(k4,k4)
!$$$  
!$$$            tmp = vbqW(sp4,eW)
!$$$  
!$$$            if (abs(k4sq) > propcut) then 
!$$$               tmp = ci/k4sq*tmp
!$$$            else 
!$$$               tmp = czero 
!$$$  	     write(*,*)'czerotmp'
!$$$            endif
!$$$  
!$$$            res = res + tmp   ! #7
!$$$  
!$$$         endif
!$$$  
!$$$      endif
!$$$  
!$$$      ! -- store current 
!$$$      if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)
!$$$  
!$$$  
!$$$    end function fW_bffbf_2W
!$$$  
!$$$  
!$$$  
!$$$  
!$$$  
!$$$    recursive function gW_sbsfbf(e,k,sp,p,fll,eW,kW,&
!$$$         &ng1,ng2,ng3,ng4,sw,giarray,qiarray,Wid,pol_int) result(res)
!$$$      double complex, intent(in) :: e(:,:), k(:,:)
!$$$      double complex, intent(in) :: sp(:,:), p(:,:)
!$$$      integer, intent(in) ::  ng1,ng2, ng3,ng4
!$$$      character, intent(in) :: fll(:)*3
!$$$      integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
!$$$      ! -----------------------------------------------------------------------
!$$$      double complex, intent(in) :: eW(:),kW(:)
!$$$      integer, intent(in) :: sw
!$$$      character :: fl1*3,fl2*3,fl3*3,fl4*3
!$$$  !    integer             :: m1,m2,m3, ms1a, ms2a,Dv!,m
!$$$      integer :: ngluon, ng5,Dv!,  ngL
!$$$      integer, parameter :: Ndumm = 0
!$$$      double complex             :: res(size(e,dim=1))
!$$$      double complex             :: tmp(size(e,dim=1))
!$$$      double complex             :: k1(size(p,dim=1))
!$$$      double complex             :: k2(size(p,dim=1))
!$$$      double complex             :: sp1(size(sp,dim=1))
!$$$      double complex             :: sp2(size(sp,dim=1))
!$$$      double complex             :: e1(size(e,dim=1))
!$$$      double complex             :: e2(size(e,dim=1))
!$$$      double complex             :: kdumm(size(k,dim=1),Ndumm)
!$$$      double complex             :: edumm(size(e,dim=1),Ndumm)
!$$$      double complex  :: k1sq,k2sq!,k3sq!,k4sq
!$$$      !real(dp) :: mass
!$$$      logical                   :: done 
!$$$  
!$$$      !if (verbose) write(*,*) 'entering gW_sbsfbf'
!$$$  
!$$$      done = .false. 
!$$$  !    if (present(giarray)) then 
!$$$  !       !if (size(qiarray) /= 4) stop 'gW_sbsfbf: wrong size qiarray'
!$$$  !       !if (size(e,dim=2) /= size(giarray)) stop 'gW_sbsfbf: ng= size(giarray)' 
!$$$  !       call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!$$$  !       if (done) return 
!$$$  !    else
!$$$  !       if (i_warn < max_warn) then 
!$$$  !          write(*,*) 'gW_sbsfbf: giarray missing', i_warn 
!$$$  !          i_warn = i_warn+1
!$$$  !       endif
!$$$  !    endif
!$$$  
!$$$      !mass = mt
!$$$      Dv = size(e,dim=1)
!$$$      ngluon = size(e,dim=2)
!$$$  
!$$$      if ((sw.ne.2).and.(sw.ne.4)) print *, 'W position incorrect'
!$$$  
!$$$      ng5 = ngluon - ng1 - ng2-ng3 - ng4
!$$$  
!$$$  
!$$$      fl1 = fll(1)
!$$$      fl2 = fll(2)
!$$$      fl3 = fll(3)
!$$$      fl4 = fll(4)
!$$$  
!$$$  
!$$$  
!$$$     if ((fl1.ne.fl2).and.(fl3.ne.fl4) .and. (WWqqqq .eqv. .false.)) then 
!$$$  !        if ((fl1.ne.fl2).and.(fl3.ne.fl4)) then 
!$$$         print *, 'flavors not consistent, gw_sbsfbf'
!$$$  
!$$$         stop 
!$$$  
!$$$      endif
!$$$  
!$$$      if (ng5 < 0) write(*,*) 'ERROR IN CURRENT G'
!$$$  
!$$$  
!$$$      if (ngluon == 0) then 
!$$$  
!$$$         res = czero
!$$$  
!$$$  
!$$$  
!$$$         if (sw.eq.2) then 
!$$$  
!$$$            e1 = gW_fbf(edumm,kdumm,sp(:,1),p(:,1),fl1,&
!$$$                 &sp(:,2),p(:,2),fl2,eW,kW,0,0,&
!$$$                 &giarray,qiarray(1:2),Wid,pol_int)
!$$$            k1 = p(:,1) + p(:,2)+kW 
!$$$            k1sq = sc(k1,k1)
!$$$  
!$$$            e2 = g_fbf(edumm,kdumm,sp(:,3),p(:,3),fl3,&
!$$$                 &sp(:,4),p(:,4),fl4,0,0,&
!$$$                 &giarray,qiarray(3:4),pol_int)
!$$$            k2 = p(:,3) + p(:,4)
!$$$            k2sq = sc(k2,k2)
!$$$  
!$$$            tmp = (-ci/k1sq)*(-ci/k2sq)*vggg(e1,k1,e2,k2)
!$$$  
!$$$         endif
!$$$  
!$$$  
!$$$         if (sw.eq.4) then 
!$$$  
!$$$            e1 = g_fbf(edumm,kdumm,sp(:,1),p(:,1),fl1,&
!$$$                 &sp(:,2),p(:,2),fl2,0,0,giarray,qiarray(1:2),pol_int)
!$$$            k1 = p(:,1) + p(:,2) 
!$$$            k1sq = sc(k1,k1)
!$$$  
!$$$            e2 = gW_fbf(edumm,kdumm,sp(:,3),p(:,3),fl3,&
!$$$                 &sp(:,4),p(:,4),fl4,eW,kW,0,0,&
!$$$                 &giarray,qiarray(3:4),Wid,pol_int)
!$$$            k2 = p(:,3) + p(:,4)+kW
!$$$            k2sq = sc(k2,k2)
!$$$  
!$$$            tmp = (-ci/k1sq)*(-ci/k2sq)*vggg(e1,k1,e2,k2)
!$$$  
!$$$         endif
!$$$  
!$$$  
!$$$         res = res + tmp  !# 1
!$$$  
!$$$  
!$$$         if ((sw.eq.2).or.(sw.eq.4)) then 
!$$$  
!$$$            sp2 = fW_bffbf(edumm,kdumm,sp(:,2:4),p(:,2:4),&
!$$$                 &fll(2:4),fl1,eW,kW,0,0,0,sw-1,&
!$$$                 &giarray,qiarray(2:4),Wid,pol_int)
!$$$            k2 = p(:,2) + p(:,3) + p(:,4)+kW
!$$$            k2sq = sc(k2,k2)
!$$$            sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$            sp2 = ci/k2sq*sp2
!$$$  
!$$$            sp1 = sp(:,1)
!$$$  
!$$$            tmp = -cone*vbqq(Dv,sp2,sp1)
!$$$  
!$$$            res = res + tmp       ! # 2
!$$$  
!$$$  
!$$$  
!$$$         endif
!$$$  
!$$$  
!$$$         if ((sw.eq.2)) then 
!$$$  
!$$$  ! XXXX
!$$$            if (WWqqqq .eqv. .true.) then
!$$$            
!$$$               sp2 = f_bffbf_2(edumm,kdumm,sp(:,2:4),p(:,2:4),&
!$$$                    &fll(2:4),fl2,0,0,0,&
!$$$                    &giarray,qiarray(2:4),pol_int)
!$$$            else
!$$$               sp2 = f_bffbf(edumm,kdumm,sp(:,2:4),p(:,2:4),&
!$$$                    &fll(2:4),fl2,0,0,0,&
!$$$                    &giarray,qiarray(2:4),pol_int)
!$$$            endif
!$$$  
!$$$            k2 = p(:,2) + p(:,3) + p(:,4) 
!$$$            k2sq = sc(k2,k2)
!$$$  
!$$$            sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$            sp2 = ci/k2sq*sp2
!$$$  
!$$$  ! XXXX        
!$$$            if (case_a4) then
!$$$       
!$$$               sp1 = bfW(edumm,kdumm,sp(:,1),p(:,1),fl1,fl4,&
!$$$                    &eW,kW,0,giarray,qiarray(1:1),Wid,pol_int)
!$$$            else
!$$$               sp1 = bfW(edumm,kdumm,sp(:,1),p(:,1),fl1,fl2,&
!$$$                    &eW,kW,0,giarray,qiarray(1:1),Wid,pol_int)
!$$$            endif
!$$$         
!$$$            k1 = -p(:,1)-kW
!$$$            sp1 = spi2(k1,sp1) !+ mass*sp1
!$$$            k1sq = sc(k1,k1)
!$$$  
!$$$            sp1 = ci/k1sq*sp1
!$$$  
!$$$            tmp = -cone*vbqq(Dv,sp2,sp1)
!$$$  
!$$$            res = res + tmp       ! # 2
!$$$  
!$$$         endif
!$$$  
!$$$  
!$$$         if ((sw.eq.2).or.(sw.eq.4)) then 
!$$$  
!$$$  
!$$$  
!$$$            sp2 = bfW_fbff(edumm,kdumm,sp(:,1:3),p(:,1:3),&
!$$$                 &fll(1:3),fl4,eW,kW,0,0,0,sw,&
!$$$                 &giarray,qiarray(1:3),Wid,pol_int)
!$$$  
!$$$            k2 = -p(:,1) - p(:,2) - p(:,3)-kW
!$$$            k2sq = sc(k2,k2)
!$$$  
!$$$            sp2 = spi2(k2,sp2) !+ mass*sp2
!$$$            sp2 = ci/k2sq*sp2
!$$$  
!$$$            sp1 = sp(:,4)
!$$$  
!$$$            tmp = -cone*vbqq(Dv,sp1,sp2)
!$$$  
!$$$            res = res + tmp       ! # 3
!$$$  
!$$$         endif
!$$$  
!$$$         if (sw.eq.4) then 
!$$$  
!$$$            sp2 = bf_fbff(edumm,kdumm,sp(:,1:3),p(:,1:3),&
!$$$                 &fll(1:3),fl3,0,0,0,&
!$$$                 &giarray,qiarray(1:3),pol_int)
!$$$  
!$$$            k2 = -p(:,1) - p(:,2) - p(:,3)
!$$$            k2sq = sc(k2,k2)
!$$$  
!$$$            sp2 = spi2(k2,sp2) !+ mass*sp2
!$$$            sp2 = ci/k2sq*sp2
!$$$  
!$$$  
!$$$            sp1 = fW(edumm,kdumm,sp(:,4),p(:,4),fl4,fl3,&
!$$$                 &eW,kW,0,giarray,qiarray(4:4),Wid,pol_int)
!$$$            k1 = kW+p(:,4)
!$$$            k1sq=sc(k1,k1)
!$$$            sp1 = spb2(sp1,k1) !+ mass*sp1
!$$$            sp1 = ci/k1sq*sp1
!$$$  
!$$$            tmp = -cone*vbqq(Dv,sp1,sp2)
!$$$  
!$$$            res = res + tmp       ! # 3
!$$$  
!$$$         endif
!$$$  
!$$$  
!$$$      else  ! for ngluon > 0
!$$$  
!$$$  
!$$$         print *, 'gW_sbsfbf ngl > 0 not done'
!$$$  
!$$$      endif   !end condition for ngluon 
!$$$  
!$$$      ! -- store current 
!$$$      if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)
!$$$  
!$$$  
!$$$    end function gW_sbsfbf
!$$$  
!$$$  
!$$$  
!$$$    !-----this is a new function, for four-fermion amplitudes 
!$$$    !---- note, that no sw here
!$$$    recursive function gW_sbsfbf_1(e,k,sp,p,fll,eW,kW,&
!$$$         &ng1,ng2,ng3,ng4,giarray,qiarray,Wid,pol_int) result(res)
!$$$      double complex, intent(in) :: e(:,:), k(:,:)
!$$$      double complex, intent(in) :: sp(:,:), p(:,:)
!$$$      integer, intent(in) ::  ng1,ng2, ng3,ng4
!$$$      character, intent(in) :: fll(:)*3
!$$$      integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
!$$$      ! -----------------------------------------------------------------------
!$$$      double complex, intent(in) :: eW(:),kW(:)
!$$$      character :: fl1*3,fl2*3,fl3*3,fl4*3
!$$$  
!$$$      integer :: ngluon, ng5,Dv
!$$$      integer, parameter :: Ndumm = 0
!$$$      double complex             :: res(size(e,dim=1))
!$$$      double complex             :: tmp(size(e,dim=1))
!$$$      double complex             :: k1(size(p,dim=1))
!$$$      double complex             :: k2(size(p,dim=1))
!$$$      double complex             :: sp1(size(sp,dim=1))
!$$$      double complex             :: sp2(size(sp,dim=1))
!$$$      double complex             :: kdumm(size(k,dim=1),Ndumm)
!$$$      double complex             :: edumm(size(e,dim=1),Ndumm)
!$$$      double complex  :: k1sq,k2sq
!$$$      !real(dp) :: mass
!$$$      logical                   :: done 
!$$$  
!$$$      !if (verbose) write(*,*) 'entering gW_sbsfbf_1'
!$$$  
!$$$      done = .false. 
!$$$  !    if (present(giarray)) then 
!$$$  !       !if (size(qiarray) /= 4) stop 'gW_sbsfbf_1: wrong size qiarray'
!$$$  !       !if (size(e,dim=2) /= size(giarray)) stop 'gW_sbsfbf_1: ng= size(giarray)' 
!$$$  !       call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!$$$  !       if (done) return 
!$$$  !    else
!$$$  !       if (i_warn < max_warn) then 
!$$$  !          write(*,*) 'gW_sbsfbf_1: giarray missing', i_warn 
!$$$  !          i_warn = i_warn+1
!$$$  !       endif
!$$$  !    endif
!$$$  
!$$$      !mass = mt
!$$$  
!$$$      Dv = size(e,dim=1)
!$$$  
!$$$      ngluon = size(e,dim=2)
!$$$  
!$$$      ng5 = ngluon - ng1 - ng2-ng3 - ng4
!$$$  
!$$$  
!$$$      fl1 = fll(1)
!$$$      fl2 = fll(2)
!$$$      fl3 = fll(3)
!$$$      fl4 = fll(4)
!$$$  
!$$$      if ((fl2.ne.'str').or.(fl3.ne.'str')) then 
!$$$  
!$$$         print *, 'flavors not consistent'
!$$$  
!$$$         stop 
!$$$  
!$$$      endif
!$$$  
!$$$      if (ng5 < 0) write(*,*) 'ERROR IN CURRENT G'
!$$$  
!$$$  
!$$$      if (ngluon == 0) then 
!$$$  
!$$$         res = czero
!$$$  
!$$$  
!$$$         sp2 = fW_bffbf(edumm,kdumm,sp(:,2:4),p(:,2:4),&
!$$$              &fll(2:4),fl1,eW,kW,0,0,0,4,&
!$$$              &giarray,qiarray(2:4),Wid,pol_int)
!$$$         k2 = p(:,2) + p(:,3) + p(:,4)+kW
!$$$         k2sq = sc(k2,k2)
!$$$         sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$         sp2 = ci/k2sq*sp2
!$$$  
!$$$         sp1 = sp(:,1)
!$$$  
!$$$         tmp = -cone*vbqq(Dv,sp2,sp1)
!$$$  
!$$$  
!$$$  
!$$$         res = res + tmp       ! # 1
!$$$  
!$$$         tmp = czero
!$$$  
!$$$  
!$$$         sp2 = f_bffbf(edumm,kdumm,sp(:,2:4),p(:,2:4),&
!$$$              &fll(2:4),fl4,0,0,0,&
!$$$              &giarray,qiarray(2:4),pol_int)
!$$$         k2 = p(:,2) + p(:,3) + p(:,4) 
!$$$         k2sq = sc(k2,k2)
!$$$  
!$$$         sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$         sp2 = ci/k2sq*sp2
!$$$  
!$$$  
!$$$         sp1 = bfW(edumm,kdumm,sp(:,1),p(:,1),fl1,fl4,&
!$$$              &eW,kW,0,giarray,qiarray(1:1),Wid,pol_int)
!$$$         k1 = -p(:,1)-kW
!$$$         sp1 = spi2(k1,sp1) !+ mass*sp1
!$$$         k1sq = sc(k1,k1)
!$$$  
!$$$         sp1 = ci/k1sq*sp1
!$$$  
!$$$         tmp = -cone*vbqq(Dv,sp2,sp1)
!$$$  
!$$$  
!$$$         res = res + tmp       ! # 2
!$$$  
!$$$  
!$$$  
!$$$  
!$$$         sp2 = bfW_fbff(edumm,kdumm,sp(:,1:3),p(:,1:3),&
!$$$              &fll(1:3),fl4,eW,kW,0,0,0,1,&
!$$$              &giarray,qiarray(1:3),Wid,pol_int)
!$$$  
!$$$  
!$$$         k2 = -p(:,1) - p(:,2) - p(:,3)-kW
!$$$         k2sq = sc(k2,k2)
!$$$  
!$$$         sp2 = spi2(k2,sp2) !+ mass*sp2
!$$$         sp2 = ci/k2sq*sp2
!$$$  
!$$$         sp1 = sp(:,4)
!$$$  
!$$$         tmp = -cone*vbqq(Dv,sp1,sp2)
!$$$  
!$$$         res = res + tmp       ! # 3
!$$$  
!$$$  
!$$$  
!$$$         sp2 = bf_fbff(edumm,kdumm,sp(:,1:3),p(:,1:3),&
!$$$              &fll(1:3),fl1,0,0,0,&
!$$$              &giarray,qiarray(1:3),pol_int)
!$$$  
!$$$         k2 = -p(:,1) - p(:,2) - p(:,3)
!$$$         k2sq = sc(k2,k2)
!$$$  
!$$$         sp2 = spi2(k2,sp2) !+ mass*sp2
!$$$         sp2 = ci/k2sq*sp2
!$$$  
!$$$  
!$$$         sp1 = fW(edumm,kdumm,sp(:,4),p(:,4),fl4,fl1,&
!$$$              &eW,kW,0,&
!$$$              &giarray,qiarray(4:4),Wid,pol_int)
!$$$         k1 = kW+p(:,4)
!$$$         k1sq=sc(k1,k1)
!$$$         sp1 = spb2(sp1,k1) !+ mass*sp1
!$$$         sp1 = ci/k1sq*sp1
!$$$  
!$$$         tmp = -cone*vbqq(Dv,sp1,sp2)
!$$$  
!$$$         res = res + tmp       ! # 4
!$$$  
!$$$  
!$$$  
!$$$  
!$$$      else  ! for ngluon > 0
!$$$  
!$$$  
!$$$         print *, 'gW_sbsfbf_1 ngl > 0 not done'
!$$$  
!$$$      endif   !end condition for ngluon 
!$$$  
!$$$      ! -- store current 
!$$$      if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)
!$$$  
!$$$  
!$$$    end function gW_sbsfbf_1
!$$$  
!$$$  
!$$$  
!$$$  
!$$$    recursive function gW_sbsfbf_2(e,k,sp,p,fll,eW,kW,&
!$$$         &ng1,ng2,ng3,ng4,giarray,qiarray,Wid,pol_int) result(res)
!$$$      double complex, intent(in) :: e(:,:), k(:,:)
!$$$      double complex, intent(in) :: sp(:,:), p(:,:)
!$$$      integer, intent(in) ::  ng1,ng2, ng3,ng4
!$$$      character, intent(in) :: fll(:)*3
!$$$      double complex, intent(in) :: eW(:),kW(:)
!$$$      integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
!$$$      ! -----------------------------------------------------------------------
!$$$      character :: fl1*3,fl2*3,fl3*3,fl4*3
!$$$      integer :: ngluon, ng5,Dv
!$$$      integer, parameter :: Ndumm = 0
!$$$      double complex             :: res(size(e,dim=1))
!$$$      double complex             :: tmp(size(e,dim=1))
!$$$      double complex             :: k2(size(p,dim=1))
!$$$      double complex             :: sp1(size(sp,dim=1))
!$$$      double complex             :: sp2(size(sp,dim=1))
!$$$      double complex             :: kdumm(size(k,dim=1),Ndumm)
!$$$      double complex             :: edumm(size(e,dim=1),Ndumm)
!$$$      double complex             :: k2sq
!$$$      logical                   :: done 
!$$$      !real(dp) :: mass 
!$$$  
!$$$      !if (verbose) write(*,*) 'entering gW_sbsfbf_1'
!$$$  
!$$$      done = .false. 
!$$$  !    if (present(giarray)) then 
!$$$  !       !if (size(qiarray) /= 4) stop 'gW_sbsfbf_1: wrong size qiarray'
!$$$  !       !if (size(e,dim=2) /= size(giarray)) stop 'gW_sbsfbf_1: ng= size(giarray)' 
!$$$  !       call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!$$$  !       if (done) return 
!$$$  !    else
!$$$  !       if (i_warn < max_warn) then 
!$$$  !          write(*,*) 'gW_sbsfbf_1: giarray missing', i_warn 
!$$$  !          i_warn = i_warn+1
!$$$  !       endif
!$$$  !    endif
!$$$  
!$$$      !mass = mt
!$$$  
!$$$      Dv = size(e,dim=1)
!$$$  
!$$$      ngluon = size(e,dim=2)
!$$$  
!$$$  
!$$$      ng5 = ngluon - ng1 - ng2-ng3 - ng4
!$$$  
!$$$  
!$$$      fl1 = fll(1)
!$$$      fl2 = fll(2)
!$$$      fl3 = fll(3)
!$$$      fl4 = fll(4)
!$$$  
!$$$      if ((fl1.ne.'str').or.(fl4.ne.'str')) then 
!$$$  
!$$$         print *, 'flavors not consistent'
!$$$  
!$$$         stop 
!$$$  
!$$$      endif
!$$$  
!$$$      if (ng5 < 0) write(*,*) 'ERROR IN CURRENT G'
!$$$  
!$$$  
!$$$      if (ngluon == 0) then 
!$$$  
!$$$         res = czero
!$$$  
!$$$  
!$$$         sp2 = fW_bffbf(edumm,kdumm,sp(:,2:4),p(:,2:4),&
!$$$              &fll(2:4),fl1,eW,kW,0,0,0,2,&
!$$$              &giarray,qiarray(2:4),Wid,pol_int)
!$$$         k2 = p(:,2) + p(:,3) + p(:,4)+kW
!$$$         k2sq = sc(k2,k2)
!$$$         sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$         sp2 = ci/k2sq*sp2
!$$$  
!$$$         sp1 = sp(:,1)
!$$$  
!$$$         tmp = -cone*vbqq(Dv,sp2,sp1)
!$$$  
!$$$  
!$$$         res = res + tmp       ! # 1
!$$$  
!$$$         sp2 = bfW_fbff(edumm,kdumm,sp(:,1:3),p(:,1:3),&
!$$$              &fll(1:3),fl4,eW,kW,0,0,0,3,&
!$$$              &giarray,qiarray(1:3),Wid,pol_int)
!$$$  
!$$$  
!$$$         k2 = -p(:,1) - p(:,2) - p(:,3)-kW
!$$$         k2sq = sc(k2,k2)
!$$$  
!$$$         sp2 = spi2(k2,sp2) !+ mass*sp2
!$$$         sp2 = ci/k2sq*sp2
!$$$  
!$$$         sp1 = sp(:,4)
!$$$  
!$$$         tmp = -cone*vbqq(Dv,sp1,sp2)
!$$$  
!$$$         res = res + tmp       ! # 3
!$$$  
!$$$  
!$$$      else  ! for ngluon > 0
!$$$  
!$$$  
!$$$         print *, 'gW_sbsfbf_2 ngl > 0 not done'
!$$$  
!$$$      endif   !end condition for ngluon 
!$$$  
!$$$      ! -- store current 
!$$$      if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)
!$$$  
!$$$    end function gW_sbsfbf_2
!$$$  
!$$$  
!$$$  
!$$$  !----------------------------------
!$$$  recursive function gW_sbsfbf_3(e,k,sp,p,fll,eW,kW,&
!$$$         &ng1,ng2,ng3,ng4,sw,giarray,qiarray,Wid,pol_int) result(res)
!$$$      double complex, intent(in) :: e(:,:), k(:,:)
!$$$      double complex, intent(in) :: sp(:,:), p(:,:)
!$$$      integer, intent(in) ::  ng1,ng2, ng3,ng4,sw
!$$$      character, intent(in) :: fll(:)*3
!$$$      double complex, intent(in) :: eW(:),kW(:)
!$$$      integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
!$$$      ! -----------------------------------------------------------------------
!$$$      character :: fl1*3,fl2*3,fl3*3,fl4*3
!$$$      integer :: ngluon, ng5,Dv
!$$$      integer, parameter :: Ndumm = 0
!$$$      double complex             :: res(size(e,dim=1))
!$$$      double complex             :: tmp(size(e,dim=1))
!$$$      double complex             :: k1(size(p,dim=1))
!$$$      double complex             :: k2(size(p,dim=1))
!$$$      double complex             :: k3(size(p,dim=1))
!$$$      double complex             :: k4(size(p,dim=1))
!$$$      double complex             :: sp1(size(sp,dim=1))
!$$$      double complex             :: sp2(size(sp,dim=1))
!$$$      double complex             :: sp3(size(sp,dim=1))
!$$$      double complex             :: sp4(size(sp,dim=1))
!$$$      double complex             :: kdumm(size(k,dim=1),Ndumm)
!$$$      double complex             :: edumm(size(e,dim=1),Ndumm)
!$$$      double complex  :: k2sq,k3sq,k1sq,k4sq
!$$$      logical                   :: done 
!$$$      !real(dp) :: mass 
!$$$  
!$$$      !if (verbose) write(*,*) 'entering gW_sbsfbf_1'
!$$$  
!$$$      done = .false. 
!$$$      if (present(giarray)) then 
!$$$         !if (size(qiarray) /= 4) stop 'gW_sbsfbf_1: wrong size qiarray'
!$$$         !if (size(e,dim=2) /= size(giarray)) stop 'gW_sbsfbf_1: ng= size(giarray)' 
!$$$         ! XXX 
!$$$         call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!$$$         if (done) return
!$$$      else
!$$$         write(*,*) 'giarray missing'
!$$$      endif
!$$$  
!$$$      !mass = mt
!$$$  
!$$$      Dv = size(e,dim=1)
!$$$  
!$$$      ngluon = size(e,dim=2)
!$$$  
!$$$  
!$$$      ng5 = ngluon - ng1 - ng2-ng3 - ng4
!$$$  
!$$$  
!$$$      fl1 = fll(1)
!$$$      fl2 = fll(2)
!$$$      fl3 = fll(3)
!$$$      fl4 = fll(4)
!$$$  
!$$$  
!$$$      if ((fl2.ne.'str').or.(fl3.ne.'str')) then 
!$$$         res = czero
!$$$         write(*,*) 'Warning: gW_sbsfbf_3 = 0 as fl2 /= str or fl3 /= str'
!$$$         stop
!$$$         return
!$$$  
!$$$      endif
!$$$  
!$$$      if (ng5 < 0) write(*,*) 'ERROR IN CURRENT G'
!$$$  
!$$$  
!$$$      if (ngluon == 0) then 
!$$$  
!$$$         res = czero
!$$$         
!$$$         if (sw == 2) then
!$$$            
!$$$            sp1 = fW_bffbf_2(e,k,sp(:,2:4),p(:,2:4),fll(2:4),fl1,eW,kW,&
!$$$                 &ng1,ng2,ng3,1,giarray,qiarray(2:4),Wid,pol_int)
!$$$         
!$$$            k1 = p(:,2)+p(:,3)+p(:,4)+kW
!$$$            k1sq = sc(k1,k1)
!$$$  
!$$$            sp1 = spb2(sp1,k1)!+mass*sp1
!$$$       
!$$$            sp2 = sp(:,1)
!$$$            tmp = -cone*vbqq(Dv,sp1,sp2)
!$$$           
!$$$        
!$$$            if (abs(k1sq) > propcut) then
!$$$               tmp = ci*tmp/k1sq
!$$$            else
!$$$               tmp = czero
!$$$            endif
!$$$            
!$$$            res = res + tmp
!$$$          
!$$$            sp2 = f_bffbf_3(e,k,sp(:,2:4),p(:,2:4),fll(2:4),fl4,ng1,ng2,ng3,&
!$$$                 &giarray,qiarray(2:4),pol_int)
!$$$            k2 = p(:,2)+p(:,3)+p(:,4)
!$$$            k2sq = sc(k2,k2)
!$$$            
!$$$            sp2 = spb2(sp2,k2)!+mass*sp2
!$$$  
!$$$            sp3 = bfW(edumm,kdumm,sp(:,1),p(:,1),fl1,fl4,eW,kW,ng1,giarray,&
!$$$                 &qiarray(1:1),pol_int)
!$$$  
!$$$            k3 = -p(:,1)-kW
!$$$            k3sq = sc(k3,k3)
!$$$  
!$$$            sp3 = spi2(k3,sp3)!+mass*sp3
!$$$  
!$$$            tmp = -cone*vbqq(Dv,sp2,sp3)
!$$$  
!$$$            
!$$$             if (abs(k2sq) > propcut) then
!$$$               tmp = ci*tmp/k2sq
!$$$            else
!$$$               tmp = czero
!$$$            endif
!$$$  
!$$$            if (abs(k3sq) > propcut) then
!$$$               tmp = ci*tmp/k3sq
!$$$            else
!$$$               tmp = czero
!$$$            endif
!$$$          
!$$$            res = res + tmp
!$$$         
!$$$            sp4 = bfW_fbff(e,k,sp(:,1:3),p(:,1:3),fll(1:3),fl4, eW,kW,&
!$$$                 &ng1,ng2,ng3,sw,giarray,qiarray(1:3),Wid,pol_int)
!$$$  
!$$$            k4 = -(p(:,1)+p(:,2)+p(:,3)+kW)
!$$$            k4sq = sc(k4,k4)
!$$$  
!$$$            sp4 = spi2(k4,sp4)!+mass*sp4
!$$$            sp1 = sp(:,4)
!$$$            tmp = -cone*vbqq(Dv,sp1,sp4)
!$$$            
!$$$            if (abs(k4sq) > propcut) then
!$$$               tmp = ci*tmp/k4sq
!$$$            else
!$$$               tmp = czero
!$$$            endif
!$$$  
!$$$            res = res + tmp
!$$$         
!$$$         endif
!$$$      
!$$$      else
!$$$         write(*,*) 'gw_sbsfbf_2 not defined for ngluons > 0'
!$$$         stop
!$$$         
!$$$      endif
!$$$      
!$$$      ! -- store current 
!$$$      if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)
!$$$  
!$$$    end function gW_sbsfbf_3
!$$$  
!$$$  
!$$$  
!$$$    recursive function bfW_fbff(e,k,sp,p,fll,fl0,&
!$$$         &eW,kW,ng1,ng2,ng3,sw,giarray,qiarray,Wid,pol_int) result(res)
!$$$      implicit none
!$$$      double complex, intent(in) :: e(:,:), k(:,:)
!$$$      double complex, intent(in) :: sp(:,:),p(:,:)
!$$$      double complex, intent(in) :: eW(:),kW(:)
!$$$      integer, intent(in) ::  ng1,ng2,ng3,sw
!$$$      character, intent(in) :: fll(:)*3
!$$$      character, intent(in) :: fl0*3   ! flavor off-shell f-line
!$$$      integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
!$$$      ! -----------------------------------------------------------------------
!$$$      character :: flaux*3
!$$$      character :: fl1*3,fl2*3,fl3*3
!$$$      integer :: ngluon, ng4
!$$$      integer, parameter :: Ndumm=0
!$$$      double complex             :: res(size(sp,dim=1))
!$$$      double complex             :: tmp(size(sp,dim=1))
!$$$      double complex             :: k1(size(k,dim=1))
!$$$      double complex             :: k2(size(k,dim=1))
!$$$      double complex             :: k4(size(k,dim=1))
!$$$      double complex             :: sp4(size(sp,dim=1))
!$$$      double complex             :: e2(size(e,dim=1))
!$$$      double complex             :: kdumm(size(k,dim=1),Ndumm)
!$$$      double complex             :: edumm(size(e,dim=1),Ndumm)
!$$$      double complex  :: k1sq,k2sq,k4sq!,k3sq,k4sq
!$$$      !real(dp) :: mass
!$$$      logical                   :: done 
!$$$  
!$$$      !if (verbose) write(*,*) 'entering bfW_fbff'
!$$$  
!$$$      done = .false. 
!$$$  !    if (present(giarray)) then 
!$$$  !       !if (size(qiarray) /= 3) stop 'bfW_fbff: wrong size qiarray'
!$$$  !       !if (size(e,dim=2) /= size(giarray)) stop 'bfW_fbff: ng= size(giarray)' 
!$$$  !       call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!$$$  !       if (done) return 
!$$$  !    else
!$$$  !       if (i_warn < max_warn) then 
!$$$  !          write(*,*) 'bfW_fbff: giarray missing', i_warn 
!$$$  !          i_warn = i_warn+1
!$$$  !       endif
!$$$  !    endif
!$$$  
!$$$      !mass = mt
!$$$  
!$$$      ngluon = size(e,dim=2)
!$$$      ng4 = ngluon - ng1 - ng2-ng3
!$$$  
!$$$      fl1 = fll(1)
!$$$      fl2 = fll(2)
!$$$      fl3 = fll(3)
!$$$  
!$$$  
!$$$      if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C'
!$$$  
!$$$      if (ngluon == 0) then 
!$$$  
!$$$         res = czero
!$$$  
!$$$  ! XXXXXXXXX
!$$$         if ((sw.eq.3) .and. (WWqqqq .eqv. .false.) ) then 
!$$$            e2 = gW_bff(edumm,kdumm,sp(:,2),p(:,2),fl2,&
!$$$                 &sp(:,3),p(:,3),fl3,eW,kW,0,0,&
!$$$                 &giarray,qiarray(2:3),Wid,pol_int)
!$$$            k1 = p(:,2)+p(:,3)+kW
!$$$            k1sq=sc(k1,k1)
!$$$  
!$$$            if (abs(k1sq) > propcut.and.fl0==fl1) then 
!$$$               tmp = -ci/k1sq*vbqg(sp(:,1),e2)
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            res = res + tmp  ! #1
!$$$  
!$$$         endif
!$$$  
!$$$         if (sw.eq.2) then 
!$$$  
!$$$            e2 = gW_fbf(edumm,kdumm,sp(:,1),p(:,1),fl1,&
!$$$                 &sp(:,2),p(:,2),fl2,eW,kW,0,0,&
!$$$                 &giarray,qiarray(1:2),Wid,pol_int)
!$$$            k2 = p(:,1)+p(:,2)+kW 
!$$$            k2sq = sc(k2,k2)
!$$$  
!$$$            if (abs(k2sq) > propcut.and.fl0==fl3) then 
!$$$               tmp = -ci/k2sq*vgbq(e2,sp(:,3))
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            !        here tmp was multiplied with czero
!$$$            res = res + tmp   ! #2
!$$$  
!$$$         endif
!$$$  
!$$$  ! XXXXX
!$$$         if (((sw.eq.1).or.(sw.eq.2)) .and. &
!$$$              ((case_a4 .eqv. .true.) .or. (WWqqqq .eqv. .false.))) then 
!$$$           
!$$$  
!$$$            e2 = g_bff(edumm,kdumm,sp(:,2),p(:,2),fl2,&
!$$$                 &sp(:,3),p(:,3),fl3,0,0,&
!$$$                 &giarray,qiarray(2:3),pol_int)
!$$$  
!$$$            k2 = p(:,2)+p(:,3)
!$$$            k2sq = sc(k2,k2)
!$$$  
!$$$            sp4 = bfW(edumm,kdumm,sp(:,1),p(:,1),fl1,fl0,eW,kW,0,&
!$$$                 &giarray,qiarray(1:1),Wid,pol_int)
!$$$            k4  = -p(:,1) - kW
!$$$            k4sq = sc(k4,k4)
!$$$            sp4 = spi2(k4,sp4) !+ mass*sp4
!$$$  
!$$$            tmp = vbqg(sp4,e2)
!$$$  
!$$$            if (abs(k2sq) > propcut) then 
!$$$               tmp = -ci/k2sq*tmp
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            if (abs(k4sq) > propcut) then 
!$$$               tmp = ci/k4sq*tmp 
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            res = res + tmp  ! #3
!$$$  
!$$$         endif
!$$$  
!$$$  
!$$$         if ((sw.eq.3).or.(sw.eq.4)) then 
!$$$  
!$$$            e2 = g_fbf(edumm,kdumm,sp(:,1),p(:,1),fl1,&
!$$$                 & sp(:,2),p(:,2),fl2,0,0,&
!$$$                 &giarray,qiarray(1:2),pol_int)
!$$$            k2 = p(:,1)+p(:,2)
!$$$            k2sq = sc(k2,k2)
!$$$  
!$$$  
!$$$            sp4 = bfW(edumm,kdumm,sp(:,3),p(:,3),fl3,fl0,eW,kW,0,&
!$$$                 &giarray,qiarray(3:3),Wid,pol_int)
!$$$            k4  = -p(:,3) - kW
!$$$            k4sq = sc(k4,k4)
!$$$            sp4 = spi2(k4,sp4) !+ mass*sp4
!$$$  
!$$$            tmp = vgbq(e2,sp4)   
!$$$  
!$$$            if (abs(k2sq) > propcut) then 
!$$$               tmp = -ci/k2sq*tmp
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            if (abs(k4sq) > propcut) then 
!$$$               tmp = ci/k4sq*tmp 
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$  
!$$$            res = res + tmp  !# 4
!$$$  
!$$$         endif
!$$$  
!$$$         if ((sw.eq.1).or.(sw.eq.4)) then 
!$$$  
!$$$  
!$$$            if (fl0.eq.'top') flaux = 'bot'
!$$$            if (fl0.eq.'bot') flaux = 'top'
!$$$  !XXXX
!$$$            if (WWqqqq) then 
!$$$               sp4 = bf_fbff_2(edumm,kdumm,sp,p,fll,flaux,&
!$$$                    &ng1,ng2,ng3,giarray,qiarray,pol_int)          
!$$$            else
!$$$               sp4 = bf_fbff(edumm,kdumm,sp,p,fll,flaux,&
!$$$                    &ng1,ng2,ng3,giarray,qiarray,pol_int)
!$$$            endif
!$$$            k4 = -p(:,1)-p(:,2)-p(:,3)
!$$$            sp4 = spi2(k4,sp4) !+ mass*sp4
!$$$            k4sq = sc(k4,k4)
!$$$  
!$$$            tmp = vWq(eW,sp4)
!$$$  
!$$$            if (abs(k4sq) > propcut) then 
!$$$               tmp = ci/k4sq*tmp
!$$$            else 
!$$$               tmp = czero
!$$$            endif
!$$$  
!$$$  
!$$$  
!$$$            res = res + tmp   ! #5
!$$$  
!$$$         endif
!$$$  
!$$$      else  ! -- this is for ngluon > 0
!$$$  
!$$$         print *, 'nglu > 0 not done for bfW_fbff'                     
!$$$  
!$$$  
!$$$      endif
!$$$  
!$$$      ! -- store current 
!$$$      if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)
!$$$  
!$$$  
!$$$    end function bfW_fbff
!$$$  
!$$$  
!$$$    !----- this recursive function is only needed for some subleading 
!$$$    !----- color amplitudes with strange quarks 
!$$$  
!$$$    recursive function fW_bffbf_1(e,k,sp,p,fll,fl0,&
!$$$         &eW,kW,ng1,ng2,ng3,giarray,qiarray,Wid,pol_int) result(res)
!$$$      implicit none
!$$$      double complex, intent(in) :: e(:,:), k(:,:)
!$$$      double complex, intent(in) :: sp(:,:),p(:,:)
!$$$      double complex, intent(in) :: eW(:),kW(:)
!$$$      integer, intent(in) ::  ng1,ng2,ng3
!$$$      character, intent(in) :: fll(:)*3
!$$$      character, intent(in) :: fl0*3   ! flavor off-shell f-line
!$$$      integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
!$$$      ! -----------------------------------------------------------------------
!$$$      character :: flaux*3
!$$$      character :: fl1*3,fl2*3,fl3*3
!$$$      integer :: ngluon, ng4, ngL,m
!$$$      integer, parameter :: Ndumm=0
!$$$      double complex             :: res(size(sp,dim=1))
!$$$      double complex             :: tmp(size(sp,dim=1))
!$$$      double complex             :: k1(size(k,dim=1))
!$$$      double complex             :: k2(size(k,dim=1))
!$$$      double complex             :: k4(size(k,dim=1))
!$$$      double complex             :: sp1(size(sp,dim=1))
!$$$      double complex             :: sp2(size(sp,dim=1))
!$$$      double complex             :: sp4(size(sp,dim=1))
!$$$      double complex             :: e1(size(e,dim=1))
!$$$      double complex             :: e2(size(e,dim=1))
!$$$      double complex             :: kdumm(size(k,dim=1),Ndumm)
!$$$      double complex             :: edumm(size(e,dim=1),Ndumm)
!$$$      double complex  :: k1sq,k2sq,k4sq!,k3sq,k4sq
!$$$      !real(dp) :: mass
!$$$      logical                   :: done 
!$$$  
!$$$      !if (verbose) write(*,*) 'entering  fW_bffbf_1',ng1,ng2,ng3
!$$$  
!$$$      done = .false. 
!$$$  !    if (present(giarray)) then 
!$$$  !       !if (size(qiarray) /= 3) stop ' fW_bffbf_1: wrong size qiarray'
!$$$  !       !if (size(e,dim=2) /= size(giarray)) stop ' fW_bffbf_1: ng= size(giarray)' 
!$$$  !       call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!$$$  !       if (done) return 
!$$$  !    else
!$$$  !       if (i_warn < max_warn) then 
!$$$  !          write(*,*) 'fW_bffbf_1: giarray missing', i_warn 
!$$$  !          i_warn = i_warn+1
!$$$  !       endif
!$$$  !    endif
!$$$  
!$$$      !mass = mt
!$$$  
!$$$      ngluon = size(e,dim=2)
!$$$      ng4 = ngluon - ng1 - ng2-ng3
!$$$  
!$$$      fl1 = fll(1)
!$$$      fl2 = fll(2)
!$$$      fl3 = fll(3)
!$$$  
!$$$  
!$$$      if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C'
!$$$  
!$$$      if ((fl1.ne.'str').or.(fl2.ne.'str')) then 
!$$$         print *, 'flavor inconsistent' 
!$$$      endif
!$$$  
!$$$      if (ngluon == 0) then 
!$$$  
!$$$         res = czero
!$$$  
!$$$         e2 = g_bff(edumm,kdumm,sp(:,1),p(:,1),fl1,&
!$$$              & sp(:,2),p(:,2),fl2,0,0,&
!$$$              &giarray,qiarray(1:2),pol_int)
!$$$         k2 = p(:,1)+p(:,2)
!$$$         k2sq = sc(k2,k2)
!$$$  
!$$$  
!$$$         sp4 = fW(edumm,kdumm,sp(:,3),p(:,3),fl3,fl0,eW,kW,0,&
!$$$              &giarray,qiarray(3:3),Wid,pol_int)
!$$$         k4  = p(:,3) + kW
!$$$         k4sq = sc(k4,k4)
!$$$         sp4 = spb2(sp4,k4) !+ mass*sp4
!$$$  
!$$$         tmp = -ci/k2sq*ci/k4sq*vgq(e2,sp4)   
!$$$  
!$$$         res = res + tmp  !# 4
!$$$  
!$$$  
!$$$         if (fl0.eq.'top') flaux = 'bot'
!$$$         if (fl0.eq.'bot') flaux = 'top'
!$$$  
!$$$         sp4 = f_bffbf(edumm,kdumm,sp,p,fll,flaux,&
!$$$              &0,0,0,giarray,qiarray,pol_int)
!$$$  
!$$$         k4 = p(:,1)+p(:,2)+p(:,3)
!$$$         sp4 = spb2(sp4,k4) !+ mass*sp4
!$$$         k4sq = sc(k4,k4)
!$$$  
!$$$         tmp = vbqW(sp4,eW)
!$$$  
!$$$         if (abs(k4sq) > propcut) then 
!$$$            tmp = ci/k4sq*tmp
!$$$         else 
!$$$            tmp = czero
!$$$         endif
!$$$  
!$$$         res = res + tmp   ! #2
!$$$  
!$$$  
!$$$      else  ! -- this is for ngluon > 0
!$$$  
!$$$         res = czero
!$$$  
!$$$         do m=0,ng4-1
!$$$  
!$$$            ngL = ng1+ ng2+ng3+m      
!$$$  
!$$$            sp1=fW_bffbf_1(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0,&
!$$$                 &eW,kW,ng1,ng2,ng3,&
!$$$                 &giarray(1:ngL),qiarray,Wid,pol_int)
!$$$            if (1<=ngL) then 
!$$$               k1=sum(k(:,1:ngL),dim=2)
!$$$            else
!$$$               k1 = czero 
!$$$            endif
!$$$            k1 = k1 + p(:,1)+p(:,2)+p(:,3)+kW
!$$$            k1sq = sc(k1,k1) !- mass**2
!$$$  
!$$$            sp1 = spb2(sp1,k1) !+ mass*sp1
!$$$  
!$$$            e2 = vgluon(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
!$$$                 &giarray(ngL+1:ngluon),pol_int)
!$$$            if (ngL+1<=ngluon) then 
!$$$               k2=sum(k(:,ngL+1:ngluon),dim=2)
!$$$            else
!$$$               k2 = czero 
!$$$            endif
!$$$            k2sq=sc(k2,k2)
!$$$  
!$$$            if (abs(k1sq) > propcut) then 
!$$$               tmp = ci/k1sq*vqg(sp1,e2)
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            if (m < ng4-1) then 
!$$$               if (abs(k2sq) > propcut) then 
!$$$                  tmp = -ci/k2sq*tmp
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$            endif
!$$$  
!$$$            res = res + tmp
!$$$  
!$$$  
!$$$  
!$$$            tmp = czero
!$$$  
!$$$         enddo  !#1
!$$$  
!$$$  
!$$$         do m=1,ng1
!$$$  
!$$$            e1 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)
!$$$            k1=sum(k(:,1:m),dim=2)
!$$$            k1sq=sc(k1,k1)
!$$$  
!$$$            sp2=fW_bffbf_1(e(:,m+1:ngluon),k(:,m+1:ngluon),&
!$$$                 &sp,p,fll,fl0,eW,kW,ng1-m,ng2,ng3,&
!$$$                 &giarray(m+1:ngluon),qiarray,Wid,pol_int)
!$$$            if (m+1<=ngluon) then 
!$$$               k2=sum(k(:,m+1:ngluon),dim=2)
!$$$            else
!$$$               k2 = czero 
!$$$            endif
!$$$            k2 = k2 + p(:,1)+p(:,2)+p(:,3) + kW
!$$$            k2sq = sc(k2,k2) 
!$$$            sp2 = spb2(sp2,k2)!+mass*sp2
!$$$  
!$$$            tmp = ci/k2sq*vgq(e1,sp2)
!$$$  
!$$$            if (m > 1) then 
!$$$               tmp = -ci/k1sq*tmp
!$$$            endif
!$$$  
!$$$            res = res + tmp
!$$$  
!$$$  
!$$$  
!$$$            tmp = czero
!$$$  
!$$$         enddo  !#2
!$$$  
!$$$  ! XXX
!$$$         if ((case_a2 .eqv. .false.)) then
!$$$            do m=0,ng3
!$$$               
!$$$               ngL = ng1+ ng2+m      
!$$$  
!$$$  
!$$$            e1 = g_bff(e(:,1:ngL),k(:,1:ngL),sp(:,1),p(:,1),fl1,&
!$$$                 &sp(:,2),p(:,2),fl2,ng1,ng2,&
!$$$                 &giarray(1:ngL),qiarray(1:2),pol_int)
!$$$            if (1<=ngL) then 
!$$$               k1=sum(k(:,1:ngL),dim=2)+p(:,1)+p(:,2)
!$$$            else
!$$$               k1=p(:,1)+p(:,2)
!$$$            endif
!$$$            k1sq=sc(k1,k1)
!$$$  
!$$$  
!$$$  
!$$$            sp2=fW(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
!$$$                 &sp(:,3),p(:,3),fl3,fl0,eW,kW,ng3-m,&
!$$$                 &giarray(ngL+1:ngluon),qiarray(3:3),Wid,pol_int)
!$$$            if (ngL+1<=ngluon) then 
!$$$               k2=sum(k(:,ngL+1:ngluon),dim=2)
!$$$            else
!$$$               k2 = czero 
!$$$            endif
!$$$            k2 = k2 + p(:,3) + kW
!$$$            k2sq = sc(k2,k2) !- mass**2
!$$$  
!$$$  
!$$$            sp2 = spb2(sp2,k2)!+mass*sp2
!$$$  
!$$$            tmp = -ci/k1sq*ci/k2sq*vgq(e1,sp2)
!$$$  
!$$$            res = res + tmp
!$$$  
!$$$  
!$$$  
!$$$            tmp = czero
!$$$  
!$$$  
!$$$         enddo  !#3
!$$$      endif
!$$$  
!$$$         if (fl0.eq.'top') flaux = 'bot'
!$$$         if (fl0.eq.'bot') flaux = 'top'
!$$$  
!$$$  
!$$$         sp4 = f_bffbf(e,k,sp,p,fll,flaux,&
!$$$              &ng1,ng2,ng3,giarray,qiarray,pol_int)
!$$$         if (1<=ngluon) then 
!$$$            k4 = sum(k(:,1:ngluon),dim=2)+p(:,1)+p(:,2)+p(:,3)
!$$$         else
!$$$            k4 = p(:,1)+p(:,2)+p(:,3)
!$$$         endif
!$$$         sp4 = spb2(sp4,k4) !+ mass*sp4
!$$$         k4sq = sc(k4,k4)
!$$$  
!$$$         tmp = ci/k4sq*vbqW(sp4,eW)
!$$$  
!$$$  
!$$$         res = res + tmp   ! #4
!$$$  
!$$$  
!$$$  
!$$$         tmp = czero
!$$$  
!$$$      endif
!$$$  
!$$$      ! -- store current 
!$$$      if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)
!$$$  
!$$$  
!$$$    end function fW_bffbf_1
!$$$  
!$$$  
!$$$  
!$$$  !---------------------
!$$$  ! This function is needed forcase A2 and A4, where there is a need to 
!$$$  ! have quarks 1-2 on the same line, but also keep track of where the W is.
!$$$  !-------------------------
!$$$  
!$$$    recursive function fW_bffbf_2(e,k,sp,p,fll,fl0,&
!$$$         &eW,kW,ng1,ng2,ng3,sw,giarray,qiarray,Wid,pol_int) result(res)
!$$$      double complex, intent(in) :: e(:,:), k(:,:)
!$$$      double complex, intent(in) :: sp(:,:),p(:,:)
!$$$      double complex, intent(in) :: eW(:),kW(:)
!$$$      integer, intent(in) ::  ng1,ng2,ng3,sw
!$$$      character, intent(in) :: fll(:)*3
!$$$      character, intent(in) :: fl0*3   ! flavor off-shell f-line
!$$$      integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
!$$$      ! -----------------------------------------------------------------------
!$$$      character :: flaux*3
!$$$      character :: fl1*3,fl2*3,fl3*3
!$$$  !    integer             :: m2,m3, ms1a, ms2a
!$$$      integer :: ngluon, ng4, ngL,m!,m1
!$$$      integer, parameter :: Ndumm=0
!$$$      double complex             :: res(size(sp,dim=1))
!$$$      double complex             :: tmp(size(sp,dim=1))
!$$$      double complex             :: k1(size(k,dim=1))
!$$$      double complex             :: k2(size(k,dim=1))
!$$$      double complex             :: k4(size(k,dim=1))
!$$$      double complex             :: sp1(size(sp,dim=1))
!$$$      double complex             :: sp2(size(sp,dim=1))
!$$$      double complex             :: sp4(size(sp,dim=1))
!$$$      double complex             :: e1(size(e,dim=1))
!$$$      double complex             :: e2(size(e,dim=1))
!$$$      double complex             :: kdumm(size(k,dim=1),Ndumm)
!$$$      double complex             :: edumm(size(e,dim=1),Ndumm)
!$$$      double complex  :: k1sq,k2sq,k4sq
!$$$      !real(dp) :: mass,mass2
!$$$      logical                   :: done 
!$$$  
!$$$      !if (verbose) write(*,*) 'entering fW_bffbf'
!$$$  
!$$$      done = .false. 
!$$$      if (present(giarray)) then 
!$$$         !if (size(qiarray) /= 3) stop 'fW_bffbf: wrong size qiarray'
!$$$         !if (size(e,dim=2) /= size(giarray)) stop 'fW_bffbf: ng= size(giarray)' 
!$$$         ! XXX 
!$$$        call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!$$$         if (done) return 
!$$$      else
!$$$         write(*,*) 'giarray missing'
!$$$      endif
!$$$  
!$$$      !mass = mt
!$$$      !mass2 = mass**2
!$$$  
!$$$      ngluon = size(e,dim=2)
!$$$      ng4 = ngluon - ng1 - ng2-ng3
!$$$      
!$$$  
!$$$      fl1 = fll(1)
!$$$      fl2 = fll(2)
!$$$      fl3 = fll(3)
!$$$  
!$$$  
!$$$      if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C: fW_bffbf', ngluon, ng1, ng2, ng3
!$$$      
!$$$      if (ngluon == 0) then 
!$$$  
!$$$         res = czero
!$$$  
!$$$  
!$$$         if (sw.eq.2)  then 
!$$$  
!$$$            e2 = gW_bff(edumm,kdumm,sp(:,1),p(:,1),fl1,&
!$$$                 &sp(:,2),p(:,2),fl2,eW,kW,0,0,&
!$$$                 &giarray,qiarray(1:2),Wid,pol_int)    
!$$$            k2 = p(:,1)+p(:,2)+kW 
!$$$            k2sq = sc(k2,k2)
!$$$  
!$$$            if (abs(k2sq) > propcut.and.fl0==fl3) then 
!$$$               tmp = -ci/k2sq*vgq(e2,sp(:,3))
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            !        used to have czero here (now)
!$$$            res = res + tmp                 
!$$$  
!$$$         endif
!$$$  
!$$$  
!$$$         if ((sw.eq.3).or.(sw.eq.4)) then 
!$$$            
!$$$            e2 = g_bff(edumm,kdumm,sp(:,1),p(:,1),fl1,&
!$$$                 & sp(:,2),p(:,2),fl2,0,0,&
!$$$                 &giarray,qiarray(1:2),pol_int)    
!$$$           
!$$$            k2 = p(:,1)+p(:,2)
!$$$            k2sq = sc(k2,k2)
!$$$       
!$$$  
!$$$            sp4 = fW(e,k,sp(:,3),p(:,3),fl3,fl0,eW,kW,0,&
!$$$                 &giarray,qiarray(3:3),Wid,pol_int)    
!$$$            k4  = p(:,3) + kW
!$$$            k4sq = sc(k4,k4)
!$$$            sp4 = spb2(sp4,k4) !+ mass*sp4
!$$$  
!$$$            tmp = vgq(e2,sp4)  
!$$$         
!$$$  
!$$$            if (abs(k2sq) > propcut) then 
!$$$               tmp = -ci/k2sq*tmp
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            if (abs(k4sq) > propcut) then 
!$$$               tmp = ci/k4sq*tmp 
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            res = res + tmp  
!$$$  
!$$$         endif
!$$$  
!$$$         if ((sw.eq.1).or.(sw.eq.4)) then 
!$$$     
!$$$            if (fl0.eq.'top') flaux = 'bot'
!$$$            if (fl0.eq.'bot') flaux = 'top'
!$$$  
!$$$           
!$$$            sp4 = f_bffbf_3(edumm,kdumm,sp,p,fll,flaux,&
!$$$                 &ng1,ng2,ng3,&
!$$$                 &giarray,qiarray,pol_int)  
!$$$            
!$$$      
!$$$            k4 = p(:,1)+p(:,2)+p(:,3)
!$$$            sp4 = spb2(sp4,k4) !+ mass*sp4
!$$$            k4sq = sc(k4,k4)
!$$$  
!$$$            tmp = vbqW(sp4,eW)
!$$$            
!$$$  
!$$$            if (abs(k4sq) > propcut) then 
!$$$               tmp = ci/k4sq*tmp
!$$$            else 
!$$$               tmp = czero
!$$$            endif
!$$$  
!$$$            res = res + tmp   
!$$$         endif
!$$$         
!$$$    
!$$$         else  ! -- this is for ngluon > 0
!$$$    
!$$$            res = czero
!$$$         
!$$$            !--------- next step is valid for all sw
!$$$  
!$$$            do m=0,ng4-1
!$$$  
!$$$               ngL = ng1+ ng2+ng3+m      
!$$$            
!$$$               sp1=fW_bffbf_2(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0,&
!$$$                    &eW,kW,ng1,ng2,ng3,sw,&
!$$$                    &giarray(1:ngL),qiarray,Wid,pol_int)    
!$$$               if (1<=ngL) then 
!$$$                  k1=sum(k(:,1:ngL),dim=2)
!$$$               else
!$$$                  k1 = czero 
!$$$               endif
!$$$               k1 = k1 + p(:,1)+p(:,2)+p(:,3)+kW
!$$$               k1sq = sc(k1,k1) !- mass2
!$$$  
!$$$               sp1 = spb2(sp1,k1) !+ mass*sp1
!$$$               
!$$$               e2 = vgluon(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
!$$$                    &giarray(ngL+1:ngluon),pol_int)    
!$$$               if (ngL+1<=ngluon) then 
!$$$                  k2=sum(k(:,ngL+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2sq=sc(k2,k2)
!$$$               
!$$$               if (abs(k1sq) > propcut) then 
!$$$                  tmp = ci/k1sq*vqg(sp1,e2)
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$  
!$$$               if (m < ng4-1) then 
!$$$                  if (abs(k2sq) > propcut) then 
!$$$                     tmp = -ci/k2sq*tmp
!$$$                  else 
!$$$                     tmp = czero 
!$$$                  endif
!$$$               endif
!$$$               
!$$$               res = res + tmp
!$$$  
!$$$            enddo  !#1
!$$$  
!$$$         !-------- next step is valid for any sw
!$$$  
!$$$            do m=1,ng1
!$$$  
!$$$               e1 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)    
!$$$               k1=sum(k(:,1:m),dim=2)
!$$$               k1sq=sc(k1,k1)
!$$$               
!$$$               sp2=fW_bffbf_2(e(:,m+1:ngluon),k(:,m+1:ngluon),&
!$$$                    &sp,p,fll,fl0,eW,kW,ng1-m,ng2,ng3,sw,&
!$$$                    &giarray(m+1:ngluon),qiarray,Wid,pol_int)    
!$$$               if (m+1<=ngluon) then 
!$$$                  k2=sum(k(:,m+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,1)+p(:,2)+p(:,3)+kW
!$$$               k2sq = sc(k2,k2) !- mass2
!$$$               sp2 = spb2(sp2,k2)!+mass*sp2
!$$$  
!$$$               if (abs(k2sq) > propcut) then 
!$$$                  tmp = ci/k2sq*vgq(e1,sp2)
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$  
!$$$               if (m > 1) then 
!$$$                  if (abs(k1sq) > propcut) then 
!$$$                     tmp = -ci/k1sq*tmp
!$$$                  else 
!$$$                     tmp = czero 
!$$$                  endif
!$$$               endif
!$$$  
!$$$               res = res + tmp
!$$$               
!$$$           
!$$$            enddo  !#2
!$$$  
!$$$            if ((sw.eq.3).or.(sw.eq.4)) then 
!$$$           
!$$$               do m=0,ng3
!$$$                  
!$$$                  ngL = ng1+ ng2+m      
!$$$                  
!$$$                  e1 = g_bff(e(:,1:ngL),k(:,1:ngL),sp(:,1),p(:,1),fl1,&
!$$$                       &sp(:,2),p(:,2),fl2,ng1,ng2,&
!$$$                       &giarray(1:ngL),qiarray(1:2),pol_int)    
!$$$                  if (1<=ngL) then 
!$$$                     k1=sum(k(:,1:ngL),dim=2)+p(:,1)+p(:,2)
!$$$                  else
!$$$                     k1 = p(:,1)+p(:,2)
!$$$                  endif
!$$$                  k1sq=sc(k1,k1)
!$$$  
!$$$                  sp2=fW(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
!$$$                       &sp(:,3),p(:,3),fl3,fl0,eW,kW,ng3-m,&
!$$$                       &giarray(ngL+1:ngluon),qiarray(3:3),Wid,pol_int)    
!$$$                  if (ngL+1<=ngluon) then 
!$$$                     k2=sum(k(:,ngL+1:ngluon),dim=2)
!$$$                  else
!$$$                     k2 = czero 
!$$$                  endif
!$$$                  k2 = k2 + p(:,3)+kW
!$$$                  k2sq = sc(k2,k2) !- mass2
!$$$                  
!$$$                  sp2 = spb2(sp2,k2)!+mass*sp2
!$$$  
!$$$                  if (abs(k1sq) > propcut) then 
!$$$                     tmp = -ci/k1sq*vgq(e1,sp2)
!$$$                  else 
!$$$                     tmp = czero 
!$$$                  endif
!$$$                  
!$$$                  if (abs(k2sq) > propcut) then 
!$$$                     tmp = ci/k2sq*tmp
!$$$                  else 
!$$$                     tmp = czero 
!$$$                  endif
!$$$  
!$$$                  res = res + tmp
!$$$              
!$$$               enddo  !#3
!$$$               
!$$$            endif
!$$$    
!$$$  
!$$$            if ((sw.eq.1).or.(sw.eq.4)) then 
!$$$  
!$$$               if (fl0.eq.'top') flaux = 'bot'
!$$$               if (fl0.eq.'bot') flaux = 'top'
!$$$               
!$$$               sp4 = f_bffbf_3(e,k,sp,p,fll,flaux,&
!$$$                    &ng1,ng2,ng3,&
!$$$                    &giarray,qiarray,pol_int)    
!$$$  
!$$$               if (1<=ngluon) then 
!$$$                  k4 = sum(k(:,1:ngluon),dim=2)+p(:,1)+p(:,2)+p(:,3)
!$$$               else
!$$$                  k4 = p(:,1)+p(:,2)+p(:,3)
!$$$               endif
!$$$  
!$$$               sp4 = spb2(sp4,k4) !+ mass*sp4
!$$$               k4sq = sc(k4,k4)
!$$$               
!$$$               tmp = vbqW(sp4,eW)
!$$$               
!$$$               if (abs(k4sq) > propcut) then 
!$$$                  tmp = ci/k4sq*tmp
!$$$               else 
!$$$                  tmp = czero
!$$$               endif
!$$$               
!$$$               res = res + tmp   ! #4
!$$$      
!$$$            endif
!$$$            
!$$$         endif                   
!$$$      
!$$$      ! -- store current 
!$$$      if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)
!$$$  
!$$$  
!$$$    end function fW_bffbf_2
!$$$  
!$$$  
!$$$  
!$$$  
!$$$  
!$$$  
!$$$    !--------------------------------------------------------------------------
!$$$    !------ this is the recursive function with 
!$$$    !-------strange quarks ; for strange quark loops 
!$$$    !------ only 
!$$$  
!$$$    recursive function fsW_fbfbf(e,k,sp,p,fll,fl0,&
!$$$         &eW,kW,ng1,ng2,ng3,sw,giarray,qiarray,Wid,pol_int) result(res)
!$$$      double complex, intent(in) :: e(:,:), k(:,:)
!$$$      double complex, intent(in) :: sp(:,:),p(:,:)
!$$$      double complex, intent(in) :: eW(:),kW(:)
!$$$      integer, intent(in) ::  ng1,ng2,ng3,sw
!$$$      character, intent(in) :: fll(:)*3
!$$$      character, intent(in) :: fl0*3   ! flavor off-shell f-line
!$$$      integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
!$$$      ! -----------------------------------------------------------------------
!$$$      character :: fl1*3,fl2*3,fl3*3
!$$$      integer :: ngluon, ng4, ngL,m
!$$$      integer, parameter :: Ndumm=0
!$$$      double complex             :: res(size(sp,dim=1))
!$$$      double complex             :: res_stored(size(sp,dim=1))
!$$$      double complex             :: tmp(size(sp,dim=1))
!$$$      double complex             :: k1(size(k,dim=1))
!$$$      double complex             :: k2(size(k,dim=1))
!$$$      double complex             :: sp1(size(sp,dim=1))
!$$$      double complex             :: sp2(size(sp,dim=1))
!$$$      double complex             :: e1(size(e,dim=1))
!$$$      double complex             :: e2(size(e,dim=1))
!$$$      double complex             :: kdumm(size(k,dim=1),Ndumm)
!$$$      double complex             :: edumm(size(e,dim=1),Ndumm)
!$$$      double complex  :: k1sq,k2sq
!$$$      !real(dp) :: mass
!$$$      logical                   :: done 
!$$$  
!$$$      !if (verbose) write(*,*) 'entering  fsW_fbfbf',ng1,ng2,ng3
!$$$  
!$$$      done = .false. 
!$$$  !    if (present(giarray)) then 
!$$$  !       !if (size(qiarray) /= 3) stop ' fsW_fbfbf: wrong size qiarray'
!$$$  !       !if (size(e,dim=2) /= size(giarray)) stop ' fsW_fbfbf: ng= size(giarray)' 
!$$$  !       call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!$$$  !       if (done) return 
!$$$  !    else
!$$$  !       if (i_warn < max_warn) then 
!$$$  !          write(*,*) 'fsW_fbfbf: giarray missing', i_warn 
!$$$  !          i_warn = i_warn+1
!$$$  !       endif
!$$$  !    endif
!$$$  
!$$$      !mass = mt
!$$$  
!$$$  
!$$$      ngluon = size(e,dim=2)
!$$$      ng4 = ngluon - ng1 - ng2-ng3
!$$$  
!$$$      fl1 = fll(1)
!$$$      fl2 = fll(2)
!$$$      fl3 = fll(3)
!$$$  
!$$$  
!$$$      if (ng4 < 0) write(*,*) 'ERROR IN CURRENT S'
!$$$  
!$$$      if (ngluon == 0) then 
!$$$  
!$$$         res = czero
!$$$  
!$$$         e2 = gW_fbf(edumm,kdumm,sp(:,1),p(:,1),fl1,&
!$$$              &sp(:,2),p(:,2),fl2,eW,kW,0,0,&
!$$$              &giarray,qiarray(1:2),Wid,pol_int)
!$$$         k1 = p(:,1)+p(:,2)+kW
!$$$         k1sq=sc(k1,k1)
!$$$  
!$$$  
!$$$  
!$$$         if (abs(k1sq) > propcut.and.fl0==fl3) then 
!$$$            tmp = -ci/k1sq*vgq(e2,sp(:,3))
!$$$         else 
!$$$            tmp = czero 
!$$$         endif
!$$$  
!$$$         res = res + tmp  ! #1
!$$$  
!$$$  
!$$$      else  ! -- this is for ngluon > 0
!$$$  
!$$$  
!$$$         res = czero                 
!$$$  
!$$$         do m=0,ng4-1
!$$$  
!$$$            ngL = ng1+ ng2+ng3+m      
!$$$  
!$$$            sp1=fsW_fbfbf(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0,&
!$$$                 &eW,kW,ng1,ng2,ng3,sw,&
!$$$                 &giarray(1:ngL),qiarray,Wid,pol_int)
!$$$            if (1<=ngL) then 
!$$$               k1=sum(k(:,1:ngL),dim=2)
!$$$            else
!$$$               k1 = czero 
!$$$            endif
!$$$            k1 = k1 + p(:,1)+p(:,2)+p(:,3)+kW
!$$$            k1sq = sc(k1,k1) !-mass*2
!$$$  
!$$$            sp1 = spb2(sp1,k1) !+ mass*sp1
!$$$  
!$$$            e2 = vgluon(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
!$$$                 &giarray(ngL+1:ngluon),pol_int)
!$$$            if (ngL+1<=ngluon) then 
!$$$               k2=sum(k(:,ngL+1:ngluon),dim=2)
!$$$            else
!$$$               k2 = czero 
!$$$            endif
!$$$            k2sq=sc(k2,k2)
!$$$  
!$$$            if (abs(k1sq) > propcut) then 
!$$$               tmp = ci/k1sq*vqg(sp1,e2)
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            if (m < ng4-1) then 
!$$$               if (abs(k2sq) > propcut) then 
!$$$                  tmp = -ci/k2sq*tmp
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$            endif
!$$$  
!$$$            res = res + tmp
!$$$  
!$$$  
!$$$         enddo  !#1
!$$$  
!$$$  
!$$$         do m=1,ng1
!$$$  
!$$$            e1 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)
!$$$            if (1<=m) then 
!$$$               k1=sum(k(:,1:m),dim=2)
!$$$            else
!$$$               k1 = czero 
!$$$            endif
!$$$            k1sq=sc(k1,k1)
!$$$  
!$$$            sp2=fsW_fbfbf(e(:,m+1:ngluon),k(:,m+1:ngluon),&
!$$$                 &sp,p,fll,fl0,eW,kW,ng1-m,ng2,ng3,sw,&
!$$$                 &giarray(m+1:ngluon),qiarray,Wid,pol_int)
!$$$            if (m+1<=ngluon) then 
!$$$               k2=sum(k(:,m+1:ngluon),dim=2)
!$$$            else
!$$$               k2 = czero 
!$$$            endif
!$$$            k2 = k2 + p(:,1)+p(:,2)+p(:,3)+kW
!$$$            k2sq = sc(k2,k2) !- mass**2
!$$$            sp2 = spb2(sp2,k2)!+mass*sp2
!$$$  
!$$$  
!$$$            if (abs(k2sq) > propcut) then 
!$$$               tmp = ci/k2sq*vgq(e1,sp2)
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$  
!$$$            if (m > 1) then 
!$$$               if (abs(k1sq) > propcut) then 
!$$$                  tmp = -ci/k1sq*tmp
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$            endif
!$$$  
!$$$            res = res + tmp
!$$$  
!$$$  
!$$$         enddo  !#2
!$$$  
!$$$  
!$$$         do m=0,ng3
!$$$  
!$$$            ngL = ng1+ ng2+m      
!$$$  
!$$$            e1 = gW_fbf(e(:,1:ngL),k(:,1:ngL),sp(:,1),p(:,1),fl1,&
!$$$                 &sp(:,2),p(:,2),fl2,eW,kW,ng1,ng2,&
!$$$                 &giarray(1:ngL),qiarray(1:2),Wid,pol_int)
!$$$            if (1<=ngL) then 
!$$$               k1=sum(k(:,1:ngL),dim=2)+p(:,1)+p(:,2)+kW
!$$$            else
!$$$               k1 = p(:,1)+p(:,2)+kW
!$$$            endif
!$$$            k1sq=sc(k1,k1)
!$$$  
!$$$  
!$$$            sp2=f(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
!$$$                 &sp(:,3),p(:,3),fl3,fl0,ng3-m,&
!$$$                 &giarray(ngL+1:ngluon),qiarray(3:3),pol_int)
!$$$            if (ngL+1<=ngluon) then 
!$$$               k2=sum(k(:,ngL+1:ngluon),dim=2)
!$$$            else
!$$$               k2 = czero 
!$$$            endif
!$$$            k2 = k2 + p(:,3)
!$$$            k2sq = sc(k2,k2) !- mass**2
!$$$  
!$$$            if (ng4 > 0.or. m < ng3) then 
!$$$               sp2 = spb2(sp2,k2)!+mass*sp2
!$$$            endif
!$$$  
!$$$  
!$$$            if (abs(k1sq) > propcut) then 
!$$$               tmp = -ci/k1sq*vgq(e1,sp2)
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            if (ng4 > 0.or.m < ng3) then 
!$$$               if (abs(k2sq) > propcut) then 
!$$$                  tmp = ci/k2sq*tmp
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$            endif
!$$$  
!$$$  
!$$$            res = res + tmp
!$$$  
!$$$  
!$$$  
!$$$         enddo  !#3
!$$$  
!$$$  
!$$$      endif
!$$$  
!$$$      !if (verbose) write(*,*) 'done fsW_fbfbf',ng1,ng2,ng3,pol_int, 'g',giarray,'q',qiarray, res
!$$$  
!$$$      ! -- store current 
!$$$      if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)
!$$$  
!$$$  
!$$$    end function fsW_fbfbf
!$$$  
!$$$    ! -- new for ferm loops with Z
!$$$  
!$$$    recursive function fW_fbfbf(e,k,sp,p,fll,fl0,&
!$$$         &eW,kW,ng1,ng2,ng3,sw,giarray,qiarray,Wid,pol_int) result(res)
!$$$      double complex, intent(in) :: e(:,:), k(:,:)
!$$$      double complex, intent(in) :: sp(:,:),p(:,:)
!$$$      double complex, intent(in) :: eW(:),kW(:)
!$$$      integer, intent(in) ::  ng1,ng2,ng3,sw
!$$$      character, intent(in) :: fll(:)*3
!$$$      character, intent(in) :: fl0*3   ! flavor off-shell f-line
!$$$      integer, intent(in), optional       :: giarray(:),qiarray(:),Wid,pol_int 
!$$$      ! -----------------------------------------------------------------------
!$$$      character :: flaux*3
!$$$      character :: fl1*3,fl2*3,fl3*3
!$$$      integer :: ngluon, ng4, ngL,m!,m1
!$$$      integer, parameter :: Ndumm=0
!$$$      double complex             :: res(size(sp,dim=1))
!$$$      double complex             :: tmp(size(sp,dim=1))
!$$$      double complex             :: k1(size(k,dim=1))
!$$$      double complex             :: k2(size(k,dim=1))
!$$$      double complex             :: k4(size(k,dim=1))
!$$$      double complex             :: sp1(size(sp,dim=1))
!$$$      double complex             :: sp2(size(sp,dim=1))
!$$$      double complex             :: sp4(size(sp,dim=1))
!$$$      double complex             :: e1(size(e,dim=1))
!$$$      double complex             :: e2(size(e,dim=1))
!$$$      double complex             :: kdumm(size(k,dim=1),Ndumm)
!$$$      double complex             :: edumm(size(e,dim=1),Ndumm)
!$$$      double complex             :: k1sq,k2sq,k4sq
!$$$      !real(dp)                :: mass,mass2
!$$$      logical                   :: done 
!$$$  
!$$$      !if (verbose) write(*,*) 'entering fW_fbfbf',ng1,ng2,ng3,sw,fl0,fll
!$$$  
!$$$      done = .false. 
!$$$  !    if (present(giarray)) then 
!$$$  !       !if (size(qiarray) /= 3) stop 'fW_fbfbf: wrong size qiarray'
!$$$  !       !if (size(e,dim=2) /= size(giarray)) stop 'fW_fbfbf: ng= size(giarray)' 
!$$$  !       call memory_check(pol_int,res,done,giarray,qiarray,Wid)
!$$$  !       if (done) return 
!$$$  !    else
!$$$  !       if (i_warn < max_warn) then 
!$$$  !          write(*,*) 'fW_fbfbf: giarray missing', i_warn 
!$$$  !          i_warn = i_warn+1
!$$$  !       endif
!$$$  !    endif
!$$$  
!$$$      !mass = mt
!$$$      !mass2 = mass**2
!$$$  
!$$$      ngluon = size(e,dim=2)
!$$$      ng4 = ngluon - ng1 - ng2-ng3
!$$$      fl1 = fll(1)
!$$$      fl2 = fll(2)
!$$$      fl3 = fll(3)
!$$$  
!$$$  
!$$$      if (ng4 > 0) stop 'ERROR IN CURRENT fW_fbfbf ng4 > 0'
!$$$      if ((ng1 < 0) .or. (ng2 < 0) .or. (ng3 < 0) .or. (ng4 < 0)) &
!$$$           &stop 'ERROR IN CURRENT fW_fbfbf: some ng<0'
!$$$      if ((sw .eq. 2) .or. (sw .eq.4)) stop 'ERROR IN CURRENT fW_fbfbf sw/=1,3'
!$$$      if (fl0.eq.'top' .or. fl0.eq.'bot') stop 'ERROR IN CURRENT fW_fbfbf: wrong flavour' 
!$$$  
!$$$      if (ngluon == 0) then 
!$$$  
!$$$         res = czero
!$$$         if (sw.eq.3) then 
!$$$  
!$$$            e1 = g_fbf(edumm,kdumm,sp(:,1),p(:,1),fl1,&
!$$$                 &sp(:,2),p(:,2),fl2,0,0,&
!$$$                 &giarray,qiarray(1:2),pol_int)    
!$$$            k1 = p(:,1)+p(:,2)
!$$$            k1sq=sc(k1,k1)
!$$$            sp2 = fW(edumm,kdumm,sp(:,3),p(:,3),fl3,fl0,&
!$$$                 &eW,kW,0,giarray,qiarray(3:3),Wid,pol_int)
!$$$            k2  = p(:,3) + kW
!$$$            k2sq = sc(k2,k2)
!$$$            sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$  
!$$$            tmp = vgq(e1,sp2)
!$$$            if (abs(k1sq) > propcut) then 
!$$$               tmp = -ci/k1sq*tmp 
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            if (abs(k2sq) > propcut) then 
!$$$               tmp = ci/k2sq*tmp
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            res = res + tmp  ! #1
!$$$         endif
!$$$  
!$$$  
!$$$         if ((sw.eq.1)) then 
!$$$  
!$$$  
!$$$            sp4 = f_fbfbf_3(edumm,kdumm,sp,p,fll,fl0,&
!$$$                 &ng1,ng2,ng3,giarray,qiarray,pol_int)
!$$$            
!$$$            k4 = p(:,1)+p(:,2)+p(:,3)
!$$$            sp4 = spb2(sp4,k4) !+ mass*sp4
!$$$            k4sq = sc(k4,k4)
!$$$  
!$$$            tmp = vbqW(sp4,eW)
!$$$  
!$$$            if (abs(k4sq) > propcut) then 
!$$$               tmp = ci/k4sq*tmp
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            res = res + tmp  ! #3
!$$$  
!$$$         endif
!$$$  
!$$$  
!$$$      else  ! -- this is for ngluon > 0
!$$$  
!$$$         res = czero
!$$$  
!$$$         do m=0,ng4-1
!$$$  
!$$$            ngL = ng1+ ng2+ng3+m      
!$$$  
!$$$            sp1=fW_fbfbf(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0,&
!$$$                 &eW,kW,ng1,ng2,ng3,sw,&
!$$$                 &giarray(1:ngL),qiarray,Wid,pol_int)    
!$$$            if (1<=ngL) then 
!$$$               k1=sum(k(:,1:ngL),dim=2)
!$$$            else
!$$$               k1 = czero 
!$$$            endif
!$$$            k1 = k1 + p(:,1)+p(:,2)+p(:,3)+kW
!$$$            k1sq = sc(k1,k1) !- mass2
!$$$  
!$$$            sp1 = spb2(sp1,k1) !+ mass*sp1
!$$$  
!$$$            e2 = vgluon(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
!$$$                 &giarray(ngL+1:ngluon),pol_int)    
!$$$            if (ngL+1<=ngluon) then 
!$$$               k2=sum(k(:,ngL+1:ngluon),dim=2)
!$$$            else
!$$$               k2 = czero 
!$$$            endif
!$$$            k2sq=sc(k2,k2)
!$$$  
!$$$            if (abs(k1sq) > propcut) then 
!$$$               tmp = ci/k1sq*vqg(sp1,e2)
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            if (m < ng4-1) then 
!$$$               if (abs(k2sq) > propcut) then 
!$$$                  tmp = -ci/k2sq*tmp
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$            endif
!$$$  
!$$$            res = res + tmp
!$$$  
!$$$  
!$$$         enddo  !#2
!$$$  
!$$$  
!$$$         !-------- next step is valid for any sw
!$$$  
!$$$         do m=1,ng1
!$$$  
!$$$            e1 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)    
!$$$            k1=sum(k(:,1:m),dim=2)
!$$$            k1sq=sc(k1,k1)
!$$$  
!$$$            sp2=fW_fbfbf(e(:,m+1:ngluon),k(:,m+1:ngluon),&
!$$$                 &sp,p,fll,fl0,eW,kW,ng1-m,ng2,ng3,sw,&
!$$$                 &giarray(m+1:ngluon),qiarray,Wid,pol_int)    
!$$$            if (m+1<=ngluon) then 
!$$$               k2=sum(k(:,m+1:ngluon),dim=2)
!$$$            else
!$$$               k2 = czero 
!$$$            endif
!$$$            k2 = k2 + p(:,1)+p(:,2)+p(:,3)+kW
!$$$            k2sq = sc(k2,k2) !- mass2
!$$$            sp2 = spb2(sp2,k2)!+mass*sp2
!$$$  
!$$$  
!$$$            if (abs(k2sq) > propcut) then 
!$$$               tmp = ci/k2sq*vgq(e1,sp2)
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$  
!$$$            if (m > 1) then 
!$$$               if (abs(k1sq) > propcut) then 
!$$$                  tmp = -ci/k1sq*tmp
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$            endif
!$$$  
!$$$            res = res + tmp
!$$$  
!$$$  
!$$$         enddo  !#3
!$$$  
!$$$  
!$$$         if (sw.eq.3) then 
!$$$  
!$$$            do m=0,ng3
!$$$  
!$$$               ngL = ng1+ ng2+m      
!$$$  
!$$$               e1 = g_fbf(e(:,1:ngL),k(:,1:ngL),sp(:,1),p(:,1),fl1,&
!$$$                    &sp(:,2),p(:,2),fl2,ng1,ng2,&
!$$$                    &giarray(1:ngL),qiarray(1:2),pol_int)    
!$$$               if (1<=ngL) then 
!$$$                  k1=sum(k(:,1:ngL),dim=2)+p(:,1)+p(:,2)
!$$$               else
!$$$                  k1 = p(:,1)+p(:,2)
!$$$               endif
!$$$               k1sq=sc(k1,k1)
!$$$  
!$$$  
!$$$               sp2=fW(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
!$$$                    &sp(:,3),p(:,3),fl3,fl0,eW,kW,ng3-m,&
!$$$                    &giarray(ngL+1:ngluon),qiarray(3:3),Wid,pol_int)    
!$$$               if (ngL+1<=ngluon) then 
!$$$                  k2=sum(k(:,ngL+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,3)+kW
!$$$               k2sq = sc(k2,k2) !- mass2
!$$$  
!$$$  
!$$$               sp2 = spb2(sp2,k2)!+mass*sp2
!$$$  
!$$$  
!$$$               if (abs(k1sq) > propcut) then 
!$$$                  tmp = -ci/k1sq*vgq(e1,sp2)
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$  
!$$$  
!$$$  
!$$$  
!$$$               if (abs(k2sq) > propcut) then 
!$$$                  tmp = ci/k2sq*tmp
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$  
!$$$               !----------used to have czero here (now)
!$$$               res = res + tmp
!$$$  
!$$$            enddo  !#4
!$$$  
!$$$         endif
!$$$  
!$$$  
!$$$  
!$$$         if ((sw.eq.1)) then 
!$$$  
!$$$            if (fl0.eq.'chr') flaux = 'chr'
!$$$  
!$$$            sp4 = f_fbfbf_3(e,k,sp,p,fll,flaux,&
!$$$                 &ng1,ng2,ng3,&
!$$$                 &giarray,qiarray,pol_int)    
!$$$  
!$$$  
!$$$            if (1<=ngluon) then 
!$$$               k4 = sum(k(:,1:ngluon),dim=2)+p(:,1)+p(:,2)+p(:,3)
!$$$            else
!$$$               k4 = p(:,1)+p(:,2)+p(:,3)
!$$$            endif
!$$$  
!$$$            sp4 = spb2(sp4,k4) !+ mass*sp4
!$$$            k4sq = sc(k4,k4)
!$$$  
!$$$            tmp = vbqW(sp4,eW)
!$$$  
!$$$            if (abs(k4sq) > propcut) then 
!$$$               tmp = ci/k4sq*tmp
!$$$            else 
!$$$               tmp = czero
!$$$            endif
!$$$  
!$$$            res = res + tmp   ! #7
!$$$  
!$$$         endif
!$$$  
!$$$      endif
!$$$  
!$$$      ! -- store current 
!$$$      if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)
!$$$  
!$$$  
!$$$    end function fW_fbfbf
!$$$  
!$$$  recursive function fW_bffbffbf(e,k,sp,p,fll,fl0,&
!$$$         &ng1,ng2,ng3,ng4,ng5,sw, eW, kW,giarray,qiarray,Wid,pol_int) result(res)
!$$$      double complex, intent(in) :: e(:,:), k(:,:), eW(:), kW(:)
!$$$      double complex, intent(in) :: sp(:,:),p(:,:)
!$$$      integer, intent(in) ::  ng1,ng2,ng3,ng4,ng5, sw
!$$$      character, intent(in) :: fll(:)*3
!$$$      character, intent(in) :: fl0*3   ! flavor off-shell f-line
!$$$      integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int, Wid 
!$$$      ! -----------------------------------------------------------------------
!$$$      character :: fl1*3,fl2*3,fl3*3, fl4*3, fl5*3, flaux*3, flaux0*3
!$$$  !    integer             :: m1,m2,m3, ms1a, ms2a!,m
!$$$      integer :: ngluon, ng6,m1,i!, ngL
!$$$      integer, parameter :: Ndumm=0
!$$$      double complex             :: res(size(sp,dim=1))
!$$$      double complex             :: tmp(size(sp,dim=1))
!$$$      double complex             :: k1(size(k,dim=1))
!$$$      double complex             :: k2(size(k,dim=1))
!$$$      double complex             :: k3(size(k,dim=1))
!$$$      double complex             :: k4(size(k,dim=1)), k5(size(k,dim=1))
!$$$      double complex             :: sp1(size(sp,dim=1))
!$$$      double complex             :: sp2(size(sp,dim=1))
!$$$      double complex             :: sp3(size(sp,dim=1))
!$$$      double complex             :: sp4(size(sp,dim=1))
!$$$      double complex             :: sp5(size(sp,dim=1))
!$$$      double complex             :: e1(size(e,dim=1))
!$$$      double complex             :: e2(size(e,dim=1))
!$$$      double complex             :: e3(size(e,dim=1))
!$$$      double complex             :: e4(size(e,dim=1))
!$$$      double complex             :: e5(size(e,dim=1))
!$$$      double complex             :: kdumm(size(k,dim=1),Ndumm)
!$$$      double complex             :: edumm(size(e,dim=1),Ndumm)
!$$$      double complex  :: k1sq,k2sq,k3sq,k4sq, k5sq
!$$$      !real(dp) :: mass
!$$$      logical                   :: done 
!$$$  
!$$$      !if (verbose) write(*,*) 'entering fW_bffbffbf'
!$$$  
!$$$      done = .false. 
!$$$      if (present(giarray)) then 
!$$$         !if (size(qiarray) /= 5) stop 'fW_bffbffbf: wrong size qiarray'
!$$$         !if (size(e,dim=2) /= size(giarray)) stop 'fW_bffbffbf: ng= size(giarray)' 
!$$$         ! XXX 
!$$$         call memory_check(pol_int,res,done,giarray,qiarray)
!$$$         if (done) return 
!$$$      else
!$$$         write(*,*) 'giarray missing'
!$$$      endif
!$$$      res = czero
!$$$      !mass = mt
!$$$  
!$$$      ngluon = size(e,dim=2)
!$$$      ng6 = ngluon - ng1 - ng2-ng3 - ng4 - ng5
!$$$  
!$$$      if (ng6 < 0) write(*,*) 'ERROR IN CURRENT fW_bffbffbf'
!$$$  
!$$$      fl1 = fll(1)
!$$$      fl2 = fll(2)
!$$$      fl3 = fll(3)
!$$$      fl4 = fll(4)
!$$$      fl5 = fll(5)
!$$$  
!$$$   if (fl1 /= fl2) then 
!$$$  
!$$$  
!$$$         if (ngluon == 0) then 
!$$$  
!$$$            res = czero
!$$$            if (fl0 .eq. 'top') flaux = 'bot'
!$$$            if (fl0 .eq. 'bot') flaux = 'top'
!$$$  
!$$$            if (sw ==1) then
!$$$               
!$$$               sp1 = f_bffbffbf(e,k,sp,p,fll,flaux,ng1,ng2,ng3,ng4,ng5,&
!$$$                    &giarray,qiarray, pol_int)
!$$$  
!$$$               k1 = p(:,1)+p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$               k1sq = sc(k1,k1)
!$$$               sp1 = spb2(sp1,k1)!+mass*sp1
!$$$               
!$$$               tmp = vbqW(sp1,eW)
!$$$               
!$$$               if (abs(k1sq) > propcut) then
!$$$                  tmp = ci*tmp/k1sq
!$$$               else
!$$$                  tmp = czero
!$$$               endif
!$$$  
!$$$               res = res + tmp                !#1
!$$$  
!$$$               sp2 = fW(e,k,sp(:,1),p(:,1),fl0,fl1,eW,kW, 0, giarray,&
!$$$                    &qiarray(1:1),Wid,pol_int)
!$$$               
!$$$               k2 = p(:,1)+kW
!$$$               k2sq = sc(k2,k2)
!$$$               sp2 = spb2(sp2,k2)!+mass*sp2
!$$$  
!$$$               e2 = g_sbsfbf(edumm,kdumm,sp(:,2:5),p(:,2:5),fll(2:5),ng1,&
!$$$                    &ng2,ng3,ng4,giarray,qiarray(2:5),pol_int)
!$$$  
!$$$               k3 = p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$               k3sq = sc(k3,k3)
!$$$               
!$$$               tmp = vqg(sp2,e2)
!$$$  
!$$$               if (abs(k2sq) > propcut) then
!$$$                  tmp = ci*tmp/k2sq
!$$$               else
!$$$                  tmp = czero
!$$$               endif
!$$$  
!$$$               if (abs(k3sq) > propcut) then
!$$$                  tmp = -ci*tmp/k3sq
!$$$               else
!$$$                  tmp = czero
!$$$               endif
!$$$               
!$$$               res = res +  tmp                !#2
!$$$  
!$$$              sp4 = fW_bffbf(e,k,sp(:,1:3),p(:,1:3),fll(1:3), fl0, eW,kW,&
!$$$                   &ng1,ng2,ng3,1,giarray,qiarray(1:3),Wid, pol_int)
!$$$  
!$$$              k4 = p(:,1)+p(:,2)+p(:,3)+kW
!$$$              k4sq = sc(k4,k4)
!$$$              sp4 = spb2(sp4, k4)!+mass*sp4
!$$$  
!$$$              e4 = g_fbf(edumm,kdumm,sp(:,4),p(:,4),fl4,sp(:,5),p(:,5),fl5,&
!$$$                   &ng1,ng2,giarray,qiarray(4:5), pol_int)
!$$$              
!$$$              k5 = p(:,4)+p(:,5)
!$$$              k5sq =sc(k5,k5)
!$$$  
!$$$              tmp = vqg(sp4,e4)
!$$$  
!$$$              if (abs(k4sq) > propcut) then
!$$$                 tmp = ci*tmp/k4sq
!$$$              else
!$$$                 tmp  =czero
!$$$              endif
!$$$              
!$$$              if (abs(k5sq) > propcut) then
!$$$                 tmp = -ci*tmp/k5sq
!$$$              else
!$$$                 tmp = czero
!$$$              endif
!$$$  
!$$$             res = res + tmp               !#3
!$$$        
!$$$  !---------------------------
!$$$  
!$$$           elseif (sw == 2) then
!$$$              
!$$$              e1 = gW_sbsfbf(edumm,kdumm,sp(:,2:5),p(:,2:5), fll(2:5), eW,kW,&
!$$$                   &ng1,ng2,ng3,ng4,2,giarray,qiarray(2:5),Wid,pol_int)
!$$$              
!$$$              k1 = p(:,2)+p(:,3)+p(:,4)+p(:,5)+kW
!$$$              k1sq = sc(k1,k1)
!$$$  
!$$$              tmp = vqg(sp(:,1),e1)
!$$$  
!$$$              if (abs(k1sq) > propcut) then
!$$$                 tmp = -ci*tmp/k1sq
!$$$              else
!$$$                 tmp = czero
!$$$              endif
!$$$              
!$$$              res = res+ tmp                 !#1
!$$$  
!$$$              sp2 = fW_bffbf(e,k,sp(:,1:3),p(:,1:3),fll(1:3),fl0,eW,kW,&
!$$$                   &ng1,ng2,ng3,3,giarray,qiarray(1:3),Wid, pol_int)
!$$$  
!$$$              k2 = p(:,1)+p(:,2)+p(:,3)+kW
!$$$              k2sq = sc(k2,k2)
!$$$              sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$  
!$$$              e2 = g_fbf(edumm,kdumm,sp(:,4),p(:,4),fl4,sp(:,5),p(:,5),fl5,&
!$$$                   &ng1,ng2,giarray,qiarray(4:5),pol_int)
!$$$              k3 = p(:,4)+p(:,5)
!$$$              k3sq = sc(k3,k3)
!$$$              
!$$$              tmp = vqg(sp2,e2)
!$$$  
!$$$              if (abs(k2sq) > propcut) then
!$$$                 tmp = ci*tmp/k2sq
!$$$              else
!$$$                 tmp = czero
!$$$              endif
!$$$  
!$$$              if (abs(k3sq) > propcut) then
!$$$                 tmp = -ci*tmp/k3sq
!$$$              else
!$$$                 tmp = czero
!$$$              endif
!$$$  
!$$$              res = res + tmp                       !#2
!$$$  
!$$$  !--------------------------            
!$$$  
!$$$           elseif (sw == 3) then
!$$$              
!$$$              e1 = gW_sbsfbf(edumm, kdumm, sp(:,2:5),p(:,2:5), fll(2:5),eW,kW,&
!$$$                   &ng1,ng2,ng3,ng4,4,giarray,qiarray(2:5),Wid,pol_int)
!$$$      
!$$$  
!$$$              k1 = p(:,2)+p(:,3)+p(:,4)+p(:,5)+kW
!$$$              k1sq = sc(k1,k1)
!$$$  
!$$$              tmp = vqg(sp(:,1),e1)
!$$$  
!$$$              if (abs(k1sq) > propcut) then
!$$$                 tmp = -ci*tmp/k1sq
!$$$              else
!$$$                 tmp = czero
!$$$              endif
!$$$  
!$$$            res = res + tmp                !#1
!$$$  
!$$$              sp2 = f_bffbf_2(e,k,sp(:,1:3),p(:,1:3),fll(1:3),fl0, ng1,ng2, ng3,&
!$$$                   &giarray,qiarray(1:3),pol_int)
!$$$              k2 = p(:,1)+p(:,2)+p(:,3)
!$$$              k2sq = sc(k2,k2)
!$$$              sp2 = spb2(sp2,k2)!+mass*sp2
!$$$              
!$$$              e2 = gW_fbf(edumm,kdumm, sp(:,4), p(:,4),fl4, sp(:,5),p(:,5), fl5,&
!$$$                   &eW,kW,ng1, ng2, giarray, qiarray(4:5), Wid, pol_int)
!$$$    
!$$$              k3 = p(:,4)+p(:,5)+kW
!$$$              k3sq = sc(k3,k3)
!$$$  
!$$$              tmp = vqg(sp2, e2)
!$$$  
!$$$              if (abs(k2sq) > propcut) then
!$$$                 tmp = ci*tmp/k2sq
!$$$              else
!$$$                 tmp = czero
!$$$              endif
!$$$              
!$$$              if (abs(k3sq) >  propcut) then
!$$$                 tmp = -ci*tmp/k3sq
!$$$              else
!$$$                 tmp = czero
!$$$              endif
!$$$  
!$$$              res = res+tmp              !#2
!$$$  
!$$$  !-------------------------------
!$$$           elseif (sw == 4) then
!$$$  
!$$$              sp1 = fW_bffbf(e,k,sp(:,1:3),p(:,1:3),fll(1:3), fl0,eW,kW,&
!$$$                   ng1,ng2,ng3,2, giarray,qiarray(1:3),Wid,pol_int)
!$$$  
!$$$              k1 = p(:,1)+p(:,2)+p(:,3)+kW
!$$$              k1sq = sc(k1,k1)
!$$$              sp1 = spb2(sp1,k1)!+mass*sp1
!$$$  
!$$$              e1 = g_fbf(edumm,kdumm,sp(:,4),p(:,4),fl4,sp(:,5),p(:,5),fl5,&
!$$$                   &ng1,ng2,giarray,qiarray(4:5),pol_int)
!$$$  
!$$$              k2 = p(:,4)+p(:,5)
!$$$              k2sq = sc(k2,k2)
!$$$  
!$$$              tmp = vqg(sp1, e1)
!$$$  
!$$$              if (abs(k1sq) > propcut) then
!$$$                 tmp = ci*tmp/k1sq
!$$$              else
!$$$                 tmp = czero
!$$$              endif
!$$$              
!$$$              if (abs(k2sq) >  propcut) then
!$$$                 tmp = -ci*tmp/k2sq
!$$$              else
!$$$                 tmp = czero
!$$$              endif
!$$$      
!$$$              res = res+tmp          !#1
!$$$  
!$$$              
!$$$              sp2 = f_bffbf(e,k,sp(:,3:5),p(:,3:5),fll(3:5),fl0,ng1,ng2,ng3,&
!$$$                   &giarray,qiarray(3:5),pol_int)
!$$$           
!$$$              k2 = p(:,3)+p(:,4)+p(:,5)
!$$$              k2sq = sc(k2,k2)
!$$$              sp2 = spb2(sp2,k2)!+mass*sp2
!$$$  
!$$$              e2 = gW_bff(edumm,kdumm,sp(:,1),p(:,1),fl1,sp(:,2),p(:,2),fl2,&
!$$$                   &eW,kW,ng1, ng2,giarray,qiarray(1:2),Wid,pol_int)
!$$$  
!$$$  
!$$$              k3 = p(:,1)+p(:,2)+kW
!$$$              k3sq = sc(k3,k3)
!$$$              
!$$$              tmp = vgq(e2, sp2)
!$$$              
!$$$              if (abs(k2sq) > propcut) then
!$$$                 tmp = ci*tmp/k2sq
!$$$              else
!$$$                 tmp = czero
!$$$              endif
!$$$              
!$$$              if (abs(k3sq) >  propcut) then
!$$$                 tmp = -ci*tmp/k3sq
!$$$              else
!$$$                 tmp = czero
!$$$              endif
!$$$      
!$$$              res = res+tmp          !#2
!$$$  
!$$$              
!$$$             
!$$$  
!$$$        ! -----------------------      
!$$$           endif
!$$$        
!$$$  
!$$$        else 
!$$$           write(*,*) 'Function fW_bffbffbf only written for nlguons = 0'
!$$$        endif
!$$$  
!$$$     elseif (fl1 .eq. fl2) then
!$$$      
!$$$           if (ngluon ==0) then
!$$$              if (fl0 == 'top') flaux0 = 'bot'
!$$$              if (fl0 == 'bot') flaux0 = 'top'
!$$$              
!$$$  
!$$$              if (sw == 1) then
!$$$                 
!$$$                 sp1 = f_bffbffbf(e,k,sp,p,fll,flaux0,ng1,ng2,ng3,ng4,ng5,&
!$$$                    &giarray,qiarray, pol_int)
!$$$                 WRITE(*,*) '---'
!$$$                 k1 = p(:,1)+p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$                 k1sq = sc(k1,k1)
!$$$                 sp1 = spb2(sp1,k1)!+mass*sp1
!$$$              
!$$$                 tmp = vbqW(sp1,eW)
!$$$              
!$$$                 if (abs(k1sq) > propcut) then
!$$$                    tmp = ci*tmp/k1sq
!$$$                 else
!$$$                    tmp = czero
!$$$                 endif
!$$$            
!$$$                 res = res + tmp 
!$$$                 
!$$$                 sp2 = fW_bffbf_2(e,k,sp(:,1:3), p(:,1:3),fll(1:3), fl0,eW,kW,&
!$$$                      &ng1,ng2,ng3,1,giarray,qiarray(1:3),Wid,pol_int)
!$$$              
!$$$                 k2 = p(:,1)+p(:,2)+p(:,3)+kW
!$$$                 k2sq =sc(k2,k2)
!$$$                 sp2 = spb2(sp2,k2)!+mass*sp2
!$$$              
!$$$                 e2 = g_fbf(edumm,kdumm,sp(:,4),p(:,4),fl4,sp(:,5),p(:,5),fl5,&
!$$$                      &ng1,ng2,giarray,qiarray(4:5),pol_int)
!$$$                 
!$$$                 k3 = p(:,4)+p(:,5)
!$$$                 k3sq =sc(k3,k3)
!$$$                 
!$$$                 tmp = vqg(sp2,e2)
!$$$                 if (abs(k2sq) > propcut) then 
!$$$                    tmp = ci/k2sq*tmp
!$$$                 else 
!$$$                    tmp = czero 
!$$$                 endif
!$$$              
!$$$                 if (abs(k3sq) > propcut) then 
!$$$                    tmp = -ci/k3sq*tmp 
!$$$                 else 
!$$$                    tmp = czero
!$$$                 endif
!$$$        
!$$$                res = res+tmp 
!$$$             
!$$$              elseif (sw == 3) then
!$$$         
!$$$                 sp1 = fW_bffbf(e,k,sp(:,3:5),p(:,3:5),fll(3:5), fl3, eW,kW, &
!$$$                      &ng1, ng2,ng3,3,giarray, qiarray(3:5), Wid,pol_int)
!$$$              
!$$$                 k1 = p(:,3)+p(:,4)+p(:,5)+kW
!$$$                 k1sq = sc(k1,k1)
!$$$                 sp1 = spb2(sp1, k1)!+mass*sp1
!$$$                 
!$$$                 e1 = g_bff(edumm,kdumm,sp(:,1),p(:,1), fl1, sp(:,2),p(:,2), fl2, &
!$$$                      &ng1, ng2, giarray, qiarray(1:2), pol_int)
!$$$                 
!$$$                 k2 = p(:,1)+p(:,2)
!$$$                 
!$$$                 k2sq =sc(k2,k2)
!$$$  
!$$$                 tmp = vgq(e1,sp1)   
!$$$  
!$$$                 if (abs(k2sq) > propcut) then 
!$$$                    tmp = -ci/k2sq*tmp
!$$$                 else 
!$$$                    tmp = czero 
!$$$                 endif
!$$$  
!$$$                 if (abs(k1sq) > propcut) then 
!$$$                    tmp = ci/k1sq*tmp 
!$$$                 else 
!$$$                    tmp = czero
!$$$                 endif
!$$$        
!$$$            res = res+tmp
!$$$  
!$$$            !! RR - not sure about this at all...
!$$$  
!$$$            sp1 = f_bffbf(e,k,sp,p,fll(1:3), fl0, ng1, ng2,ng3,&
!$$$                 &giarray, qiarray(1:3), pol_int)
!$$$            k1 = p(:,3)+p(:,2)+p(:,1)
!$$$            k1sq = sc(k1,k1)
!$$$            sp1 = spb2(sp1, k1)!+mass*sp1
!$$$            
!$$$            e1 = gW_fbf(edumm,kdumm,sp(:,4),p(:,4), fl4, sp(:,5),p(:,5), fl5, &
!$$$                 &eW,kW,ng1, ng2, giarray, qiarray(4:5),Wid, pol_int)
!$$$  
!$$$            k2 = p(:,4)+p(:,5)+kW
!$$$  
!$$$            k2sq =sc(k2,k2)
!$$$  
!$$$            
!$$$            tmp  = vqg(sp1,e1)
!$$$            if (abs(k2sq) > propcut) then 
!$$$               tmp = -ci/k2sq*tmp
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            if (abs(k1sq) > propcut) then 
!$$$               tmp = ci/k1sq*tmp 
!$$$            else 
!$$$               tmp = czero
!$$$            endif
!$$$        
!$$$            res = res+tmp
!$$$         elseif (sw ==2) then
!$$$            
!$$$            if (fl1 == 'str') then
!$$$               sp1 = fW_bffbf(e,k,sp(:,3:5),p(:,3:5),fll(3:5),fl0,eW,kW,&
!$$$                    &ng1,ng2,ng3,1,giarray,qiarray(3:5),Wid,pol_int)
!$$$        
!$$$               k1 = p(:,3)+p(:,4)+p(:,5)+kW
!$$$               k1sq = sc(k1,k1)
!$$$               sp1 = spb2(sp1,k1)!+mass*sp1
!$$$               
!$$$               e1 = g_bff(edumm,kdumm,sp(:,1),p(:,1),fl1,sp(:,2),p(:,2),fl2,&
!$$$                    &ng1,ng2,giarray,qiarray(1:2),pol_int)
!$$$   
!$$$               k2 = p(:,1)+p(:,2)
!$$$               k2sq = sc(k2,k2)
!$$$               
!$$$               tmp = vgq(e1,sp1)
!$$$            
!$$$               if (abs(k1sq)>propcut) then
!$$$                  tmp = ci*tmp/k1sq
!$$$               else
!$$$                  tmp = czero
!$$$               endif
!$$$            
!$$$               if (abs(k2sq) > propcut) then 
!$$$                  tmp = -ci/k2sq*tmp
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$            
!$$$              
!$$$               res = res+tmp
!$$$  
!$$$          
!$$$               sp4 = fW_bffbf_2(e,k,sp(:,1:3),p(:,1:3),fll(1:3),fl0,eW,kW,&
!$$$                    &ng1,ng2,ng3,3,giarray,qiarray(1:3),Wid, pol_int)
!$$$               
!$$$               k4 = p(:,1)+p(:,2)+p(:,3)+kW
!$$$               k4sq = sc(k4,k4)
!$$$               sp4 = spb2(sp4,k4)!+mass*sp4
!$$$               
!$$$               e3 = g_fbf(edumm,kdumm,sp(:,4),p(:,4),fl4,sp(:,5),p(:,5),fl5,&
!$$$                    &ng1,ng2,giarray,qiarray(4:5),pol_int)
!$$$               
!$$$               k3 = p(:,4)+p(:,5)
!$$$               k3sq = sc(k3,k3)
!$$$               
!$$$               tmp = vqg(sp4,e3)
!$$$            
!$$$             
!$$$            
!$$$               if (abs(k4sq)>propcut) then
!$$$                  tmp = ci*tmp/k4sq
!$$$               else
!$$$                  tmp = czero
!$$$               endif
!$$$            
!$$$               if (abs(k3sq) > propcut) then 
!$$$                  tmp = -ci/k3sq*tmp
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$  
!$$$               res = res+tmp
!$$$  
!$$$               
!$$$            elseif (fl1 /= 'str') then
!$$$  
!$$$               sp5 = f(e,k,sp(:,1),p(:,1),fl1,fl0,ng1,giarray,qiarray(1:1),pol_int)
!$$$               
!$$$               e5 = gW_sbsfbf_3(edumm,kdumm, sp(:,2:5),p(:,2:5),fll(2:5),eW,kW,&
!$$$                    &ng1,ng2,ng3,ng4,2,giarray,qiarray(2:5),Wid,pol_int)
!$$$               
!$$$               k5 = p(:,2)+p(:,3)+p(:,4)+p(:,5)+kW
!$$$               k5sq = sc(k5,k5)
!$$$               
!$$$               tmp = vqg(sp5,e5)
!$$$               
!$$$               if (abs(k5sq) > propcut) then 
!$$$                  tmp = -ci/k5sq*tmp
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$  
!$$$               res = res+tmp
!$$$            endif
!$$$         
!$$$         else
!$$$            stop 'fw_bffbffbf: sw/=1 and sw /= 3 and sw /= 2'
!$$$         endif
!$$$  
!$$$      else
!$$$         write(*,*) 'Function fW_bffbffbf not written for ngluons /=0'
!$$$      endif
!$$$  
!$$$   endif
!$$$  
!$$$      ! -- store current 
!$$$      if (present(giarray)) call store_result(pol_int,res,giarray,qiarray,Wid)
!$$$  
!$$$  
!$$$    end function fW_bffbffbf

end module recurrenceC
