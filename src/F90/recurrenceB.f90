  !------------------------------------------------------------------------!
  ! Authors: Giulia Zanderighi...	                                   !
  !------------------------------------------------------------------------!
!! File generated automatically by autogen.pl from 
!! general precision template PRECfiles/genPREC/PRECrecurrenceGbitstwo.f90.

module recurrenceB
  use consts_dp
  use spinfns
  use recurrenceA
  implicit none

  public :: g_bff,g_fbf,f_bffbf
  public :: g_sbsfbf, bf_fbff, bf_fbff_2
  public :: f_bffbf_2, f_bffbf_3
!$$  public :: f_fbfbf_3, f_bffbffbf

  logical :: verbose = .false. 
  
  private
  
  contains 

  ! ----- Gluon current, with additional qbq  pair
  ! ----- first fermion is outgoing line, second (bar-fermion) is 
  ! ----- incoming fermion line
  ! ------ ms1 gives the number of gluon lines to the left of 
  ! ------ outgoing fermion line; ms2 gives the number of fermion 
  ! ------ lines between the two fermion lines
  recursive function g_bff(e,k,sp1,p1,fl1,sp2,p2,fl2,ms1,ms2,&
       &giarray,qiarray,pol_int) result(res)
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp1(:), p1(:), sp2(:), p2(:)
    integer, intent(in) ::  ms1,ms2
    character, intent(in) :: fl1*3, fl2*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int 
    ! -----------------------------------------------------------------------
    integer             :: m1,m2, ms1a, ms2a!,m,m3
    integer :: ngluon, ng1, ng2, ng3, ngL,Dv
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
!    real(dp) :: mass,mass2
    logical                   :: done 

    if (fl1.ne.fl2) then 
       res = czero 
       return 
    endif

    if (verbose) write(*,*) 'entering g_bff',size(e,dim=2),ms1,ms2
    if (verbose) write(*,*) 'entering g_bff:qi',qiarray
    if (verbose .and. present(giarray)) write(*,*) 'entering g_bff:gi',giarray
    done = .false. 
!    if (present(giarray)) then 
!       if (size(qiarray) /= 2) stop 'g_bff: wrong size qiarray'
!       if (size(e,dim=2) /= size(giarray)) stop 'g_bff: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray)
!       if (done) return 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'g_bff: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif


!    mass = mt
!    mass2 = mt**2

    if (size(sp1) == 4)  Dv=4
    if (size(sp1) == 8)  Dv=6
    if (size(sp1) == 16)  Dv=8

    ngluon = size(e,dim=2)
    ng1 = ms1
    ng2 = ms2
    ng3 = ngluon - ms1 - ms2


    if (ng3 < 0) write(*,*) 'ERROR IN CURRENT A:g_bff'

    if (ngluon == 0) then 

       res = vbqq(Dv,sp1,sp2)

    else

       res = czero


       do m1=1,ng1

          e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)
          k1=sum(k(:,1:m1),dim=2)
          k1sq = sc(k1,k1)

          ms1a = ng1-m1
          e2=g_bff(e(:,m1+1:ngluon),k(:,m1+1:ngluon),&
               &sp1,p1,fl1,sp2,p2,fl2,ms1a,ng2,giarray(m1+1:ngluon),qiarray,&
               &pol_int)
          if (m1+1<=ngluon) then 
             k2 = sum(k(:,m1+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2 = k2 + p1 + p2
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

          res = res + tmp  !#1


          do m2 = 0,ng3-1    

             ms1a=ng1-m1
             e2=g_bff(e(:,m1+1:ng1+ng2+m2),k(:,m1+1:ng1+ng2+m2),&
                  &sp1,p1,fl1,sp2,p2,fl2,ms1a,ng2,&
                  &giarray(m1+1:ng1+ng2+m2),qiarray(1:2),pol_int)
             if (m1+1<=ng1+ng2+m2) then 
                k2=sum(k(:,m1+1:ng1+ng2+m2),dim=2)
             else
                k2 = czero 
             endif
             k2 = k2 + p1 + p2
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

             res = res + tmp  !#4

          enddo


          if (m1.lt.ng1) then 

             do m2 = m1+1,ng1   

                e2=vgluon(e(:,m1+1:m2),k(:,m1+1:m2),&
                     &giarray(m1+1:m2),pol_int)
                if (m1+1<=m2) then 
                   k2=sum(k(:,m1+1:m2),dim=2)
                else
                   k2 = czero 
                endif
                k2sq=sc(k2,k2)


                ms1a=ng1-m2
                e3=g_bff(e(:,m2+1:ngluon),k(:,m2+1:ngluon),&
                     &sp1,p1,fl1,sp2,p2,fl2,ms1a,ng2,&
                     &giarray(m2+1:ngluon),qiarray(1:2),pol_int)
                if (m2+1<=ngluon) then 
                   k3=sum(k(:,m2+1:ngluon),dim=2)
                else
                   k3 = czero 
                endif 
                k3 = k3 + p1 + p2
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

                res = res + tmp  !#3

             enddo

          endif

       enddo                          !#1 & #4 & #3

       do m1=0,ng3-1 

          e1 = g_bff(e(:,1:ng1+ng2+m1),k(:,1:ng1+ng2+m1),&
               &sp1,p1,fl1,sp2,p2,fl2,ng1,ng2,&
               &giarray(1:ng1+ng2+m1),qiarray(1:2),pol_int)
          if (1<=ng1+ng2+m1) then 
             k1=sum(k(:,1:ng1+ng2+m1),dim=2)
          else
             k1 = czero 
          endif
          k1 = k1 + p1 + p2
          k1sq=sc(k1,k1)

          e2=vgluon(e(:,ng1+ng2+m1+1:ngluon),k(:,ng1+ng2+m1+1:ngluon),&
               &giarray(ng1+ng2+m1+1:ngluon),pol_int)
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

          res = res + tmp  !#2


          if (m1.lt.ng3-2) then 

             ngL=ng1+ng2+m1

             do m2=ngL+1,ngluon-1

                e2=vgluon(e(:,ngL+1:m2),k(:,ngL+1:m2),&
                     &giarray(ngL+1:m2),pol_int)
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

                res = res + tmp  !#5

             enddo

          endif

       enddo   !#2 & #5




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
          k3sq = sc(k3,k3)!-mass2

          if (ng1 > 0.or.m1 > 0) sp3 = spb2(sp3,k3)!+mass*sp3

          ms2a=ng2-m1
!          write(*,*) 'calling bf from g_bff' 
          sp4=bf(e(:,ms1a+1:ngluon),k(:,ms1a+1:ngluon),&
               &sp2,p2,fl2,fl2,ms2a,&
               &giarray(ms1a+1:ngluon),qiarray(2:2),pol_int)
          if (ms1a+1<=ngluon) then 
             k4=sum(k(:,ms1a+1:ngluon),dim=2)
          else
             k4 = czero 
          endif
          k4 = - k4 - p2
          k4sq = sc(k4,k4)!-mass2

          if (ng3 > 0.or.ng2-m1>0) sp4 = spi2(k4,sp4)!+mass*sp4

          tmp = vbqq(Dv,sp3,sp4)


          if (ng1 > 0.or.m1 > 0) then  
             if (abs(k3sq) > propcut) then 
                tmp= ci/k3sq*tmp
             else
                tmp = czero 
             endif
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


    endif   !end condition for ngluon

    ! -- store current 
!    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray)


  end function g_bff



  recursive function g_fbf(e,k,sp1,p1,fl1,sp2,p2,fl2,ng1,ng2,&
       &giarray,qiarray,pol_int) result(res)
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp1(:), p1(:), sp2(:), p2(:)
    integer, intent(in)       ::  ng1,ng2
    character, intent(in)     :: fl1*3, fl2*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int 
    ! -----------------------------------------------------------------------
    integer                   :: m1,m2, ms1a, ms2a,Dv!,m,m3
    integer                   :: ngluon,  ng3, ngL
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
    double complex             :: k1sq,k2sq,k3sq,k4sq
!    real(dp)                :: mass,mass2
    logical                   :: done 

    if (fl1.ne.fl2) then 
       res = czero 
       return 
    endif

    if (size(sp1) == 4)  Dv=4
    if (size(sp1) == 8)  Dv=6
    if (size(sp1) == 16)  Dv=8

    !if (size(e,dim=2) == 0) then 
    !   res = -cone*vbqq(Dv,sp2,sp1)
    !   return 
    !endif

    if (verbose) write(*,*) 'entering g_fbf',size(e,dim=2),ng1,ng2
    done = .false. 
!    if (present(giarray)) then 
!       if (size(qiarray) /= 2) stop 'g_fbf: wrong size qiarray'
!       if (size(e,dim=2) /= size(giarray)) stop 'g_fbf: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray)
!       if (done) return 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'g_fbf: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif


!    mass = mt
!    mass2 = mt**2



    ngluon = size(e,dim=2)
    ng3 = ngluon - ng1 - ng2

    if (verbose) write(*,*) 'in function g_fbf', ng1, ngluon

    if (ng3 < 0) write(*,*) 'ERROR IN CURRENT B:g_fbf'


    if (ngluon == 0) then 

       res = -cone*vbqq(Dv,sp2,sp1)
    else

       res = czero

       do m1=1,ng1

          e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)
          k1=sum(k(:,1:m1),dim=2) 
          k1sq = sc(k1,k1)

          ms1a = ng1-m1
          e2=g_fbf(e(:,m1+1:ngluon),k(:,m1+1:ngluon),&
               &sp1,p1,fl1,sp2,p2,fl2,ms1a,ng2,giarray(m1+1:ngluon),qiarray,&
               &pol_int)
          if (m1+1<=ngluon) then 
             k2 = sum(k(:,m1+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2 = k2 + p1 + p2
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

          e1 = g_fbf(e(:,1:ng1+ng2+m1),k(:,1:ng1+ng2+m1),&
               &sp1,p1,fl1,sp2,p2,fl2,ng1,ng2,giarray(1:ng1+ng2+m1),qiarray,&
               &pol_int)
          if (1<=ng1+ng2+m1) then 
             k1=sum(k(:,1:ng1+ng2+m1),dim=2)
          else
             k1 = czero 
          endif
          k1 = k1 + p1 + p2
          k1sq=sc(k1,k1)

          e2=vgluon(e(:,ng1+ng2+m1+1:ngluon),k(:,ng1+ng2+m1+1:ngluon),&
               &giarray(ng1+ng2+m1+1:ngluon),pol_int)
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

             e2=vgluon(e(:,m1+1:m2),k(:,m1+1:m2),&
                  &giarray(m1+1:m2),pol_int)
             if (m1+1<=m2) then 
                k2=sum(k(:,m1+1:m2),dim=2)
             else
                k2 = czero
             endif
             k2sq=sc(k2,k2)


             ms1a=ng1-m2
             e3=g_fbf(e(:,m2+1:ngluon),k(:,m2+1:ngluon),&
                  &sp1,p1,fl1,sp2,p2,fl2,ms1a,ng2,&
                  &giarray(m2+1:ngluon),qiarray,pol_int)
             if (m2+1<=ngluon) then 
                k3=sum(k(:,m2+1:ngluon),dim=2)
             else
                k3 = czero 
             endif
             k3 = k3 + p1 + p2
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
             e2=g_fbf(e(:,m1+1:ng1+ng2+m2),k(:,m1+1:ng1+ng2+m2),&
                  &sp1,p1,fl1,sp2,p2,fl2,ms1a,ng2,&
                  &giarray(m1+1:ng1+ng2+m2),qiarray,pol_int)
             if (m1+1<=ng1+ng2+m2) then 
                k2=sum(k(:,m1+1:ng1+ng2+m2),dim=2)
             else
                k2 = czero 
             endif
             k2 = k2 + p1 + p2
             k2sq=sc(k2,k2)


             e3=vgluon(e(:,ng1+ng2+m2+1:ngluon),k(:,ng1+ng2+m2+1:ngluon),&
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

          e1=g_fbf(e(:,1:ngL),k(:,1:ngL),&
               &sp1,p1,fl1,sp2,p2,fl2,ng1,ng2,&
               &giarray(1:ngL),qiarray,pol_int)
          if (1<=ngL) then 
             k1=sum(k(:,1:ngL),dim=2)
          else
             k1 = czero 
          endif
          k1 = k1 + p1 + p2
          k1sq=sc(k1,k1)

          do m2=ngL+1,ngluon-1

             e2=vgluon(e(:,ngL+1:m2),k(:,ngL+1:m2),&
                  &giarray(ngL+1:m2),pol_int)
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



          if (ng1 > 0.or.m1 > 0) then 
             if (abs(k3sq) > propcut) then 
                tmp= ci/k3sq*tmp
             else 
                tmp = czero 
             endif
          endif

          if (ng3 > 0.or.ng2-m1>0) then 
             if (abs(k4sq) > propcut) then 
                tmp= ci/k4sq*tmp 
             else 
                tmp = czero 
             endif
          endif
!          write(*,*) 'tmp',(tmp)
          res = res + tmp

       enddo  !#6

    endif   !end condition for ngluon 

    ! -- store current 
!    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray)
    


  end function g_fbf


  recursive function f_bffbf(e,k,sp,p,fll,fl0,ng1,ng2,ng3,&
       &giarray,qiarray,pol_int) result(res)
    implicit none
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp(:,:),p(:,:)
    integer, intent(in) ::  ng1,ng2,ng3
    character, intent(in) :: fll(:)*3
    character, intent(in) :: fl0*3   ! flavor off-shell f-line
    integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int 
    ! -----------------------------------------------------------------------
    character :: fl1*3,fl2*3,fl3*3
!    integer             :: m1,m2,m3, ms1a, ms2a,m
    integer :: ngluon, ng4, ngL,m
    integer, parameter :: Ndumm=0
    double complex             :: res(size(sp,dim=1))
    double complex             :: res_stored(size(sp,dim=1))
    double complex             :: tmp(size(sp,dim=1))
    double complex             :: k1(size(k,dim=1))
    double complex             :: k2(size(k,dim=1))
!    double complex             :: k3(size(k,dim=1))
!    double complex             :: k4(size(k,dim=1))
    double complex             :: sp1(size(sp,dim=1))
    double complex             :: sp2(size(sp,dim=1))
!    double complex             :: sp3(size(sp,dim=1))
!    double complex             :: sp4(size(sp,dim=1))
    double complex             :: e1(size(e,dim=1))
    double complex             :: e2(size(e,dim=1))
!    double complex             :: e3(size(e,dim=1))
    double complex             :: kdumm(size(k,dim=1),Ndumm)
    double complex             :: edumm(size(e,dim=1),Ndumm)
    double complex  :: k1sq,k2sq!,k3sq!,k4sq
!    real(dp) :: mass,mass2
    logical                   :: done, onshell  

    if (verbose) write(*,*) 'entering f_bffbf',ng1,ng2,ng3
    done = .false. 
!    if (present(giarray)) then 
!       if (size(qiarray) /= 3) stop 'f_bffbf: wrong size qiarray'
!       if (size(e,dim=2) /= size(giarray)) stop 'bf_fbf: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray)
!       if (done) return ! XXX 
!!       if (done) res_stored = res 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'b_bffbf: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif

!    mass = mt
!    mass2 = mt**2

    ngluon = size(e,dim=2)
    ng4 = ngluon - ng1 - ng2-ng3

    !if (verbose) write(*,*) 'in function f_bffbf', ng1, ngluon

    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)


    if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C:f_bffbf'

!$$$ TM FOR NOW
!$$$      if ( abs(pmass(p(:,1)+p(:,2)+sum(k,dim=2)) - mz) < tol) then 
!$$$         onshell =.true. 
!$$$      else
!$$$         onshell = .false.
!$$$      endif

    if (ngluon == 0) then 

       res = czero

       e2 = g_fbf(edumm,kdumm,sp(:,2),p(:,2),fl2,sp(:,3),p(:,3),fl3,0,0,&
            &giarray,qiarray(2:3),pol_int)

       k1 = p(:,2)+p(:,3)
       k1sq=sc(k1,k1)

       if (abs(k1sq) > propcut.and.fl0==fl1) then 
          tmp = -ci/k1sq*vqg(sp(:,1),e2)
       else 
          tmp = czero 
       endif

       res = res + tmp 

       e2 = g_bff(edumm,kdumm,sp(:,1),p(:,1),fl1,sp(:,2),p(:,2),fl2,0,0,&
            &giarray,qiarray(1:2),pol_int)
       k2 = p(:,1)+p(:,2) 
       k2sq = sc(k2,k2)



       if (abs(k2sq) > propcut.and.fl0==fl3) then 
          tmp = -ci/k2sq*vgq(e2,sp(:,3))
       else 
          tmp = czero 
       endif

       res = res + tmp

    else                          !for #gluons > 0

       res = czero

       do m=0,ng2

          ! need to comment this out to get extra_fermion1 to work...        
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

          if (ng1 > 0.or.m>0) then 
             if (abs(k1sq) > propcut) then 
                tmp = ci/k1sq*tmp
             else 
                tmp = czero 
             endif
          endif


          res = res + tmp

       enddo  !#1

       do m=0,ng4-1

          ngL = ng1+ ng2+ng3+m      

          sp1=f_bffbf(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0,&
               &ng1,ng2,ng3,&
               &giarray(1:ngL),qiarray,pol_int)
          if (1<=ngL) then 
             k1=sum(k(:,1:ngL),dim=2)
          else
             k1 = czero
          endif 
          k1 = k1 + p(:,1)+p(:,2)+p(:,3)
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


       do m=1,ng1

          e1 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)
          k1=sum(k(:,1:m),dim=2)
          k1sq=sc(k1,k1)

          sp2=f_bffbf(e(:,m+1:ngluon),k(:,m+1:ngluon),&
               &sp,p,fll,fl0,ng1-m,ng2,ng3,&
               &giarray(m+1:ngluon),qiarray,pol_int)
          if (m+1<=ngluon) then 
             k2=sum(k(:,m+1:ngluon),dim=2)
          else
             k2 = czero 
          endif 
          k2 = k2 + p(:,1) + p(:,2) + p(:,3)
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



       do m=0,ng3

          ngL = ng1+ ng2+m      

          e1 = g_bff(e(:,1:ngL),k(:,1:ngL),sp(:,1),p(:,1),fl1,&
               &sp(:,2),p(:,2),fl2,ng1,ng2,giarray(1:ngL),qiarray(1:2),pol_int)
          if (1<=ngL) then 
             k1=sum(k(:,1:ngL),dim=2)+p(:,1)+p(:,2)
          else
             k1 = p(:,1)+p(:,2)
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

          res = res + tmp

       enddo  !#4



    endif

    ! -- store current 
!    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray)

    if (verbose) write(*,*) 'done f_bffbf ',done, ' ', fl0,fl1,fl2,fl3,ng1,ng2,ng3,pol_int, 'g',giarray,'q',qiarray, res

!    if (done .and. any(res-res_stored /= czero)) then 
!       write(*,*) ':stored',res
!       write(*,*) ':stored',res_stored
!       write(*,*) ':stored',res-res_stored
!    endif


  end function f_bffbf


  recursive function f_bffbf_2(e,k,sp,p,fll,fl0,ng1,ng2,ng3,&
       &giarray,qiarray,pol_int) result(res)
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp(:,:),p(:,:)
    integer, intent(in) ::  ng1,ng2,ng3
    character, intent(in) :: fll(:)*3
    character, intent(in) :: fl0*3   ! flavor off-shell f-line
    integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int 
    ! -----------------------------------------------------------------------
    character :: fl1*3,fl2*3,fl3*3
!    integer             :: m1,m2,m3, ms1a, ms2a
    integer :: ngluon, ng4, ngL,m
    integer, parameter :: Ndumm=0
    double complex             :: res(size(sp,dim=1))
    double complex             :: res_stored(size(sp,dim=1))
    double complex             :: tmp(size(sp,dim=1))
    double complex             :: k1(size(k,dim=1))
    double complex             :: k2(size(k,dim=1))
    double complex             :: sp1(size(sp,dim=1))
    double complex             :: sp2(size(sp,dim=1))
    double complex             :: e1(size(e,dim=1))
    double complex             :: e2(size(e,dim=1))
    double complex             :: kdumm(size(k,dim=1),Ndumm)
    double complex             :: edumm(size(e,dim=1),Ndumm)
    double complex  :: k1sq,k2sq!,k3sq!,k4sq
!    real(dp) :: mass,mass2
    logical                   :: done 

    if (verbose) write(*,*) 'entering f_bffbf_2',ng1,ng2,ng3
    done = .false. 
!    if (present(giarray)) then 
!       if (size(qiarray) /= 3) stop 'f_bffbf: wrong size qiarray'
!       if (size(e,dim=2) /= size(giarray)) stop 'bf_fbf: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray)
!       if (done) return 
!       !res_stored = res 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'f_bffbf_2: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif

!    mass = mt
!    mass2 = mt**2

    ngluon = size(e,dim=2)
    ng4 = ngluon - ng1 - ng2-ng3

    !if (verbose) write(*,*) 'in function f_bffbf', ng1, ngluon

    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)


    if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C:f_bffbf_2'

    if (ngluon == 0) then 

       res = czero

       e2 = g_fbf(edumm,kdumm,sp(:,2),p(:,2),fl2,sp(:,3),p(:,3),fl3,0,0,&
            &giarray,qiarray(2:3),pol_int)

       k1 = p(:,2)+p(:,3)
       k1sq=sc(k1,k1)

       if (abs(k1sq) > propcut.and.fl0==fl1) then 
          tmp = -ci/k1sq*vqg(sp(:,1),e2)
       else 
          tmp = czero 
       endif

       res = res + tmp 

    else                          !for #gluons > 0

       res = czero

       do m=0,ng2

          ! ???
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

          if (ng1 > 0.or.m>0) then 
             if (abs(k1sq) > propcut) then 
                tmp = ci/k1sq*tmp
             else 
                tmp = czero 
             endif
          endif


          res = res + tmp

       enddo  !#1

       do m=0,ng4-1

          ngL = ng1+ ng2+ng3+m      

          sp1=f_bffbf_2(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0,&
               &ng1,ng2,ng3,&
               &giarray(1:ngL),qiarray,pol_int)
          if (1<=ngL) then 
             k1=sum(k(:,1:ngL),dim=2)
          else
             k1 = czero 
          endif 
          k1 = k1 + p(:,1)+p(:,2)+p(:,3)
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


       do m=1,ng1

          e1 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)
          k1=sum(k(:,1:m),dim=2)
          k1sq=sc(k1,k1)

          sp2=f_bffbf_2(e(:,m+1:ngluon),k(:,m+1:ngluon),&
               &sp,p,fll,fl0,ng1-m,ng2,ng3,&
               &giarray(m+1:ngluon),qiarray,pol_int)
          if (m+1<=ngluon) then 
             k2=sum(k(:,m+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2 = k2 + p(:,1) + p(:,2) + p(:,3)
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


    endif

    !! -- store current 
!    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray)

    if (verbose) write(*,*) 'done f_bffbf_2',ng1,ng2,ng3,pol_int, 'g',giarray,'q',qiarray, res
    if (verbose) write(*,*) 'done f_bffbf_2 ',done, ' ', fl0,fl1,fl2,fl3,ng1,ng2,ng3,pol_int, 'g',giarray,'q',qiarray, res

!    if (done .and. any(res-res_stored /= czero)) then 
!       write(*,*) ':stored',res
!       write(*,*) ':stored',res_stored
!       write(*,*) ':stored',res-res_stored
!    endif


  end function f_bffbf_2


  recursive function f_bffbf_3(e,k,sp,p,fll,fl0,ng1,ng2,ng3,&
       &giarray,qiarray,pol_int) result(res)
    implicit none
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp(:,:),p(:,:)
    integer, intent(in) ::  ng1,ng2,ng3
    character, intent(in) :: fll(:)*3
    character, intent(in) :: fl0*3   ! flavor off-shell f-line
    integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int 
    ! -----------------------------------------------------------------------
    character :: fl1*3,fl2*3,fl3*3
!    integer             :: m1,m2,m3, ms1a, ms2a
    integer :: ngluon, ng4, ngL,m
    integer, parameter :: Ndumm=0
    double complex             :: res(size(sp,dim=1))
    double complex             :: res_stored(size(sp,dim=1))
    double complex             :: tmp(size(sp,dim=1))
    double complex             :: k1(size(k,dim=1))
    double complex             :: k2(size(k,dim=1))
!    double complex             :: k3(size(k,dim=1))
!    double complex             :: k4(size(k,dim=1))
    double complex             :: sp1(size(sp,dim=1))
    double complex             :: sp2(size(sp,dim=1))
!    double complex             :: sp3(size(sp,dim=1))
!    double complex             :: sp4(size(sp,dim=1))
    double complex             :: e1(size(e,dim=1))
    double complex             :: e2(size(e,dim=1))
!    double complex             :: e3(size(e,dim=1))
    double complex             :: kdumm(size(k,dim=1),Ndumm)
    double complex             :: edumm(size(e,dim=1),Ndumm)
    double complex  :: k1sq,k2sq!,k3sq!,k4sq
!    real(dp) :: mass,mass2
    logical                   :: done

    if (verbose) write(*,*) 'entering f_bffbf_3',ng1,ng2,ng3
    if (verbose .and. present(giarray)) write(*,*) 'entering f_bffbf_3:g',giarray
    if (verbose .and. present(giarray)) write(*,*) 'entering f_bffbf_3:q',qiarray
    done = .false. 
!    if (present(giarray)) then 
!       if (size(qiarray) /= 3) stop 'f_bffbf_3: wrong size qiarray'
!       if (size(e,dim=2) /= size(giarray)) stop 'bf_fbf: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray)
!       if (done) return 
!!       if (done) res_stored = res 
!   else
!      if (i_warn < max_warn) then 
!          write(*,*) 'f_bffbf_3: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif

!    mass = mt
!    mass2 = mt**2

    ngluon = size(e,dim=2)
    ng4 = ngluon - ng1 - ng2-ng3

    !if (verbose) write(*,*) 'in function f_bffbf_3', ng1, ngluon

    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)


    if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C:f_bffbf_3'

    if (ngluon == 0) then 

       res = czero
!       write(*,*) 'call 1',giarray,fl1,fl2,qiarray
       e2 = g_bff(edumm,kdumm,sp(:,1),p(:,1),fl1,&
            &sp(:,2),p(:,2),fl2,0,0,giarray,qiarray(1:2),pol_int)
       k2 = p(:,1)+p(:,2) 
       k2sq = sc(k2,k2)


       if (abs(k2sq) > propcut.and.fl0==fl3) then 
          tmp = -ci/k2sq*vgq(e2,sp(:,3))
       else 
          tmp = czero 
       endif

       res = res + tmp

    else                          !for #gluons > 0

       res = czero


       do m=0,ng4-1

          ngL = ng1+ ng2+ng3+m      

          sp1=f_bffbf_3(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0,&
               &ng1,ng2,ng3,giarray(1:ngL),qiarray,pol_int)
          if (1<=ngL) then 
             k1=sum(k(:,1:ngL),dim=2)
          else
             k1 = czero 
          end if
          k1 = k1 + p(:,1)+p(:,2)+p(:,3)
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


       do m=1,ng1

          e1 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)
          k1=sum(k(:,1:m),dim=2)
          k1sq=sc(k1,k1)

          sp2=f_bffbf_3(e(:,m+1:ngluon),k(:,m+1:ngluon),&
               &sp,p,fll,fl0,ng1-m,ng2,ng3,&
               &giarray(m+1:ngluon),qiarray,pol_int)
          if (m+1<=ngluon) then 
             k2=sum(k(:,m+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2 = k2 + p(:,1) + p(:,2) + p(:,3)
          k2sq = sc(k2,k2) !- mass2
          sp2 = spb2(sp2,k2) !+mass*sp2

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



       do m=0,ng3

          ngL = ng1+ ng2+m      

!       write(*,*) 'call 2',ng1,ng2,m,ng3
          e1 = g_bff(e(:,1:ngL),k(:,1:ngL),sp(:,1),p(:,1),fl1,&
               &sp(:,2),p(:,2),fl2,ng1,ng2,giarray(1:ngL),qiarray(1:2),pol_int)
          if (1<=ngL) then 
             k1=sum(k(:,1:ngL),dim=2)+p(:,1)+p(:,2)
          else
             k1 = p(:,1)+p(:,2)
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

          res = res + tmp

       enddo  !#4



    endif

!    if (done .and. any(res-res_stored /= czero)) then 
!       write(*,*) 'res_3:stored',res
!       write(*,*) 'res_3:stored',res_stored
!       write(*,*) 'res_3:stored',res-res_stored
!       write(*,*) 'stored',ng1,ng2,ng3,ngluon-ng1-ng2-ng3
!       if (present(giarray)) write(*,*) 'giarray',giarray
!       write(*,*) 'qiarray',qiarray
!    endif


    !! -- store current 
!    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray)


  end function f_bffbf_3


  recursive function bf_fbff(e,k,sp,p,fll,fl0,ng1,ng2,ng3,&
       &giarray,qiarray,pol_int) result(res)
    implicit none
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp(:,:),p(:,:)
    integer, intent(in) ::  ng1,ng2,ng3
    character, intent(in) :: fll(:)*3
    character, intent(in) :: fl0*3   ! flavor off-shell f-line
    integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int 
    ! -----------------------------------------------------------------------
    character :: fl1*3,fl2*3,fl3*3
!    integer             :: m1,m2,m3, ms1a, ms2a,m
    integer :: ngluon, ng4, ngL,m
    integer, parameter :: Ndumm=0
    double complex             :: res(size(sp,dim=1))
    double complex             :: tmp(size(sp,dim=1))
    double complex             :: k1(size(k,dim=1))
    double complex             :: k2(size(k,dim=1))
    double complex             :: sp1(size(sp,dim=1))
    double complex             :: sp2(size(sp,dim=1))
    double complex             :: e1(size(e,dim=1))
    double complex             :: e2(size(e,dim=1))
    double complex             :: kdumm(size(k,dim=1),Ndumm)
    double complex             :: edumm(size(e,dim=1),Ndumm)
    double complex             :: k1sq,k2sq
!    real(dp) :: mass,mass2
    logical                   :: done 

    if (verbose) write(*,*) 'entering bf_fbff'
    done = .false. 
!    if (present(giarray)) then 
!       if (size(qiarray) /= 3) stop 'bf_fbff: wrong size qiarray'
!       if (size(e,dim=2) /= size(giarray)) stop 'bf_fbf: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray)
!       if (done) return 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'bf_fbff: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif

!    mass = mt
!    mass2 = mt**2

    ngluon = size(e,dim=2)
    ng4 = ngluon - ng1 - ng2-ng3

    !if (verbose) write(*,*) 'in function bf_fbff', ng1, ngluon

    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)


    if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C:bf_fbff'

    if (ngluon == 0) then 

       res = czero

       !write(*,*) 'call 3'
       e2 = g_bff(edumm,kdumm,sp(:,2),p(:,2),fl2,&
            &sp(:,3),p(:,3),fl3,0,0,giarray,qiarray(2:3),pol_int)

       k1 = p(:,2)+p(:,3)
       k1sq=sc(k1,k1)


       if (abs(k1sq) > propcut.and.fl0==fl1) then 
          tmp = -ci/k1sq*vbqg(sp(:,1),e2)
       else 
          tmp = czero 
       endif

       res = res + tmp 
       e2 = g_fbf(edumm,kdumm,sp(:,1),p(:,1),fl1,&
            &sp(:,2),p(:,2),fl2,0,0,giarray,qiarray(1:2),pol_int)
       k2 = p(:,1)+p(:,2) 
       k2sq = sc(k2,k2)


       if (abs(k2sq) > propcut.and.fl0==fl3) then 
          tmp = -ci/k2sq*vgbq(e2,sp(:,3))
       else 
          tmp = czero 
       endif

       res = res + tmp

    else                          !for #gluons > 0


       res = czero

       do m=0,ng2
          !write(*,*) 'call bf from _3' 
          sp1=bf(e(:,1:ng1+m),k(:,1:ng1+m),sp(:,1),p(:,1),&
               &fl1,fl0,ng1,giarray(1:ng1+m),qiarray(1:1),pol_int)
          if (1<=ng1+m) then 
             k1=sum(k(:,1:ng1+m),dim=2)
          else
             k1 = czero 
          endif
          k1 = -k1 - p(:,1)
          k1sq = sc(k1,k1) !- mass**2

          if (ng1 > 0.or.m>0) then 
             sp1 = spi2(k1,sp1)!+mass*sp1
          endif

          e2 = g_bff(e(:,ng1+m+1:ngluon),k(:,ng1+m+1:ngluon),&
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
             tmp = -ci/k2sq*vbqg(sp1,e2)
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

       enddo  !#1

       do m=0,ng4-1

          ngL = ng1+ ng2+ng3+m      

          sp1=bf_fbff(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0,&
               &ng1,ng2,ng3,giarray(1:ngL),qiarray,pol_int)
          if (1<=ngL) then 
             k1=sum(k(:,1:ngL),dim=2)
          else
             k1 = czero 
          endif
          k1 = -k1 - p(:,1)-p(:,2)-p(:,3)
          k1sq = sc(k1,k1) !- mass2

          sp1 = spi2(k1,sp1) !+ mass*sp1

          e2 = vgluon(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
               &giarray(ngL+1:ngluon),pol_int)
          if (ngL+1<=ngluon) then 
             k2=sum(k(:,ngL+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2sq=sc(k2,k2)

          if (abs(k1sq) > propcut) then 
             tmp = ci/k1sq*vbqg(sp1,e2)
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

       do m=1,ng1

          e1 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)
          k1=sum(k(:,1:m),dim=2)
          k1sq=sc(k1,k1)

          sp2=bf_fbff(e(:,m+1:ngluon),k(:,m+1:ngluon),&
               &sp,p,fll,fl0,ng1-m,ng2,ng3,&
               &giarray(m+1:ngluon),qiarray,pol_int)
          if (m+1<=ngluon) then 
             k2=sum(k(:,m+1:ngluon),dim=2)
          else
             k2 = czero 
          endif
          k2 = -k2 - p(:,1) - p(:,2) - p(:,3)
          k2sq = sc(k2,k2) !- mass2
          sp2 = spi2(k2,sp2)! + mass*sp2

          if (abs(k2sq) > propcut) then 
             tmp = ci/k2sq*vgbq(e1,sp2)
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

       do m=0,ng3

          ngL = ng1+ ng2+m      

          e1 = g_fbf(e(:,1:ngL),k(:,1:ngL),sp(:,1),p(:,1),fl1,&
               &sp(:,2),p(:,2),fl2,ng1,ng2,giarray(1:ngL),qiarray(1:2),pol_int)
          if (1<=ngL) then 
             k1= sum(k(:,1:ngL),dim=2)+p(:,1)+p(:,2)
          else
             k1 = p(:,1)+p(:,2)
          endif
          k1sq=sc(k1,k1)


          sp2=bf(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
               &sp(:,3),p(:,3),fl3,fl0,ng3-m,&
               &giarray(ngL+1:ngluon),qiarray(3:3),pol_int)
          if (ngL+1<=ngluon) then 
             k2=sum(k(:,ngL+1:ngluon),dim=2)
          else
             k2 = czero
          endif
          k2 = -k2 - p(:,3)
          k2sq = sc(k2,k2) !- mass2

          if (ng4 > 0.or. m < ng3) then 
             sp2 = spi2(k2,sp2)!+mass*sp2
          endif


          if (abs(k1sq) > propcut) then 
             tmp = -ci/k1sq*vgbq(e1,sp2)
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

          res = res + tmp

       enddo  !#4




    endif

    !! -- store current 
!    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray)


  end function bf_fbff

!--------------------------------


recursive function bf_fbff_2(e,k,sp,p,fll,fl0,ng1,ng2,ng3,&
       &giarray,qiarray,pol_int) result(res)
    implicit none
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp(:,:),p(:,:)
    integer, intent(in) ::  ng1,ng2,ng3
    character, intent(in) :: fll(:)*3
    character, intent(in) :: fl0*3   ! flavor off-shell f-line
    integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int 
    ! -----------------------------------------------------------------------
    character :: fl1*3,fl2*3,fl3*3
!    integer             :: m1,m2,m3, ms1a, ms2a,m
    integer :: ngluon, ng4, ngL,m
    integer, parameter :: Ndumm=0
    double complex             :: res(size(sp,dim=1))
    double complex             :: tmp(size(sp,dim=1))
    double complex             :: k1(size(k,dim=1))
    double complex             :: k2(size(k,dim=1))
    double complex             :: sp1(size(sp,dim=1))
    double complex             :: sp2(size(sp,dim=1))
    double complex             :: e1(size(e,dim=1))
    double complex             :: e2(size(e,dim=1))
    double complex             :: kdumm(size(k,dim=1),Ndumm)
    double complex             :: edumm(size(e,dim=1),Ndumm)
    double complex  :: k1sq,k2sq
!    real(dp) :: mass,mass2
    logical                   :: done 

    if (verbose) write(*,*) 'entering bf_fbff'
    done = .false. 
!    if (present(giarray)) then 
!       if (size(qiarray) /= 3) stop 'bf_fbff: wrong size qiarray'
!       if (size(e,dim=2) /= size(giarray)) stop 'bf_fbf: ng= size(giarray)' 
!       ! XXX 
!       call memory_check(pol_int,res,done,giarray,qiarray)
!       if (done) return 
!    else
!       write(*,*) 'giarray missing'
!    endif

!    mass = mt
!    mass2 = mt**2

    ngluon = size(e,dim=2)
    ng4 = ngluon - ng1 - ng2-ng3

    !if (verbose) write(*,*) 'in function bf_fbff', ng1, ngluon

    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)


    if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C:bf_fbff'

    if (ngluon == 0) then 

       res = czero
 
       e2 = g_fbf(edumm,kdumm,sp(:,1),p(:,1),fl1,&
            &sp(:,2),p(:,2),fl2,0,0,giarray,qiarray(1:2),pol_int)
       k2 = p(:,1)+p(:,2) 
       k2sq = sc(k2,k2)


       if (abs(k2sq) > propcut.and.fl0==fl3) then 
          tmp = -ci/k2sq*vgbq(e2,sp(:,3))
       else 
          tmp = czero 
       endif

       res = res + tmp

    else                          !for #gluons > 0
       write(*,*) 'error in bf_fbff: current not written for ngluons > 0'
    endif
    
end function bf_fbff_2


  ! ---- this is the current for gluon splitting into 
  ! ----- two fermion pairs, of different flavors
  recursive function g_sbsfbf(e,k,sp,p,fll,&
       &ng1,ng2,ng3,ng4,&
       &giarray,qiarray,pol_int) result(res)
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp(:,:), p(:,:)
    integer, intent(in) ::  ng1,ng2, ng3,ng4
    character, intent(in) :: fll(:)*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int 
    ! -----------------------------------------------------------------------
    character :: fl1*3,fl2*3,fl3*3,fl4*3
    integer             :: m1,m2!,m3, ms1a, ms2a,Dv!,m
    integer :: ngluon,ng5,Dv,ms1a!,ngL
    integer, parameter :: Ndumm = 0
    double complex             :: res(size(e,dim=1))
    double complex             :: tmp(size(e,dim=1))
    double complex             :: k1(size(p,dim=1))
    double complex             :: k2(size(p,dim=1))
    double complex             :: k3(size(p,dim=1))
!    double complex             :: k4(size(p,dim=1))
    double complex             :: sp1(size(sp,dim=1))
    double complex             :: sp2(size(sp,dim=1))
!    double complex             :: sp3(size(sp,dim=1))
!    double complex             :: sp4(size(sp,dim=1))
    double complex             :: e1(size(e,dim=1))
    double complex             :: e2(size(e,dim=1))
    double complex             :: e3(size(e,dim=1))
    double complex             :: kdumm(size(k,dim=1),Ndumm)
    double complex             :: edumm(size(e,dim=1),Ndumm)
    double complex  :: k1sq,k2sq,k3sq!,k4sq
!    real(dp) :: mass,mass2
    logical                   :: done 

    if (verbose) write(*,*) 'entering g_sbsfbf'
    done = .false. 
!    if (present(giarray)) then 
!       if (size(qiarray) /= 4) stop 'g_sbsfbf: wrong size qiarray'
!       if (size(e,dim=2) /= size(giarray)) stop 'g_sbsfbf: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray)
!       if (done) return 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'g_sbsfbf: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif



!    mass = mt
!    mass2= mass**2

    Dv = size(e,dim=1)

    ngluon = size(e,dim=2)

    ng5 = ngluon - ng1 - ng2-ng3 - ng4

    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)
    fl4 = fll(4)



    if (ng5 < 0) write(*,*) 'ERROR IN CURRENT G'


    if ((fl1.eq.fl2).and.(fl3.eq.fl4)) then 

       if (ngluon == 0) then 

          res = czero

          e1 = g_fbf(edumm,kdumm,sp(:,1),p(:,1),fl1,&
               &sp(:,2),p(:,2),fl2,0,0,giarray,qiarray(1:2),pol_int)
          k1 = p(:,1) + p(:,2) 
          k1sq = sc(k1,k1)

          e2 = g_fbf(edumm,kdumm,sp(:,3),p(:,3),fl3,&
               &sp(:,4),p(:,4),fl4,0,0,giarray,qiarray(3:4),pol_int)
          k2 = p(:,3) + p(:,4)
          k2sq = sc(k2,k2)

          tmp = (-ci/k1sq)*(-ci/k2sq)*vggg(e1,k1,e2,k2)

          res = res + tmp  !# 1


          sp2 = f_bffbf(edumm,kdumm,sp(:,2:4),p(:,2:4),&
               &fll(2:4),fl1,0,0,0,giarray,qiarray(2:4),pol_int)
          k2 = p(:,2) + p(:,3) + p(:,4)
          k2sq = sc(k2,k2)

          sp2 = spb2(sp2,k2) !+ mass*sp2
          sp2 = ci/k2sq*sp2

          sp1 = sp(:,1)

          tmp = -cone*vbqq(Dv,sp2,sp1)

          res = res + tmp       ! # 2

          sp2 = bf_fbff(edumm,kdumm,sp(:,1:3),p(:,1:3),&
               &fll(1:3),fl4,0,0,0,giarray,qiarray(1:3),pol_int)

          k2 = -p(:,1) - p(:,2) - p(:,3)
          k2sq = sc(k2,k2)

          sp2 = spi2(k2,sp2) !+ mass*sp2
          sp2 = ci/k2sq*sp2

          sp1 = sp(:,4)

          tmp = -cone*vbqq(Dv,sp1,sp2)

          res = res + tmp       ! # 3


       else  ! for ngluon > 0

          res = czero

          do m1=1,ng1

             e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)
             k1=sum(k(:,1:m1),dim=2)
             k1sq = sc(k1,k1)

             ms1a = ng1-m1
             e2=g_sbsfbf(e(:,m1+1:ngluon),k(:,m1+1:ngluon),&
                  &sp,p,fll,ms1a,ng2,ng3,ng4,giarray(m1+1:ngluon),qiarray,pol_int)
             if (m1+1<=ngluon) then 
                k2 = sum(k(:,m1+1:ngluon),dim=2)
             else
                k2 = czero 
             endif
             k2 = k2 + p(:,1) + p(:,2) + p(:,3)+p(:,4)
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

          do m1=0,ng3

             e1=g_fbf(e(:,1:ng1+ng2+m1),k(:,1:ng1+ng2+m1),&
                  &sp(:,1),p(:,1),fl1,sp(:,2),p(:,2),fl2,ng1,ng2,&
                  &giarray(1:ng1+ng2+m1),qiarray(1:2),pol_int)
             if (1<=ng1+ng2+m1) then 
                k1=sum(k(:,1:ng1+ng2+m1),dim=2)
             else
                k1 = czero 
             endif
             k1 = k1 + p(:,1)+p(:,2)
             k1sq = sc(k1,k1)

             e2=g_fbf(e(:,ng1+ng2+m1+1:ngluon),&
                  &k(:,ng1+ng2+m1+1:ngluon),    &
                  &sp(:,3),p(:,3),fl3,sp(:,4),p(:,4),fl4,ng3-m1,ng4,&
                  &giarray(ng1+ng2+m1+1:ngluon),qiarray(3:4),pol_int)
             if (ng1+ng2+m1+1<=ngluon) then 
                k2 = sum(k(:,ng1+ng2+m1+1:ngluon),dim=2)
             else
                k2 = czero 
             endif
             k2 = k2 + p(:,3) + p(:,4) 
             k2sq = sc(k2,k2)


             tmp = (-ci/k2sq)*(-ci/k1sq)*vggg(e1,k1,e2,k2)

             res = res + tmp



          enddo  !#2


          do m1=1,ng5

             e1 = g_sbsfbf(e(:,1:ng1+ng2+ng3+ng4+m1-1),&
                  &k(:,1:ng1+ng2+ng3+ng4+m1-1),sp,p,fll,ng1,ng2,ng3,ng4,&
                  &giarray(1:ng1+ng2+ng3+ng4+m1-1),qiarray,pol_int)
             if (1<=ng1+ng2+ng3+ng4+m1-1) then 
                k1 = sum(k(:,1:ng1+ng2+ng3+ng4+m1-1),dim=2)
             else
                k1 = czero 
             endif
             k1 = k1 + p(:,1)+p(:,2)+p(:,3) + p(:,4)
             k1sq = sc(k1,k1)

             e2 = vgluon(e(:,ng1+ng2+ng3+ng4+m1:ngluon),&
                  &k(:,ng1+ng2+ng3+ng4+m1:ngluon),&
                  &giarray(ng1+ng2+ng3+ng4+m1:ngluon),pol_int)
             if (ng1+ng2+ng3+ng4+m1<=ngluon) then 
                k2 =sum(k(:,ng1+ng2+ng3+ng4+m1:ngluon),dim=2) 
             else
                k2 = czero 
             endif
             k2sq=sc(k2,k2)

             tmp = -ci/k1sq*vggg(e1,k1,e2,k2)

             if (ng1+ng2+ng3+ng4+m1.ne.ngluon) then 
                tmp = -ci/k2sq*tmp
             endif

             res = res + tmp



          enddo   !   #3



          do m1=1,ng1-1

             e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)
             k1=sum(k(:,1:m1),dim=2)
             k1sq = sc(k1,k1)


             do m2= m1+1,ng1

                e2=vgluon(e(:,m1+1:m2),k(:,m1+1:m2),giarray(m1+1:m2),pol_int)
                if (m1+1<=m2) then 
                   k2=sum(k(:,m1+1:m2),dim=2)
                else
                   k2 = czero 
                endif
                k2sq = sc(k2,k2)



                e3 = g_sbsfbf(e(:,m2+1:ngluon),&
                     &k(:,m2+1:ngluon),sp,p,fll,ng1-m2,ng2,ng3,ng4,&
                     &giarray(m2+1:ngluon),qiarray,pol_int)
                if (m2+1<=ngluon) then 
                   k3 = sum(k(:,m2+1:ngluon),dim=2)
                else
                   k3 = czero 
                endif
                k3 = k3 + p(:,1)+p(:,2)+p(:,3) + p(:,4)
                k3sq = sc(k3,k3)


                tmp = -ci/k3sq*vgggg(e1,e2,e3)

                if (m1 > 1) then 
                   tmp = -ci/k1sq*tmp
                endif

                if (m2 > m1+1) then 
                   tmp = -ci/k2sq*tmp
                endif

                res = res + tmp



             enddo

          enddo   !   #4




          do m1=1,ng1

             e1=vgluon(e(:,1:m1),k(:,1:m1))
             k1=sum(k(:,1:m1),dim=2)
             k1sq = sc(k1,k1)


             do m2= 0,ng3

                e2=g_fbf(e(:,m1+1:ng1+ng2+m2),k(:,m1+1:ng1+ng2+m2),&
                     &sp(:,1),p(:,1),fl1,sp(:,2),p(:,2),fl2,ng1-m1,ng2,&
                     &giarray(m1+1:ng1+ng2+m2),qiarray(1:2),pol_int)
                if (m1+1<=ng1+ng2+m2) then 
                   k2=sum(k(:,m1+1:ng1+ng2+m2),dim=2)
                else
                   k2 = czero 
                endif
                k2 = k2 + p(:,1) + p(:,2)
                k2sq = sc(k2,k2)


                e3 = g_fbf(e(:,ng1+ng2+m2+1:ngluon),&
                     &k(:,ng1+ng2+m2+1:ngluon),sp(:,3),p(:,3),fl3,&
                     &sp(:,4),p(:,4),fl4,ng3 - m2,ng4,&
                     &giarray(ng1+ng2+m2+1:ngluon),qiarray(3:4),pol_int)
                if (ng1+ng2+m2+1<=ngluon) then 
                   k3 = sum(k(:,ng1+ng2+m2+1:ngluon),dim=2)
                else
                   k3 = czero 
                endif
                k3 = k3 + p(:,3) + p(:,4)
                k3sq = sc(k3,k3)


                tmp = -ci/k3sq*(-ci)/k2sq*vgggg(e1,e2,e3)

                if (m1 > 1) then 
                   tmp = -ci/k1sq*tmp
                endif

                res = res + tmp



             enddo

          enddo   !   #5


          do m1=1,ng1

             e1=vgluon(e(:,1:m1),k(:,1:m1))
             k1=sum(k(:,1:m1),dim=2)
             k1sq = sc(k1,k1)


             do m2= 1,ng5-1

                e2=g_sbsfbf(e(:,m1+1:ng1+ng2+ng3+ng4+m2),&
                     &k(:,m1+1:ng1+ng2+ng3+ng4+m2),&
                     &sp,p,fll,ng1-m1,ng2,ng3,ng4,&
                     &giarray(m1+1:ng1+ng2+ng3+ng4+m2),qiarray,pol_int)
                if (m1+1<=ng1+ng2+ng3+ng4+m2) then 
                   k2=sum(k(:,m1+1:ng1+ng2+ng3+ng4+m2),dim=2)
                else
                   k2 = czero 
                endif
                k2 = k2 + p(:,1)+p(:,2)+p(:,3)+p(:,4)
                k2sq = sc(k2,k2)



                e3 = vgluon(e(:,ng1+ng2+ng3+ng4+m2+1:ngluon),&
                     &k(:,ng1+ng2+ng3+ng4+m2+1:ngluon),&
                     &giarray(ng1+ng2+ng3+ng4+m2+1:ngluon),pol_int)
                if (ng1+ng2+ng3+ng4+m2+1<=ngluon) then 
                   k3 = sum(k(:,ng1+ng2+ng3+ng4+m2+1:ngluon),dim=2)
                else
                   k3 = czero 
                endif
                k3sq = sc(k3,k3)


                tmp = -ci/k2sq*vgggg(e1,e2,e3)

                if (m1 > 1) then 
                   tmp = -ci/k1sq*tmp
                endif

                if (m2 < ng5-1) then 
                   tmp = -ci/k3sq*tmp 
                endif

                res = res + tmp

             enddo

          enddo   !   #6


          do m1=0,ng3-1

             e1=g_fbf(e(:,1:ng1+ng2+m1),k(:,1:ng1+ng2+m1),&
                  &sp(:,1),p(:,1),fl1,sp(:,2),p(:,2),fl2,ng1,ng2,&
                  &giarray(1:ng1+ng2+m1),qiarray(1:2),pol_int)
             if (1<=ng1+ng2+m1) then 
                k1=sum(k(:,1:ng1+ng2+m1),dim=2)
             else
                k1 = czero
             endif
             k1 = k1 + p(:,1)+p(:,2)
             k1sq = sc(k1,k1)

             do m2= m1+1,ng3

                e2=vgluon(e(:,ng1+ng2+m1+1:ng1+ng2+m2),&
                     &k(:,ng1+ng2+m1+1:ng1+ng2+m2),&
                     &giarray(ng1+ng2+m1+1:ng1+ng2+m2),pol_int)
                if (ng1+ng2+m1+1<=ng1+ng2+m2) then 
                   k2=sum(k(:,ng1+ng2+m1+1:ng1+ng2+m2),dim=2)
                else
                   k2 = czero 
                endif
                k2sq = sc(k2,k2)



                e3 = g_fbf(e(:,ng1+ng2+m2+1:ngluon),&
                     &k(:,ng1+ng2+m2+1:ngluon),sp(:,3),p(:,3),fl3,&
                     &sp(:,4),p(:,4),fl4,ng3-m2,ng4,&
                     &giarray(ng1+ng2+m2+1:ngluon),qiarray(3:4),pol_int)
                if (ng1+ng2+m2+1<=ngluon) then 
                   k3 = sum(k(:,ng1+ng2+m2+1:ngluon),dim=2)
                else
                   k3 = czero 
                endif
                k3 = k3 + p(:,3)+p(:,4)
                k3sq = sc(k3,k3)

                tmp = -ci/k3sq*(-ci)/k1sq*vgggg(e1,e2,e3)


                if (m2 > m1+1) then 
                   tmp = -ci/k2sq*tmp 
                endif

                res = res + tmp

             enddo

          enddo   !   #7


          do m1=0,ng3

             e1=g_fbf(e(:,1:ng1+ng2+m1),k(:,1:ng1+ng2+m1),&
                  &sp(:,1),p(:,1),fl1,sp(:,2),p(:,2),fl2,ng1,ng2,&
                  &giarray(1:ng1+ng2+m1),qiarray(1:2),pol_int)
             if (1<=ng1+ng2+m1) then 
                k1=sum(k(:,1:ng1+ng2+m1),dim=2)
             else
                k1 = czero 
             endif
             k1 = k1 + p(:,1)+p(:,2)
             k1sq = sc(k1,k1)

             do m2= 0,ng5-1

                e2=g_fbf(e(:,ng1+ng2+m1+1:ng1+ng2+ng3+ng4+m2),&
                     &k(:,ng1+ng2+m1+1:ng1+ng2+ng3+ng4+m2),sp(:,3),p(:,3),fl3,&
                     &sp(:,4),p(:,4),fl4,ng3-m1,ng4,&
                     &giarray(ng1+ng2+m1+1:ng1+ng2+ng3+ng4+m2),&
                     &qiarray(3:4),pol_int)
                if (ng1+ng2+m1+1<=ng1+ng2+ng3+ng4+m2) then 
                   k2=sum(k(:,ng1+ng2+m1+1:ng1+ng2+ng3+ng4+m2),dim=2)
                else
                   k2 = czero 
                endif
                k2 = k2 + p(:,3)+p(:,4)
                k2sq = sc(k2,k2)


                e3 = vgluon(e(:,ng1+ng2+ng3+ng4+m2+1:ngluon),&
                     &k(:,ng1+ng2+ng3+ng4+m2+1:ngluon),&
                     &giarray(ng1+ng2+ng3+ng4+m2+1:ngluon),pol_int)
                if (ng1+ng2+ng3+ng4+m2+1<=ngluon) then 
                   k3 = sum(k(:,ng1+ng2+ng3+ng4+m2+1:ngluon),dim=2)
                else
                   k3 = czero 
                endif
                k3sq = sc(k3,k3)

                tmp = -ci/k1sq*(-ci)/k2sq*vgggg(e1,e2,e3)


                if (m2 < ng5-1) then 
                   tmp = -ci/k3sq*tmp 
                endif

                res = res + tmp

             enddo

          enddo   !   #8



          do m1=0,ng5-2

             e1=g_sbsfbf(e(:,1:ng1+ng2+ng3+ng4+m1),&
                  &k(:,1:ng1+ng2+ng3+ng4+m1),&
                  &sp,p,fll,ng1,ng2,ng3,ng4,&
                  &giarray(1:ng1+ng2+ng3+ng4+m1),qiarray,pol_int)
             if (1<=ng1+ng2+ng3+ng4+m1) then 
                k1=sum(k(:,1:ng1+ng2+ng3+ng4+m1),dim=2)
             else
                k1 = czero 
             endif
             k1 = k1 + p(:,1)+p(:,2)+p(:,3)+p(:,4)
             k1sq = sc(k1,k1)

             do m2= m1+1,ng5-1

                e2=vgluon(e(:,ng1+ng2+ng3+ng4+m1+1:ng1+ng2+ng3+ng4+m2),&
                     &k(:,ng1+ng2+m1+1:ng1+ng2+ng3+ng4+m2),&
                     &giarray(ng1+ng2+ng3+ng4+m1+1:ng1+ng2+ng3+ng4+m2),pol_int)
                if (ng1+ng2+m1+1<=ng1+ng2+ng3+ng4+m2) then 
                   k2=sum(k(:,ng1+ng2+m1+1:ng1+ng2+ng3+ng4+m2),dim=2)
                else
                   k2 = czero 
                endif
                k2sq = sc(k2,k2)

                e3 = vgluon(e(:,ng1+ng2+ng3+ng4+m2+1:ngluon),&
                     &k(:,ng1+ng2+ng3+ng4+m2+1:ngluon),&
                     &giarray(ng1+ng2+ng3+ng4+m2+1:ngluon),pol_int)
                if (ng1+ng2+ng3+ng4+m2+1<=ngluon) then 
                   k3 = sum(k(:,ng1+ng2+ng3+ng4+m2+1:ngluon),dim=2)
                else
                   k3 = czero 
                endif
                k3sq = sc(k3,k3)

                tmp = -ci/k1sq*vgggg(e1,e2,e3)

                if (m2 > m1+1) then 
                   tmp = -ci/k2sq*tmp
                endif

                if (m2 < ng5-1) then 
                   tmp = -ci/k3sq*tmp 
                endif

                res = res + tmp

             enddo

          enddo   !   #9


          !-------here 
          do m1=0,ng4

             sp2 = bf_fbff(e(:,1:ng1+ng2+ng3+m1),&
                  &k(:,1:ng1+ng2+ng3+m1),&
                  &sp(:,1:3),p(:,1:3),fll(1:3),fl4,ng1,ng2,ng3,&
                  &giarray(1:ng1+ng2+ng3+m1),qiarray(1:3),pol_int)
             if (1<=ng1+ng2+ng3+m1) then 
                k2=sum(k(:,1:ng1+ng2+ng3+m1),dim=2)
             else
                k2 = czero 
             endif
             k2 = -k2 -p(:,1) - p(:,2) - p(:,3)
             k2sq = sc(k2,k2)

             sp2 = spi2(k2,sp2) !+ mass*sp2
             sp2 = ci/k2sq*sp2

             sp1 = f(e(:,ng1+ng2+ng3+m1+1:ngluon),&
                  &k(:,ng1+ng2+ng3+m1+1:ngluon),sp(:,4),p(:,4),&
                  &fl4,fl3,ng4-m1,&
                  &giarray(ng1+ng2+ng3+m1+1:ngluon),qiarray(4:4),pol_int)

             if (ng1+ng2+ng3+m1+1<=ngluon) then 
                k1=sum(k(:,ng1+ng2+ng3+m1+1:ngluon),dim=2)
             else
                k1 = czero 
             endif
             k1 = k1+p(:,4)
             k1sq = sc(k1,k1)

             if (ng5>0.or.m1< ng4) then 
                sp1 = spb2(sp1,k1) !+ mass*sp1
                sp1 = ci/k1sq*sp1
             endif

             tmp = -cone*vbqq(Dv,sp1,sp2)



             res = res + tmp       ! #10

          enddo



          do m1=0,ng2

             sp2 = f_bffbf(e(:,ng1+m1+1:ngluon),k(:,ng1+m1+1:ngluon),&
                  &sp(:,2:4),p(:,2:4),&
                  &fll(2:4),fl1,ng2-m1,ng3,ng4,&
                  &giarray(ng1+m1+1:ngluon),qiarray(2:4),pol_int)
             if (ng1+m1+1<=ngluon) then 
                k2=sum(k(:,ng1+m1+1:ngluon),dim=2)
             else
                k2 = czero 
             endif
             k2 = k2 + p(:,2) + p(:,3) + p(:,4)
             k2sq = sc(k2,k2)

             sp2 = spb2(sp2,k2) !+ mass*sp2
             sp2 = ci/k2sq*sp2

             sp1 = bf(e(:,1:ng1+m1),k(:,1:ng1+m1),sp(:,1),p(:,1),&
                  &fl1,fl2,ng1,&
                  &giarray(1:ng1+m1),qiarray(1:1),pol_int)
             if (1<=ng1+m1) then 
                k1 = sum(k(:,1:ng1+m1),dim=2)
             else
                k1 = czero 
             endif
             k1 = -k1 - p(:,1)
             k1sq = sc(k1,k1)

             if (ng1>0.or.m1>0) then
                sp1 = spi2(k1,sp1) !+ mass*sp1
                sp1 = ci/k1sq*sp1
             endif

             tmp = -cone*vbqq(Dv,sp2,sp1)



             res = res + tmp       ! #11


          enddo

       endif   !end condition for ngluon 

    endif !(for a flavor combination ) 



    if ((fl1.eq.fl4).and.(fl2.eq.fl3)) then 



       if (ngluon == 0) then 

          res = czero


          sp2 = f_bffbf(edumm,kdumm,sp(:,2:4),p(:,2:4),&
               &fll(2:4),fl1,0,0,0,&
               &giarray,qiarray(2:4),pol_int)
          k2 = p(:,2) + p(:,3) + p(:,4)
          k2sq = sc(k2,k2)

          sp2 = spb2(sp2,k2) !+ mass*sp2
          sp2 = ci/k2sq*sp2

          sp1 = sp(:,1)

          tmp = -cone*vbqq(Dv,sp2,sp1)

          res = res + tmp       ! # 2

          sp2 = bf_fbff(edumm,kdumm,sp(:,1:3),p(:,1:3),&
               &fll(1:3),fl4,0,0,0,&
               &giarray,qiarray(1:3),pol_int)

          k2 = -p(:,1) - p(:,2) - p(:,3)
          k2sq = sc(k2,k2)

          sp2 = spi2(k2,sp2) !+ mass*sp2
          sp2 = ci/k2sq*sp2

          sp1 = sp(:,4)

          tmp = -cone*vbqq(Dv,sp1,sp2)

          res = res + tmp       ! # 3


       else  ! for ngluon > 0


          res = czero      

          do m1=1,ng1

             e1=vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)
             k1=sum(k(:,1:m1),dim=2)
             k1sq = sc(k1,k1)

             ms1a = ng1-m1
             e2=g_sbsfbf(e(:,m1+1:ngluon),k(:,m1+1:ngluon),&
                  &sp,p,fll,ms1a,ng2,ng3,ng4,&
                  &giarray(m1+1:ngluon),qiarray,pol_int)
             if (m1+1<=ngluon) then 
                k2 = sum(k(:,m1+1:ngluon),dim=2)
             else
                k2 = czero 
             endif
             k2 = k2 + p(:,1) + p(:,2) + p(:,3)+p(:,4)
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




          do m1=1,ng5

             e1 = g_sbsfbf(e(:,1:ng1+ng2+ng3+ng4+m1-1),&
                  &k(:,1:ng1+ng2+ng3+ng4+m1-1),sp,p,fll,ng1,ng2,ng3,ng4,&
                  &giarray(1:ng1+ng2+ng3+ng4+m1-1),qiarray,pol_int)
             if (1<=ng1+ng2+ng3+ng4+m1-1) then 
                k1 = sum(k(:,1:ng1+ng2+ng3+ng4+m1-1),dim=2)
             else
                k1 = czero 
             endif
             k1 = k1 + p(:,1)+p(:,2)+p(:,3) + p(:,4)
             k1sq = sc(k1,k1)

             e2 = vgluon(e(:,ng1+ng2+ng3+ng4+m1:ngluon),&
                  &k(:,ng1+ng2+ng3+ng4+m1:ngluon),&
                  &giarray(ng1+ng2+ng3+ng4+m1:ngluon),pol_int)
             if (ng1+ng2+ng3+ng4+m1<=ngluon) then 
                k2 =sum(k(:,ng1+ng2+ng3+ng4+m1:ngluon),dim=2) 
             else
                k2 = czero 
             endif
             k2sq=sc(k2,k2)

             tmp = -ci/k1sq*vggg(e1,k1,e2,k2)

             if (ng1+ng2+ng3+ng4+m1.ne.ngluon) then 
                tmp = -ci/k2sq*tmp
             endif

             res = res + tmp

          enddo   !   #3




          do m1=1,ng1-1

             e1=vgluon(e(:,1:m1),k(:,1:m1))
             k1=sum(k(:,1:m1),dim=2)
             k1sq = sc(k1,k1)


             do m2= m1+1,ng1

                e2=vgluon(e(:,m1+1:m2),k(:,m1+1:m2),giarray(m1+1:m2),pol_int)
                if (m1+1<=m2) then 
                   k2=sum(k(:,m1+1:m2),dim=2)
                else
                   k2 = czero 
                endif
                k2sq = sc(k2,k2)



                e3 = g_sbsfbf(e(:,m2+1:ngluon),&
                     &k(:,m2+1:ngluon),sp,p,fll,ng1-m2,ng2,ng3,ng4,&
                     &giarray(m2+1:ngluon),qiarray,pol_int)
                if (m2+1<=ngluon) then 
                   k3 = sum(k(:,m2+1:ngluon),dim=2)
                else
                   k3 = czero 
                endif
                k3 = k3 + p(:,1)+p(:,2)+p(:,3) + p(:,4)
                k3sq = sc(k3,k3)


                tmp = -ci/k3sq*vgggg(e1,e2,e3)

                if (m1 > 1) then 
                   tmp = -ci/k1sq*tmp
                endif

                if (m2 > m1+1) then 
                   tmp = -ci/k2sq*tmp
                endif

                res = res + tmp



             enddo

          enddo   !   #4





          do m1=1,ng1

             e1=vgluon(e(:,1:m1),k(:,1:m1))
             k1=sum(k(:,1:m1),dim=2)
             k1sq = sc(k1,k1)


             do m2= 1,ng5-1

                e2=g_sbsfbf(e(:,m1+1:ng1+ng2+ng3+ng4+m2),&
                     &k(:,m1+1:ng1+ng2+ng3+ng4+m2),&
                     &sp,p,fll,ng1-m1,ng2,ng3,ng4,&
                     &giarray(m1+1:ng1+ng2+ng3+ng4+m2),qiarray,pol_int)
                if (m1+1<=ng1+ng2+ng3+ng4+m2) then 
                   k2=sum(k(:,m1+1:ng1+ng2+ng3+ng4+m2),dim=2)
                else
                   k2 = czero 
                endif
                k2 = k2 + p(:,1)+p(:,2)+p(:,3)+p(:,4)
                k2sq = sc(k2,k2)



                e3 = vgluon(e(:,ng1+ng2+ng3+ng4+m2+1:ngluon),&
                     &k(:,ng1+ng2+ng3+ng4+m2+1:ngluon),&
                     &giarray(ng1+ng2+ng3+ng4+m2+1:ngluon),pol_int)
                if (ng1+ng2+ng3+ng4+m2+1<=ngluon) then 
                   k3 = sum(k(:,ng1+ng2+ng3+ng4+m2+1:ngluon),dim=2)
                else
                   k3 = czero 
                endif
                k3sq = sc(k3,k3)


                tmp = -ci/k2sq*vgggg(e1,e2,e3)

                if (m1 > 1) then 
                   tmp = -ci/k1sq*tmp
                endif

                if (m2 < ng5-1) then 
                   tmp = -ci/k3sq*tmp 
                endif

                res = res + tmp

             enddo

          enddo   !   #6



          do m1=0,ng5-2

             e1=g_sbsfbf(e(:,1:ng1+ng2+ng3+ng4+m1),&
                  &k(:,1:ng1+ng2+ng3+ng4+m1),&
                  &sp,p,fll,ng1,ng2,ng3,ng4,&
                  &giarray(1:ng1+ng2+ng3+ng4+m1),qiarray,pol_int)
             if (1<=ng1+ng2+ng3+ng4+m1) then 
                k1=sum(k(:,1:ng1+ng2+ng3+ng4+m1),dim=2)
             else
                k1 = czero 
             endif
             k1 = k1 + p(:,1)+p(:,2)+p(:,3)+p(:,4)
             k1sq = sc(k1,k1)

             do m2= m1+1,ng5-1

                e2=vgluon(e(:,ng1+ng2+ng3+ng4+m1+1:ng1+ng2+ng3+ng4+m2),&
                     &k(:,ng1+ng2+m1+1:ng1+ng2+ng3+ng4+m2),&
                     &giarray(ng1+ng2+ng3+ng4+m1+1:ng1+ng2+ng3+ng4+m2),pol_int)
                if (ng1+ng2+m1+1<=ng1+ng2+ng3+ng4+m2) then 
                   k2=sum(k(:,ng1+ng2+m1+1:ng1+ng2+ng3+ng4+m2),dim=2)
                else
                   k2 = czero
                endif
                k2sq = sc(k2,k2)

                e3 = vgluon(e(:,ng1+ng2+ng3+ng4+m2+1:ngluon),&
                     &k(:,ng1+ng2+ng3+ng4+m2+1:ngluon),&
                     &giarray(ng1+ng2+ng3+ng4+m2+1:ngluon),pol_int)
                if (ng1+ng2+ng3+ng4+m2+1<=ngluon) then 
                   k3 = sum(k(:,ng1+ng2+ng3+ng4+m2+1:ngluon),dim=2)
                else
                   k3 = czero 
                endif
                k3sq = sc(k3,k3)

                tmp = -ci/k1sq*vgggg(e1,e2,e3)

                if (m2 > m1+1) then 
                   tmp = -ci/k2sq*tmp
                endif

                if (m2 < ng5-1) then 
                   tmp = -ci/k3sq*tmp 
                endif

                res = res + tmp

             enddo

          enddo   !   #9



          !-------here 

          do m1=0,ng4

             sp2 = bf_fbff(e(:,1:ng1+ng2+ng3+m1),&
                  &k(:,1:ng1+ng2+ng3+m1),&
                  &sp(:,1:3),p(:,1:3),fll(1:3),fl4,ng1,ng2,ng3,&
                  &giarray(1:ng1+ng2+ng3+m1),qiarray(1:3),pol_int)

             if (1<=ng1+ng2+ng3+m1) then 
                k2=sum(k(:,1:ng1+ng2+ng3+m1),dim=2)
             else
                k2 = czero 
             endif
             k2 = -k2 -p(:,1) - p(:,2) - p(:,3)
             k2sq = sc(k2,k2)

             sp2 = spi2(k2,sp2) !+ mass*sp2
             sp2 = ci/k2sq*sp2

             sp1 = f(e(:,ng1+ng2+ng3+m1+1:ngluon),&
                  &k(:,ng1+ng2+ng3+m1+1:ngluon),sp(:,4),p(:,4),&
                  &fl4,fl4,ng4-m1,&
                  &giarray(ng1+ng2+ng3+m1+1:ngluon),qiarray(4:4),pol_int)

             if (ng1+ng2+ng3+m1+1<=ngluon) then 
                k1=sum(k(:,ng1+ng2+ng3+m1+1:ngluon),dim=2)
             else
                k1 = czero 
             endif
             k1 = k1+p(:,4)
             k1sq = sc(k1,k1)

             if (ng5>0.or.m1< ng4) then 
                sp1 = spb2(sp1,k1) !+ mass*sp1
                sp1 = ci/k1sq*sp1
             endif

             tmp = -cone*vbqq(Dv,sp1,sp2)



             res = res + tmp       ! #10



          enddo


          do m1=0,ng2

             sp2 = f_bffbf(e(:,ng1+m1+1:ngluon),k(:,ng1+m1+1:ngluon),&
                  &sp(:,2:4),p(:,2:4),&
                  &fll(2:4),fl1,ng2-m1,ng3,ng4,&
                  &giarray(ng1+m1+1:ngluon),qiarray(2:4),pol_int)
             if (ng1+m1+1<=ngluon) then 
                k2=sum(k(:,ng1+m1+1:ngluon),dim=2)
             else
                k2 = czero 
             endif
             k2 = k2 + p(:,2) + p(:,3) + p(:,4)
             k2sq = sc(k2,k2)

             sp2 = spb2(sp2,k2) !+ mass*sp2
             sp2 = ci/k2sq*sp2

             sp1 = bf(e(:,1:ng1+m1),k(:,1:ng1+m1),sp(:,1),p(:,1),&
                  &fl1,fl1,ng1,&
                 & giarray(1:ng1+m1),qiarray(1:1),pol_int)
             if (1<=ng1+m1) then 
                k1 = sum(k(:,1:ng1+m1),dim=2)
             else
                k1 = czero 
             endif
             k1 = -k1 - p(:,1)
             k1sq = sc(k1,k1)

             if (ng1>0.or.m1>0) then
                sp1 = spi2(k1,sp1) !+ mass*sp1
                sp1 = ci/k1sq*sp1
             endif

             tmp = -cone*vbqq(Dv,sp2,sp1)


             res = res + tmp       ! #11

          enddo


       endif   !end condition for ngluon 

    endif ! for a particular flavor combination 

    ! -- store current 
!    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray)

  end function g_sbsfbf

  ! ---- this is the current for gluon splitting into 
  ! ----- two fermion pairs, of different flavors
  recursive function g_bssbff(e,k,sp,p,fll,&
       &ng1,ng2,ng3,ng4,&
       &giarray,qiarray,pol_int) result(res)
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp(:,:), p(:,:)
    integer, intent(in) ::  ng1,ng2, ng3,ng4
    character, intent(in) :: fll(:)*3
    integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int 
    ! -----------------------------------------------------------------------
    character :: fl1*3,fl2*3,fl3*3,fl4*3
    integer             :: m1,m2
    integer :: ngluon,ng5,Dv,ms1a
    integer, parameter :: Ndumm = 0
    double complex             :: res(size(e,dim=1))
    double complex             :: tmp(size(e,dim=1))
    double complex             :: k1(size(p,dim=1))
    double complex             :: k2(size(p,dim=1))
    double complex             :: k3(size(p,dim=1))
    double complex             :: sp1(size(sp,dim=1))
    double complex             :: sp2(size(sp,dim=1))
    double complex             :: e1(size(e,dim=1))
    double complex             :: e2(size(e,dim=1))
    double complex             :: e3(size(e,dim=1))
    double complex             :: kdumm(size(k,dim=1),Ndumm)
    double complex             :: edumm(size(e,dim=1),Ndumm)
    double complex             :: k1sq,k2sq,k3sq
!    real(dp)                :: mass,mass2
    logical                   :: done

    if (verbose) write(*,*) 'entering g_bssfbff'
    done = .false. 
!    if (present(giarray)) then 
!       if (size(qiarray) /= 4) stop 'g_bssbff: wrong size qiarray'
!       if (size(e,dim=2) /= size(giarray)) stop 'g_sbsfbf: ng= size(giarray)' 
!       call memory_check(pol_int,res,done,giarray,qiarray)
!       if (done) return 
!    else
!       if (i_warn < max_warn) then 
!          write(*,*) 'g_bssfbff: giarray missing', i_warn 
!          i_warn = i_warn+1
!       endif
!    endif

!    mass = mt
!    mass2= mass**2

    Dv = size(e,dim=1)

    ngluon = size(e,dim=2)

    if (ngluon > 0) stop 'g_bssbff: please supply current for sucessfull &
         &completion of this calculation...'

    ng5 = ngluon - ng1 - ng2-ng3 - ng4

    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)
    fl4 = fll(4)


    if ((ng1 < 0) .or. (ng2 < 0) .or. (ng3 < 0) .or. (ng4 < 0) .or. (ng5 < 0)) &
         &write(*,*) 'ERROR IN CURRENT g_bssbff'


    if ((fl1.eq.fl2).and.(fl3.eq.fl4)) then 

       if (ngluon == 0) then 

          res = czero

          e1 = g_bff(edumm,kdumm,sp(:,1),p(:,1),fl1,&
               &sp(:,2),p(:,2),fl2,0,0,giarray,qiarray(1:2),pol_int)
          k1 = p(:,1) + p(:,2) 
          k1sq = sc(k1,k1)

          e2 = g_bff(edumm,kdumm,sp(:,3),p(:,3),fl3,&
               &sp(:,4),p(:,4),fl4,0,0,giarray,qiarray(3:4),pol_int)
          k2 = p(:,3) + p(:,4)
          k2sq = sc(k2,k2)

          tmp = (-ci/k1sq)*(-ci/k2sq)*vggg(e1,k1,e2,k2)

          res = res + tmp  !# 1


          sp2 = bf_fbff(edumm,kdumm,sp(:,2:4),p(:,2:4),&
               &fll(2:4),fl1,0,0,0,giarray,qiarray(2:4),pol_int)
          k2 = -p(:,2) - p(:,3) - p(:,4)
          k2sq = sc(k2,k2)

          sp2 = spi2(k2,sp2) !+ mass*sp2
          sp2 = ci/k2sq*sp2

          sp1 = sp(:,1)

          tmp = vbqq(Dv,sp1,sp2)

          res = res + tmp       ! # 2

          sp2 = f_bffbf(edumm,kdumm,sp(:,1:3),p(:,1:3),&
               &fll(1:3),fl4,0,0,0,giarray,qiarray(1:3),pol_int)

          k2 = p(:,1) + p(:,2) + p(:,3)
          k2sq = sc(k2,k2)

          sp2 = spb2(sp2,k2) !+ mass*sp2
          sp2 = ci/k2sq*sp2

          sp1 = sp(:,4)

          tmp = vbqq(Dv,sp2,sp1)

          res = res + tmp       ! # 3




       endif   !end condition for ngluon 

    else
       stop ' g_bssbff: flavour combination not implemented' 
    endif !(for a flavor combination ) 



    ! -- store current 
!    if (present(giarray)) call store_result(pol_int,res,giarray,qiarray)

  end function g_bssbff



!$$$    recursive function f_bffbffbf(e,k,sp,p,fll,fl0,&
!$$$         &ng1,ng2,ng3,ng4,ng5,giarray,qiarray,pol_int) result(res)
!$$$      double complex, intent(in) :: e(:,:), k(:,:)
!$$$      double complex, intent(in) :: sp(:,:),p(:,:)
!$$$      integer, intent(in) ::  ng1,ng2,ng3,ng4,ng5
!$$$      character, intent(in) :: fll(:)*3
!$$$      character, intent(in) :: fl0*3   ! flavor off-shell f-line
!$$$      integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int 
!$$$      ! -----------------------------------------------------------------------
!$$$      character :: fl1*3,fl2*3,fl3*3, fl4*3, fl5*3
!$$$      integer :: ngluon, ng6,m1,ngmax
!$$$      integer, parameter :: Ndumm=0
!$$$      double complex             :: res(size(sp,dim=1))
!$$$      double complex             :: tmp(size(sp,dim=1))
!$$$      double complex             :: k1(size(k,dim=1))
!$$$      double complex             :: k2(size(k,dim=1))
!$$$      double complex             :: sp1(size(sp,dim=1))
!$$$      double complex             :: sp2(size(sp,dim=1))
!$$$      double complex             :: e1(size(e,dim=1))
!$$$      double complex             :: e2(size(e,dim=1))
!$$$      double complex             :: kdumm(size(k,dim=1),Ndumm)
!$$$      double complex             :: edumm(size(e,dim=1),Ndumm)
!$$$      double complex             :: k1sq,k2sq
!$$$  !    real(dp)                :: mass
!$$$      logical                   :: done, onshell  
!$$$  
!$$$      if (verbose)    write(*,*) 'entering f_bffbffbf',ng1,ng2,ng3,ng4,ng5
!$$$      done = .false. 
!$$$  !    if (present(giarray)) then 
!$$$  !       if (size(qiarray) /= 5) stop 'f_bffbffbf: wrong size qiarray'
!$$$  !       if (size(e,dim=2) /= size(giarray)) stop 'f_bffbffbf: ng= size(giarray)' 
!$$$  !       call memory_check(pol_int,res,done,giarray,qiarray)
!$$$  !       if (done) return 
!$$$  !    else
!$$$  !       if (i_warn < max_warn) then 
!$$$  !          write(*,*) 'f_bffbffbf: giarray missing', i_warn 
!$$$  !          i_warn = i_warn+1
!$$$  !       endif
!$$$  !    endif
!$$$  
!$$$  !    mass = mt
!$$$  
!$$$      ngluon = size(e,dim=2)
!$$$      ng6 = ngluon - ng1 - ng2-ng3 - ng4 - ng5
!$$$  
!$$$      if (ng6 < 0) write(*,*) 'ERROR IN CURRENT F'
!$$$  
!$$$      fl1 = fll(1)
!$$$      fl2 = fll(2)
!$$$      fl3 = fll(3)
!$$$      fl4 = fll(4)
!$$$      fl5 = fll(5)
!$$$  
!$$$      if ((fl0.eq.fl1).and.(fl2.eq.fl3).and.(fl4.eq.fl5)) then 
!$$$  
!$$$      if (ngluon == 0) then 
!$$$  
!$$$            res = czero
!$$$  
!$$$            if (ferm_loops_Z) then 
!$$$  
!$$$               e1 = g_sbsfbf(edumm,kdumm,sp(:,1:4),p(:,1:4),fll(1:4),&
!$$$                    &0,0,0,0,&
!$$$                    &giarray,qiarray(1:4),pol_int)
!$$$               
!$$$               k1 = p(:,1)+p(:,2)+p(:,3)+p(:,4)
!$$$               k1sq=sc(k1,k1)
!$$$               
!$$$               if (abs(k1sq) > propcut.and.fl0==fl5) then 
!$$$                  tmp = -ci/k1sq*vgbq(e1,sp(:,5))
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$  
!$$$            res = res + tmp  !#1
!$$$  
!$$$  
!$$$            e2 = g_fbf(edumm,kdumm,sp(:,4),p(:,4),fl4,&
!$$$                 &sp(:,5),p(:,5),fl5,0,0,&
!$$$                 &giarray,qiarray(4:5),pol_int)
!$$$            k2 = p(:,4) + p(:,5)
!$$$            k2sq = sc(k2,k2)
!$$$  
!$$$            sp1 = f_bffbf(edumm,kdumm,sp(:,1:3),p(:,1:3),&
!$$$                 &fll(1:3),fl1,0,0,0,giarray,qiarray(1:3),pol_int)
!$$$            k1 = p(:,1) + p(:,2) + p(:,3)
!$$$            k1sq = sc(k1,k1)
!$$$            sp1 = spb2(sp1,k1) !+ mass*sp1
!$$$  
!$$$            tmp = -ci/k2sq*ci/k1sq*vqg(sp1,e2)
!$$$  
!$$$            res = res + tmp  !#2
!$$$  
!$$$            else
!$$$            e2 = g_sbsfbf(edumm,kdumm,sp(:,2:5),p(:,2:5),fll(2:5),&
!$$$                 &0,0,0,0,&
!$$$                 &giarray,qiarray(2:5),pol_int)
!$$$  
!$$$            k2 = p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$            k2sq=sc(k2,k2)
!$$$  
!$$$            if (abs(k2sq) > propcut.and.fl0==fl1) then 
!$$$               tmp = -ci/k2sq*vqg(sp(:,1),e2)
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            res = res + tmp  !#1
!$$$  
!$$$  
!$$$            e2 = g_fbf(edumm,kdumm,sp(:,4),p(:,4),fl4,&
!$$$                 &sp(:,5),p(:,5),fl5,0,0,&
!$$$                 &giarray,qiarray(4:5),pol_int)
!$$$            k2 = p(:,4) + p(:,5)
!$$$            k2sq = sc(k2,k2)
!$$$  
!$$$            sp1 = f_bffbf(edumm,kdumm,sp(:,1:3),p(:,1:3),&
!$$$                 &fll(1:3),fl1,0,0,0,giarray,qiarray(1:3),pol_int)
!$$$            k1 = p(:,1) + p(:,2) + p(:,3)
!$$$            k1sq = sc(k1,k1)
!$$$            sp1 = spb2(sp1,k1) !+ mass*sp1
!$$$  
!$$$            tmp = -ci/k2sq*ci/k1sq*vqg(sp1,e2)
!$$$  
!$$$            res = res + tmp  !#2
!$$$         endif
!$$$  
!$$$         else                          !for #gluons > 0
!$$$  
!$$$  
!$$$            res = czero
!$$$  
!$$$            do m1 = 0, ng2
!$$$  
!$$$               sp1 = f(e(:,1:ng1+m1),k(:,1:ng1+m1),&
!$$$                    &sp(:,1),p(:,1),fl1,fl0,ng1,&
!$$$                    &giarray(1:ng1+m1),qiarray(1:1),pol_int)
!$$$               if (1<=ng1+m1) then 
!$$$                  k1 = sum(k(:,1:ng1+m1),dim=2)
!$$$               else
!$$$                  k1 = czero 
!$$$               endif 
!$$$               k1 = k1 + p(:,1)
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$               if (ng1 > 0.or.m1 > 0) then 
!$$$                  sp1 = spb2(sp1,k1)!+mass*sp1
!$$$                  sp1 = ci/k1sq*sp1
!$$$               endif
!$$$  
!$$$               e2 = g_sbsfbf(e(:,ng1+m1+1:ngluon),&
!$$$                    &k(:,ng1+m1+1:ngluon),sp(:,2:5),p(:,2:5),fll(2:5),&
!$$$                    &ng2-m1,ng3,ng4,ng5,&
!$$$                    &giarray(ng1+m1+1:ngluon),qiarray(2:5),pol_int)
!$$$               if (ng1+m1+1<=ngluon) then 
!$$$                  k2 = sum(k(:,ng1+m1+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$  
!$$$               tmp = -ci/k2sq*vqg(sp1,e2)
!$$$  
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$            enddo ! #1
!$$$  
!$$$            do m1 = 0,ng4
!$$$  
!$$$  
!$$$               sp1 = f_bffbf(e(:,1:ng1+ng2+ng3+m1),&
!$$$                    &k(:,1:ng1+ng2+ng3+m1),sp(:,1:3),p(:,1:3),&
!$$$                    &fll(1:3),fl0,ng1,ng2,ng3,&
!$$$                    &giarray(1:ng1+ng2+ng3+m1),qiarray(1:3),pol_int)
!$$$  
!$$$               if (1<=ng1+ng2+ng3+m1) then 
!$$$                  k1 = sum(k(:,1:ng1+ng2+ng3+m1),dim=2)
!$$$               else
!$$$                  k1 = czero 
!$$$               endif
!$$$               k1 = k1 + p(:,1)+p(:,2)+p(:,3)
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$               sp1 = spb2(sp1,k1) !+ mass*sp1
!$$$  
!$$$               e2 = g_fbf(e(:,ng1+ng2+ng3+m1+1:ngluon),&
!$$$                    &k(:,ng1+ng2+ng3+m1+1:ngluon),sp(:,4),p(:,4),fl4,&
!$$$                    &sp(:,5),p(:,5),fl5,ng4-m1,ng5,&
!$$$                    &giarray(ng1+ng2+ng3+m1+1:ngluon),qiarray(4:5),pol_int)
!$$$  
!$$$               if (ng1+ng2+ng3+m1+1<=ngluon) then 
!$$$                  k2 = sum(k(:,ng1+ng2+ng3+m1+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,4)+p(:,5)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$               tmp = ci/k1sq*(-ci)/k2sq*vqg(sp1,e2)
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$            enddo   ! #2
!$$$  
!$$$  
!$$$            do m1 = 1,ng1
!$$$  
!$$$               e1 = vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)
!$$$               k1= sum(k(:,1:m1),dim=2)
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$               sp2 = f_bffbffbf(e(:,m1+1:ngluon),&
!$$$                    &k(:,m1+1:ngluon),sp,p,fll,fl0,&
!$$$                    &ng1-m1,ng2,ng3,ng4,ng5,&
!$$$                    &giarray(m1+1:ngluon),qiarray,pol_int)
!$$$               if (m1+1<=ngluon) then 
!$$$                  k2= sum(k(:,m1+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,1)+p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$               sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$  
!$$$               sp2 = ci/k2sq*sp2
!$$$  
!$$$  
!$$$  
!$$$               if (m1 > 1) then 
!$$$                  e1 = -ci/k1sq*e1
!$$$               endif
!$$$  
!$$$               tmp  = vgq(e1,sp2)
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$            enddo ! #3
!$$$  
!$$$  
!$$$  
!$$$            do m1 = 1,ng6
!$$$  
!$$$               e1 = vgluon(e(:,ng1+ng2+ng3+ng4+ng5+m1:ngluon),&
!$$$                    &k(:,ng1+ng2+ng3+ng4+ng5+m1:ngluon),&
!$$$                    &giarray(ng1+ng2+ng3+ng4+ng5+m1:ngluon),pol_int)
!$$$               if (ng1+ng2+ng3+ng4+ng5+m1<=ngluon) then 
!$$$                  k1= sum(k(:,ng1+ng2+ng3+ng4+ng5+m1:ngluon),dim=2)
!$$$               else
!$$$                  k1 = czero 
!$$$               endif
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$               sp2 = f_bffbffbf(e(:,1:ng1+ng2+ng3+ng4+ng5+m1-1),&
!$$$                    &k(:,1:ng1+ng2+ng3+ng4+ng5+m1-1),sp,p,fll,fl0,&
!$$$                    &ng1,ng2,ng3,ng4,ng5,&
!$$$                    &giarray(1:ng1+ng2+ng3+ng4+ng5+m1-1),qiarray,pol_int)
!$$$               if (1<=ng1+ng2+ng3+ng4+ng5+m1-1) then 
!$$$                  k2= sum(k(:,1:ng1+ng2+ng3+ng4+ng5+m1-1),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,1)+p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$               sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$  
!$$$               sp2 = ci/k2sq*sp2
!$$$  
!$$$  
!$$$               if (m1 < ng6) then 
!$$$                  e1 = -ci/k1sq*e1
!$$$               endif
!$$$  
!$$$               tmp  = vqg(sp2,e1)
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$            enddo ! #4
!$$$  
!$$$         endif
!$$$  
!$$$      endif
!$$$  
!$$$  
!$$$  
!$$$      if ((fl0.eq.fl3).and.(fl1.eq.fl2)&
!$$$           &.and.(fl4.eq.fl5)) then 
!$$$  
!$$$  
!$$$         if (ngluon == 0) then 
!$$$  
!$$$            res = czero
!$$$  
!$$$            e1 = g_bff(edumm,kdumm,sp(:,1),p(:,1),fl1,&
!$$$                 &sp(:,2),p(:,2),fl2,0,0,&
!$$$                 &giarray,qiarray(1:2),pol_int)
!$$$            k1 = p(:,1)+p(:,2)
!$$$            k1sq = sc(k1,k1)
!$$$  
!$$$            sp2 = f_bffbf(edumm,kdumm,sp(:,3:5),p(:,3:5),&
!$$$                 &fll(3:5),fl0,0,0,0,&
!$$$                 &giarray,qiarray(3:5),pol_int)
!$$$            k2 = p(:,3)+p(:,4)+p(:,5)
!$$$            k2sq = sc(k2,k2)
!$$$            sp2 = spb2(sp2,k2)!+mass*sp2
!$$$  
!$$$            if (abs(k2sq) > propcut.and.fl3==fl0) then 
!$$$               tmp = -ci/k1sq*ci/k2sq*vgq(e1,sp2)
!$$$            else 
!$$$               tmp = czero 
!$$$            endif
!$$$  
!$$$            res = res + tmp  !#1
!$$$  
!$$$  
!$$$            e2 = g_fbf(edumm,kdumm,sp(:,4),p(:,4),fl4,&
!$$$                 &sp(:,5),p(:,5),fl5,0,0,&
!$$$                 &giarray,qiarray(4:5),pol_int)
!$$$            k2 = p(:,4) + p(:,5)
!$$$            k2sq = sc(k2,k2)
!$$$  
!$$$            sp1 = f_bffbf(edumm,kdumm,sp(:,1:3),p(:,1:3),&
!$$$                 &fll(1:3),fl0,0,0,0,&
!$$$                 &giarray,qiarray(1:3),pol_int)
!$$$            k1 = p(:,1) + p(:,2) + p(:,3)
!$$$            k1sq = sc(k1,k1)
!$$$            sp1 = spb2(sp1,k1) !+ mass*sp1
!$$$  
!$$$            tmp = -ci/k2sq*ci/k1sq*vqg(sp1,e2)
!$$$  
!$$$            res = res + tmp  !#2
!$$$  
!$$$         else       !for #gluons > 0
!$$$  
!$$$  
!$$$  
!$$$            res = czero
!$$$  
!$$$            do m1 = 0, ng3
!$$$  
!$$$               e1 = g_bff(e(:,1:ng1+ng2+m1),&
!$$$                    &k(:,1:ng1+ng2+m1),sp(:,1),p(:,1),fl1,&
!$$$                    &sp(:,2),p(:,2),fl2,&
!$$$                    &ng1,ng2,giarray(1:ng1+ng2+m1),qiarray(1:2),pol_int)
!$$$               if (1<=ng1+ng2+m1) then 
!$$$                  k1 = sum(k(:,1:ng1+ng2+m1),dim=2)
!$$$               else
!$$$                  k1 = czero 
!$$$               endif
!$$$               k1 = k1 + p(:,1)+p(:,2)
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$  
!$$$               sp2 = f_bffbf(e(:,ng1+ng2+m1+1:ngluon),&
!$$$                    &k(:,ng1+ng2+m1+1:ngluon),&
!$$$                    &sp(:,3:5),p(:,3:5),fll(3:5),fl0,ng3-m1,ng4,ng5,&
!$$$                    &giarray(ng1+ng2+m1+1:ngluon),qiarray(3:5),pol_int)
!$$$               if (ng1+ng2+m1+1<=ngluon) then 
!$$$                  k2 = sum(k(:,ng1+ng2+m1+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,3)+p(:,4)+p(:,5)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$  
!$$$               sp2 = spb2(sp2,k2)!+mass*sp2
!$$$               sp2 = ci/k2sq*sp2
!$$$  
!$$$  
!$$$               tmp = -ci/k1sq*vgq(e1,sp2)
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$            enddo ! #1
!$$$  
!$$$            do m1 = 0,ng4
!$$$  
!$$$  
!$$$               sp1 = f_bffbf(e(:,1:ng1+ng2+ng3+m1),&
!$$$                    &k(:,1:ng1+ng2+ng3+m1),sp(:,1:3),p(:,1:3),&
!$$$                    &fll(1:3),fl0,ng1,ng2,ng3,&
!$$$                    &giarray(1:ng1+ng2+ng3+m1),qiarray(1:3),pol_int)
!$$$  
!$$$               if (1<=ng1+ng2+ng3+m1) then 
!$$$                  k1 = sum(k(:,1:ng1+ng2+ng3+m1),dim=2)
!$$$               else
!$$$                  k1 = czero 
!$$$               endif
!$$$               k1 = k1 + p(:,1)+p(:,2)+p(:,3)
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$               sp1 = spb2(sp1,k1) !+ mass*sp1
!$$$  
!$$$               e2 = g_fbf(e(:,ng1+ng2+ng3+m1+1:ngluon),&
!$$$                    &k(:,ng1+ng2+ng3+m1+1:ngluon),sp(:,4),p(:,4),fl4,&
!$$$                    &sp(:,5),p(:,5),fl5,ng4-m1,ng5,&
!$$$                    &giarray(ng1+ng2+ng3+m1+1:ngluon),qiarray(4:5),pol_int)
!$$$  
!$$$               if (ng1+ng2+ng3+m1+1<=ngluon) then 
!$$$                  k2 = sum(k(:,ng1+ng2+ng3+m1+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,4)+p(:,5)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$               tmp = ci/k1sq*(-ci)/k2sq*vqg(sp1,e2)
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$            enddo   ! #2  
!$$$  
!$$$  
!$$$            do m1 = 1,ng1
!$$$  
!$$$               e1 = vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)
!$$$               k1= sum(k(:,1:m1),dim=2)
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$               sp2 = f_bffbffbf(e(:,m1+1:ngluon),&
!$$$                    &k(:,m1+1:ngluon),sp,p,fll,fl0,&
!$$$                    &ng1-m1,ng2,ng3,ng4,ng5,&
!$$$                    &giarray(m1+1:ngluon),qiarray,pol_int)
!$$$               if (m1+1<=ngluon) then 
!$$$                  k2= sum(k(:,m1+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,1)+p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$               sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$  
!$$$               sp2 = ci/k2sq*sp2
!$$$  
!$$$  
!$$$  
!$$$               if (m1 > 1) then 
!$$$                  e1 = -ci/k1sq*e1
!$$$               endif
!$$$  
!$$$               tmp  = vgq(e1,sp2)
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$            enddo ! #3
!$$$  
!$$$  
!$$$  
!$$$            do m1 = 1,ng6
!$$$  
!$$$               e1 = vgluon(e(:,ng1+ng2+ng3+ng4+ng5+m1:ngluon),&
!$$$                    &k(:,ng1+ng2+ng3+ng4+ng5+m1:ngluon),&
!$$$                    &giarray(ng1+ng2+ng3+ng4+ng5+m1:ngluon),pol_int)
!$$$               if (ng1+ng2+ng3+ng4+ng5+m1<=ngluon) then 
!$$$                  k1= sum(k(:,ng1+ng2+ng3+ng4+ng5+m1:ngluon),dim=2)
!$$$               else
!$$$                  k1 = czero 
!$$$               endif
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$               sp2 = f_bffbffbf(e(:,1:ng1+ng2+ng3+ng4+ng5+m1-1),&
!$$$                    &k(:,1:ng1+ng2+ng3+ng4+ng5+m1-1),sp,p,fll,fl0,&
!$$$                    &ng1,ng2,ng3,ng4,ng5,&
!$$$                    &giarray(1:ng1+ng2+ng3+ng4+ng5+m1-1),qiarray,pol_int)
!$$$               if (1<=ng1+ng2+ng3+ng4+ng5+m1-1) then 
!$$$                  k2= sum(k(:,1:ng1+ng2+ng3+ng4+ng5+m1-1),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,1)+p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$               sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$  
!$$$               sp2 = ci/k2sq*sp2
!$$$  
!$$$  
!$$$               if (m1 < ng6) then 
!$$$                  e1 = -ci/k1sq*e1
!$$$               endif
!$$$  
!$$$               tmp  = vqg(sp2,e1)
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$            enddo ! #4
!$$$  
!$$$  
!$$$  
!$$$         endif
!$$$  
!$$$      endif  ! for second flavor           
!$$$  
!$$$  
!$$$      if ((fl0.eq.fl1).and.(fl2.eq.fl5)&
!$$$           &.and.(fl3.eq.fl4)) then 
!$$$  
!$$$  
!$$$  
!$$$         if (ngluon == 0) then 
!$$$  
!$$$            res = czero
!$$$  
!$$$  
!$$$            e1 = g_sbsfbf(edumm,kdumm,sp(:,2:5),p(:,2:5),fll(2:5),&
!$$$                 &0,0,0,0,giarray,qiarray(2:5),pol_int)
!$$$            k1 = p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$            k1sq = sc(k1,k1)
!$$$  
!$$$  
!$$$            tmp = -ci/k1sq*vqg(sp(:,1),e1)
!$$$  
!$$$  
!$$$            res = res + tmp  !#1
!$$$  
!$$$         else   !for #gluons > 0
!$$$  
!$$$  
!$$$            res = czero
!$$$  
!$$$            do m1 = 0, ng2
!$$$  
!$$$               sp1 = f(e(:,1:ng1+m1),k(:,1:ng1+m1),&
!$$$                    &sp(:,1),p(:,1),fl1,fl0,ng1,&
!$$$                    &giarray(1:ng1+m1),qiarray(1:1),pol_int)
!$$$               if (1<=ng1+m1) then 
!$$$                  k1 = sum(k(:,1:ng1+m1),dim=2)
!$$$               else
!$$$                  k1 = czero 
!$$$               endif
!$$$               k1 = k1 + p(:,1)
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$               if (ng1 > 0.or.m1 > 0) then 
!$$$                  sp1 = spb2(sp1,k1)!+mass*sp1
!$$$                  sp1 = ci/k1sq*sp1
!$$$               endif
!$$$  
!$$$               e2 = g_sbsfbf(e(:,ng1+m1+1:ngluon),&
!$$$                    &k(:,ng1+m1+1:ngluon),sp(:,2:5),p(:,2:5),fll(2:5),&
!$$$                    &ng2-m1,ng3,ng4,ng5,&
!$$$                    &giarray(ng1+m1+1:ngluon),qiarray(2:5),pol_int)
!$$$               if (ng1+m1+1<=ngluon) then 
!$$$                  k2 = sum(k(:,ng1+m1+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$  
!$$$               tmp = -ci/k2sq*vqg(sp1,e2)
!$$$  
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$            enddo ! #1
!$$$  
!$$$  
!$$$  
!$$$            do m1 = 1,ng1
!$$$  
!$$$               e1 = vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)
!$$$               k1= sum(k(:,1:m1),dim=2)
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$               sp2 = f_bffbffbf(e(:,m1+1:ngluon),&
!$$$                    &k(:,m1+1:ngluon),sp,p,fll,fl0,&
!$$$                    &ng1-m1,ng2,ng3,ng4,ng5,&
!$$$                    &giarray(m1+1:ngluon),qiarray,pol_int)
!$$$               if (m1+1<=ngluon) then 
!$$$                  k2= sum(k(:,m1+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,1)+p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$               sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$  
!$$$               sp2 = ci/k2sq*sp2
!$$$  
!$$$  
!$$$  
!$$$               if (m1 > 1) then 
!$$$                  e1 = -ci/k1sq*e1
!$$$               endif
!$$$  
!$$$               tmp  = vgq(e1,sp2)
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$            enddo ! #3
!$$$  
!$$$  
!$$$  
!$$$            do m1 = 1,ng6
!$$$  
!$$$               e1 = vgluon(e(:,ng1+ng2+ng3+ng4+ng5+m1:ngluon),&
!$$$                    &k(:,ng1+ng2+ng3+ng4+ng5+m1:ngluon),&
!$$$                    &giarray(ng1+ng2+ng3+ng4+ng5+m1:ngluon),pol_int)
!$$$               if (ng1+ng2+ng3+ng4+ng5+m1<=ngluon) then 
!$$$                  k1= sum(k(:,ng1+ng2+ng3+ng4+ng5+m1:ngluon),dim=2)
!$$$               else
!$$$                  k1 = czero 
!$$$               endif
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$               sp2 = f_bffbffbf(e(:,1:ng1+ng2+ng3+ng4+ng5+m1-1),&
!$$$                    &k(:,1:ng1+ng2+ng3+ng4+ng5+m1-1),sp,p,fll,fl0,&
!$$$                    &ng1,ng2,ng3,ng4,ng5,&
!$$$                    &giarray(1:ng1+ng2+ng3+ng4+ng5+m1-1),qiarray,pol_int)
!$$$               if (1<=ng1+ng2+ng3+ng4+ng5+m1-1) then 
!$$$                  k2= sum(k(:,1:ng1+ng2+ng3+ng4+ng5+m1-1),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,1)+p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$               sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$  
!$$$               sp2 = ci/k2sq*sp2
!$$$  
!$$$  
!$$$               if (m1 < ng6) then 
!$$$                  e1 = -ci/k1sq*e1
!$$$               endif
!$$$  
!$$$               tmp  = vqg(sp2,e1)
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$            enddo ! #4
!$$$  
!$$$  
!$$$         endif
!$$$  
!$$$      endif  ! for second flavor           
!$$$  
!$$$  
!$$$  
!$$$      ! --------------------------------------------------------------------------      
!$$$      ! --- NEW FOR FERMION_LOOPS_Z_SBS
!$$$      if ((fl0.eq.fl5).and.(fl1.eq.fl2).and.(fl3.eq.fl4)) then 
!$$$  
!$$$  
!$$$      if ( abs(pmass(p(:,1)+p(:,2)+p(:,3)+p(:,4)+sum(k,dim=2)) - mz) < tol) then 
!$$$         onshell =.true. 
!$$$      else
!$$$         onshell = .false.
!$$$      endif
!$$$  
!$$$      if (ngluon == 0) then 
!$$$  
!$$$            res = czero
!$$$  
!$$$            if (onshell) then 
!$$$               tmp = czero 
!$$$            else
!$$$               e1 = g_bssbff(edumm,kdumm,sp(:,1:4),p(:,1:4),fll(1:4),&
!$$$                    &0,0,0,0,&
!$$$                    &giarray,qiarray(1:4),pol_int)
!$$$               
!$$$               k1 = p(:,1)+p(:,2)+p(:,3)+p(:,4)
!$$$               k1sq=sc(k1,k1)
!$$$               
!$$$               if (abs(k1sq) > propcut.and.fl0==fl5) then 
!$$$                  tmp = -ci/k1sq*vgq(e1,sp(:,5))
!$$$               else 
!$$$                  tmp = czero 
!$$$               endif
!$$$            endif
!$$$            res = res + tmp  !#1
!$$$            
!$$$            e1 = g_bff(edumm,kdumm,sp(:,1),p(:,1),fl1,&
!$$$                 &sp(:,2),p(:,2),fl2,0,0,&
!$$$                 &giarray,qiarray(1:2),pol_int)
!$$$            k1 = p(:,1) + p(:,2)
!$$$            k1sq = sc(k1,k1)
!$$$  
!$$$            sp2 = f_bffbf(edumm,kdumm,sp(:,3:5),p(:,3:5),&
!$$$                 &fll(3:5),fl0,0,0,0,giarray,qiarray(3:5),pol_int)
!$$$            k2 = p(:,3) + p(:,4) + p(:,5)
!$$$            k2sq = sc(k2,k2)
!$$$            sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$  
!$$$            tmp = -ci/k2sq*ci/k1sq*vgq(e1,sp2)
!$$$  
!$$$            res = res + tmp  !#2
!$$$  
!$$$  
!$$$         else                          !for #gluons > 0
!$$$  
!$$$  
!$$$            res = czero
!$$$  
!$$$            do m1 = 0, ng3
!$$$  
!$$$               sp1 = f_bffbf(e(:,ng1+ng2+m1+1:ngluon),k(:,ng1+ng2+m1+1:ngluon),&
!$$$                    &sp(:,3:5),p(:,3:5),fll(3:5),fl0,ng3-m1,ng4,ng5,&
!$$$                    &giarray(ng1+ng2+m1+1:ngluon),qiarray(3:5),pol_int)
!$$$               k1 = sum(k(:,ng1+ng2+m1+1:ngluon),dim=2)
!$$$               k1 = p(:,3) +p(:,4)+p(:,5)+k1
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$               sp1 = spb2(sp1,k1)!+mass*sp1
!$$$               sp1 = ci/k1sq*sp1
!$$$  
!$$$               e2 = g_bff(e(:,1:ng1+ng2+m1),&
!$$$                    &k(:,1:ng1+ng2+m1),sp(:,1),p(:,1),fll(1),sp(:,2),p(:,2),fll(2),&
!$$$                    &ng1,ng2,&
!$$$                    &giarray(1:ng1+ng2+m1),qiarray(1:2),pol_int)
!$$$               k2 = sum(k(:,1:ng1+ng2+m1),dim=2)
!$$$               k2 = k2 + p(:,1)+p(:,2)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$               tmp = -ci/k2sq*vgq(e2,sp1)
!$$$               
!$$$               res = res + tmp
!$$$  
!$$$            enddo ! #1
!$$$  
!$$$  
!$$$  
!$$$            if (onshell .and. ng6 == 0) then 
!$$$               ngmax = ng5-1 
!$$$            else
!$$$               ngmax = ng5 
!$$$            endif
!$$$  
!$$$            do m1 = 0,ngmax
!$$$  
!$$$               
!$$$               e1 = g_bssbff(e(:,1:ng1+ng2+ng3+ng4+m1),&
!$$$                    &k(:,1:ng1+ng2+ng3+ng4+m1),sp(:,1:4),p(:,1:4),&
!$$$                    &fll(1:4),ng1,ng2,ng3,ng4,&
!$$$                    &giarray(1:ng1+ng2+ng3+ng4+m1),qiarray(1:4),pol_int)
!$$$  
!$$$               if (1<=ng1+ng2+ng3+ng4+m1) then 
!$$$                  k1 = sum(k(:,1:ng1+ng2+ng3+ng4+m1),dim=2)
!$$$               else
!$$$                  k1 = czero 
!$$$               endif
!$$$               k1 = k1 + p(:,1)+p(:,2)+p(:,3)+p(:,4)
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$  
!$$$               if (ng1+ng2+ng3+ng4+m1+1<=ngluon) then 
!$$$                  sp2 = f(e(:,ng1+ng2+ng3+ng4+m1+1:ngluon),&
!$$$                       &k(:,ng1+ng2+ng3+ng4+m1+1:ngluon),sp(:,5),p(:,5),fl5,&
!$$$                       &fl0,ng5-m1,&
!$$$                       &giarray(ng1+ng2+ng3+ng4+m1+1:ngluon),qiarray(5:5),pol_int)
!$$$  
!$$$                  k2 = sum(k(:,ng1+ng2+ng3+ng4+m1+1:ngluon),dim=2)
!$$$                  k2 = k2 + p(:,5)
!$$$                  k2sq = sc(k2,k2)
!$$$                  sp2 = spb2(sp2,k2) !+ mass * sp2 
!$$$                  tmp = (-ci)/k1sq*ci/k2sq*vgq(e1,sp2)
!$$$               else
!$$$                  sp2 = sp(:,5) 
!$$$                  tmp = -ci/k1sq*vgq(e1,sp2)
!$$$               endif
!$$$  
!$$$               res = res + tmp 
!$$$               
!$$$            enddo   ! #2
!$$$  
!$$$  
!$$$            do m1 = 1,ng1
!$$$  
!$$$               e1 = vgluon(e(:,1:m1),k(:,1:m1),giarray(1:m1),pol_int)
!$$$               k1= sum(k(:,1:m1),dim=2)
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$               sp2 = f_bffbffbf(e(:,m1+1:ngluon),&
!$$$                    &k(:,m1+1:ngluon),sp,p,fll,fl0,&
!$$$                    &ng1-m1,ng2,ng3,ng4,ng5,&
!$$$                    &giarray(m1+1:ngluon),qiarray,pol_int)
!$$$               if (m1+1<=ngluon) then 
!$$$                  k2= sum(k(:,m1+1:ngluon),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,1)+p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$               sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$  
!$$$               sp2 = ci/k2sq*sp2
!$$$  
!$$$  
!$$$  
!$$$               if (m1 > 1) then 
!$$$                  e1 = -ci/k1sq*e1
!$$$               endif
!$$$  
!$$$               tmp  = vgq(e1,sp2)
!$$$  
!$$$               ! 
!$$$               res = res + tmp
!$$$  
!$$$            enddo ! #3
!$$$  
!$$$  
!$$$  
!$$$            do m1 = 1,ng6
!$$$  
!$$$               e1 = vgluon(e(:,ng1+ng2+ng3+ng4+ng5+m1:ngluon),&
!$$$                    &k(:,ng1+ng2+ng3+ng4+ng5+m1:ngluon),&
!$$$                    &giarray(ng1+ng2+ng3+ng4+ng5+m1:ngluon),pol_int)
!$$$               if (ng1+ng2+ng3+ng4+ng5+m1<=ngluon) then 
!$$$                  k1= sum(k(:,ng1+ng2+ng3+ng4+ng5+m1:ngluon),dim=2)
!$$$               else
!$$$                  k1 = czero 
!$$$               endif
!$$$               k1sq = sc(k1,k1)
!$$$  
!$$$               sp2 = f_bffbffbf(e(:,1:ng1+ng2+ng3+ng4+ng5+m1-1),&
!$$$                    &k(:,1:ng1+ng2+ng3+ng4+ng5+m1-1),sp,p,fll,fl0,&
!$$$                    &ng1,ng2,ng3,ng4,ng5,&
!$$$                    &giarray(1:ng1+ng2+ng3+ng4+ng5+m1-1),qiarray,pol_int)
!$$$               if (1<=ng1+ng2+ng3+ng4+ng5+m1-1) then 
!$$$                  k2= sum(k(:,1:ng1+ng2+ng3+ng4+ng5+m1-1),dim=2)
!$$$               else
!$$$                  k2 = czero 
!$$$               endif
!$$$               k2 = k2 + p(:,1)+p(:,2)+p(:,3)+p(:,4)+p(:,5)
!$$$               k2sq = sc(k2,k2)
!$$$  
!$$$               sp2 = spb2(sp2,k2) !+ mass*sp2
!$$$  
!$$$               sp2 = ci/k2sq*sp2
!$$$  
!$$$  
!$$$               if (m1 < ng6) then 
!$$$                  e1 = -ci/k1sq*e1
!$$$               endif
!$$$  
!$$$               tmp  = vqg(sp2,e1)
!$$$  
!$$$               res = res + tmp
!$$$  
!$$$            enddo ! #4
!$$$  
!$$$         endif
!$$$  
!$$$      endif
!$$$  
!$$$      ! -- store current 
!$$$      if (present(giarray)) call store_result(pol_int,res,giarray,qiarray)
!$$$  
!$$$  
!$$$    end function f_bffbffbf
!$$$  
!$$$  
!$$$  
!$$$
!$$$ ! -----------------for Z
!$$$ !
!$$$   recursive function f_fbfbf_3(e,k,sp,p,fll,fl0,ng1,ng2,ng3,&
!$$$        &giarray,qiarray,pol_int) result(res)
!$$$     implicit none
!$$$     double complex, intent(in) :: e(:,:), k(:,:)
!$$$     double complex, intent(in) :: sp(:,:),p(:,:)
!$$$     integer, intent(in) ::  ng1,ng2,ng3
!$$$     character, intent(in) :: fll(:)*3
!$$$     character, intent(in) :: fl0*3   ! flavor off-shell f-line
!$$$     integer, intent(in), optional       :: giarray(:),qiarray(:),pol_int 
!$$$     ! -----------------------------------------------------------------------
!$$$     character :: fl1*3,fl2*3,fl3*3
!$$$     integer :: ngluon, ng4, ngL,m,nmax 
!$$$     integer, parameter :: Ndumm=0
!$$$     double complex             :: res(size(sp,dim=1))
!$$$     double complex             :: res_stored(size(sp,dim=1))
!$$$     double complex             :: tmp(size(sp,dim=1))
!$$$     double complex             :: k1(size(k,dim=1))
!$$$     double complex             :: k2(size(k,dim=1))
!$$$     double complex             :: sp1(size(sp,dim=1))
!$$$     double complex             :: sp2(size(sp,dim=1))
!$$$     double complex             :: e1(size(e,dim=1))
!$$$     double complex             :: e2(size(e,dim=1))
!$$$     double complex             :: kdumm(size(k,dim=1),Ndumm)
!$$$     double complex             :: edumm(size(e,dim=1),Ndumm)
!$$$     double complex             :: k1sq,k2sq
!$$$ !    real(dp)                :: mass,mass2
!$$$     logical                   :: done, onshell 
!$$$ 
!$$$     if (verbose) write(*,*) 'entering f_bffbf_3',ng1,ng2,ng3
!$$$     if (verbose .and. present(giarray)) write(*,*) 'entering f_bffbf_3:g',giarray
!$$$     if (verbose .and. present(giarray)) write(*,*) 'entering f_bffbf_3:q',qiarray
!$$$     done = .false. 
!$$$ !    if (present(giarray)) then 
!$$$ !       if (size(qiarray) /= 3) stop 'f_bffbf_3: wrong size qiarray'
!$$$ !       if (size(e,dim=2) /= size(giarray)) stop 'bf_fbf: ng= size(giarray)' 
!$$$ !       call memory_check(pol_int,res,done,giarray,qiarray)
!$$$ !       if (done) return 
!$$$ !!       if (done) res_stored = res 
!$$$ !   else
!$$$ !       if (i_warn < max_warn) then 
!$$$ !          write(*,*) 'f_bffbf_3: giarray missing', i_warn 
!$$$ !          i_warn = i_warn+1
!$$$ !       endif
!$$$ !    endif
!$$$ 
!$$$ !    mass = mt
!$$$ !    mass2 = mt**2
!$$$ 
!$$$     ngluon = size(e,dim=2)
!$$$     ng4 = ngluon - ng1 - ng2-ng3
!$$$     if ((ng1 < 0) .or. (ng2 < 0) .or. (ng3 < 0) .or. (ng4 < 0)) then 
!$$$        write(*,*) 'f_bffbf_3: some ng <0',ng1,ng2,ng3,ng4
!$$$        stop 
!$$$     endif
!$$$     !if (verbose) write(*,*) 'in function f_bffbf_3', ng1, ngluon
!$$$ 
!$$$     fl1 = fll(1)
!$$$     fl2 = fll(2)
!$$$     fl3 = fll(3)
!$$$ 
!$$$ 
!$$$     if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C:f_fbfbf_3'
!$$$ 
!$$$     if ( abs(pmass(p(:,1)+p(:,2)+sum(k,dim=2)) - mz) < tol) then 
!$$$        onshell =.true. 
!$$$     else
!$$$        onshell = .false.
!$$$     endif
!$$$     if (ngluon == 0) then 
!$$$ 
!$$$        res = czero
!$$$ 
!$$$        if (onshell) return 
!$$$        e2 = g_fbf(edumm,kdumm,sp(:,1),p(:,1),fl1,&
!$$$             &sp(:,2),p(:,2),fl2,0,0,giarray,qiarray(1:2),pol_int)
!$$$        k2 = p(:,1)+p(:,2) 
!$$$        k2sq = sc(k2,k2)
!$$$ 
!$$$ 
!$$$        if (abs(k2sq) > propcut.and.fl0==fl3) then 
!$$$           tmp = -ci/k2sq*vgq(e2,sp(:,3))
!$$$        else 
!$$$           tmp = czero 
!$$$        endif
!$$$ 
!$$$        res = res + tmp
!$$$ 
!$$$     else                          !for #gluons > 0
!$$$ 
!$$$        res = czero
!$$$ 
!$$$ 
!$$$        do m=0,ng4-1
!$$$ 
!$$$           ngL = ng1+ ng2+ng3+m      
!$$$ 
!$$$           sp1=f_fbfbf_3(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0,&
!$$$                &ng1,ng2,ng3,giarray(1:ngL),qiarray,pol_int)
!$$$           if (1<=ngL) then 
!$$$              k1=sum(k(:,1:ngL),dim=2)
!$$$           else
!$$$              k1 = czero 
!$$$           end if
!$$$           k1 = k1 + p(:,1)+p(:,2)+p(:,3)
!$$$           k1sq = sc(k1,k1) !- mass2
!$$$ 
!$$$           sp1 = spb2(sp1,k1) !+ mass*sp1
!$$$ 
!$$$           e2 = vgluon(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
!$$$                &giarray(ngL+1:ngluon),pol_int)
!$$$           if (ngL+1<=ngluon) then 
!$$$              k2=sum(k(:,ngL+1:ngluon),dim=2)
!$$$           else
!$$$              k2 = czero 
!$$$           endif
!$$$           k2sq=sc(k2,k2)
!$$$ 
!$$$           if (abs(k1sq) > propcut) then 
!$$$              tmp = ci/k1sq*vqg(sp1,e2)
!$$$           else 
!$$$              tmp = czero 
!$$$           endif
!$$$ 
!$$$           if (m < ng4-1) then 
!$$$              if (abs(k2sq) > propcut) then 
!$$$                 tmp = -ci/k2sq*tmp
!$$$              else 
!$$$                 tmp = czero 
!$$$              endif
!$$$           endif
!$$$ 
!$$$           res = res + tmp
!$$$ !          write(*,*) 'tmp1',(tmp)
!$$$        enddo  !#2
!$$$ 
!$$$ 
!$$$        do m=1,ng1
!$$$ 
!$$$           e1 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)
!$$$           k1=sum(k(:,1:m),dim=2)
!$$$           k1sq=sc(k1,k1)
!$$$ 
!$$$           sp2=f_fbfbf_3(e(:,m+1:ngluon),k(:,m+1:ngluon),&
!$$$                &sp,p,fll,fl0,ng1-m,ng2,ng3,&
!$$$                &giarray(m+1:ngluon),qiarray,pol_int)
!$$$           if (m+1<=ngluon) then 
!$$$              k2=sum(k(:,m+1:ngluon),dim=2)
!$$$           else
!$$$              k2 = czero 
!$$$           endif
!$$$           k2 = k2 + p(:,1) + p(:,2) + p(:,3)
!$$$           k2sq = sc(k2,k2) !- mass2
!$$$           sp2 = spb2(sp2,k2) !+mass*sp2
!$$$ 
!$$$           if (abs(k2sq) > propcut) then 
!$$$              tmp = ci/k2sq*vgq(e1,sp2)
!$$$           else 
!$$$              tmp = czero 
!$$$           endif
!$$$ 
!$$$ 
!$$$           if (m > 1) then 
!$$$              if (abs(k1sq) > propcut) then 
!$$$                 tmp = -ci/k1sq*tmp
!$$$              else 
!$$$                 tmp = czero 
!$$$              endif
!$$$           endif
!$$$ 
!$$$           res = res + tmp
!$$$ !          write(*,*) 'tmp2',(tmp)
!$$$ 
!$$$        enddo  !#3
!$$$ 
!$$$ 
!$$$        if (onshell) then 
!$$$           nmax = ng3-1
!$$$        else
!$$$           nmax = ng3
!$$$        endif
!$$$ 
!$$$        do m=0,nmax ! ng3
!$$$ 
!$$$           ngL = ng1+ ng2+m      
!$$$ 
!$$$           e1 = g_fbf(e(:,1:ngL),k(:,1:ngL),sp(:,1),p(:,1),fl1,&
!$$$                &sp(:,2),p(:,2),fl2,ng1,ng2,giarray(1:ngL),qiarray(1:2),pol_int)
!$$$           if (1<=ngL) then 
!$$$              k1=sum(k(:,1:ngL),dim=2)+p(:,1)+p(:,2)
!$$$           else
!$$$              k1 = p(:,1)+p(:,2)
!$$$           endif
!$$$           k1sq=sc(k1,k1)
!$$$ 
!$$$ 
!$$$           sp2=f(e(:,ngL+1:ngluon),k(:,ngL+1:ngluon),&
!$$$                &sp(:,3),p(:,3),fl3,fl0,ng3-m,&
!$$$                &giarray(ngL+1:ngluon),qiarray(3:3),pol_int)
!$$$           if (ngL+1<=ngluon) then 
!$$$              k2=sum(k(:,ngL+1:ngluon),dim=2)
!$$$           else
!$$$              k2 = czero 
!$$$           endif
!$$$           k2 = k2 + p(:,3)
!$$$           k2sq = sc(k2,k2) !- mass2
!$$$ 
!$$$           if (ng4 > 0.or. m < ng3) then 
!$$$              sp2 = spb2(sp2,k2)!+mass*sp2
!$$$           endif
!$$$ 
!$$$ 
!$$$           if (abs(k1sq) > propcut) then 
!$$$              tmp = -ci/k1sq*vgq(e1,sp2)
!$$$           else 
!$$$              tmp = czero 
!$$$           endif
!$$$ 
!$$$ 
!$$$ 
!$$$           if (ng4 > 0.or.m < ng3) then 
!$$$              if (abs(k2sq) > propcut) then 
!$$$                 tmp = ci/k2sq*tmp
!$$$              else 
!$$$                 tmp = czero 
!$$$              endif
!$$$           endif
!$$$ 
!$$$           res = res + tmp
!$$$ 
!$$$        enddo  !#4
!$$$ 
!$$$ 
!$$$ 
!$$$     endif
!$$$ 
!$$$     !! -- store current 
!$$$     if (present(giarray)) call store_result(pol_int,res,giarray,qiarray)
!$$$ 
!$$$ 
!$$$   end function f_fbfbf_3



end module recurrenceB
