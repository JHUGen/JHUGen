  !------------------------------------------------------------------------!
  ! Authors: Tom Melia, Kirill Melnikov, Raoul Rontsch, Giulia Zanderighi  !
  ! Date: 25/10/2010                                                       !
  ! Used for arXiv:1007.5313 Wp Wp 2 jets                                  !
  !------------------------------------------------------------------------!
  module recurrence
    use consts_dp
    use recurrenceA
    use recurrenceB
    use recurrenceC
    use spinfns
    implicit none
    public :: fWW_bffbf
 
   contains




  recursive function fWW_bffbf(e,k,sp,p,fll,fl0,&
       &eWin,kWin,ng1,ng2,ng3,sw,giarray,qiarray,WWid,pol_int) result(res)
    double complex, intent(in) :: e(:,:), k(:,:)
    double complex, intent(in) :: sp(:,:),p(:,:)
    double complex, intent(in) :: eWin(:,:),kWin(:,:)
    integer, intent(in) ::  ng1,ng2,ng3,sw
    character, intent(in) :: fll(:)*3
    character, intent(in) :: fl0*3   ! flavor off-shell f-line
    integer, intent(in), optional       :: giarray(:),qiarray(:),WWid(:),pol_int   ! --------------------------------------------------------------------
    double complex::sp4(size(sp,dim=1)), res(size(sp,dim=1)),tmp(size(sp,dim=1))
    double complex:: e2(size(e,dim=1)),sp1(size(sp,dim=1)),sp2(size(sp,dim=1))
    double complex::kW(size(kWin,dim=1), size(kWin,dim=2)), eW(size(eWin,dim=1), size(eWin,dim=2))
    character::flaux*3,fl1*3,fl2*3,fl3*3, fllnew(3)*3
    integer::ngluons,ng4,m,ngL
    integer, parameter::Ndumm = 0
    double complex             :: kdumm(size(k,dim=1),Ndumm),k5sq
    double complex             :: edumm(size(e,dim=1),Ndumm)
    double complex::k2sq,k4sq,k1sq,k4(size(k,dim=1)), k2(size(k,dim=1)), k5(size(k,dim=1))
    double complex::k1(size(k,dim=1)), spnew(size(sp,dim=1),size(sp,dim=2)),pnew(size(p,dim=1),size(p,dim=2))
    !double precision :: mass,mass2

!    mass = mt
!    mass2 = mass**2

    fl1 = fll(1)
    fl2 = fll(2)
    fl3 = fll(3)
    kW(:,1) = kWin(:,1)
    kW(:,2) = kWin(:,2)

    eW(:,1) = eWin(:,1)
    eW(:,2) = eWin(:,2)


    ngluons = size(e,dim=2)


    ng4 = ngluons - ng1 - ng2-ng3
    res = czero
    if (ng4 < 0) write(*,*) 'ERROR IN CURRENT C: fWW_bffbf', ngluons, ng1, ng2, ng3, sw
    if (ngluons == 0) then

!       if (sw == 1) then
!
!          !--- Both W's above gluon
!
!          sp4  = fWW(edumm, kdumm, sp(:,1),p(:,1), fl1, fl0,eW, kW, 0, giarray,&
!               &qiarray(1:1), WWid, pol_int)
!          k4 = p(:,1)+kW(:,1)+kW(:,2)
!          k4sq = sc(k4,k4)
!          sp4 = spb2(sp4, k4)!+sp4*mass
!
!          e2 = g_fbf(edumm, kdumm, sp(:,2), p(:,2), fl2, sp(:,3),p(:,3),fl3,0,0,&
!               &giarray,qiarray(2:3),pol_int) 
!          k2 = p(:,2)+p(:,3)
!          k2sq = sc(k2,k2)
!
!          tmp = vqg(sp4, e2)
!
!
!          if (abs(k4sq) > propcut) then
!             tmp = ci*tmp/k4sq
!          else
!             tmp = czero
!          endif
!
!          if (abs(k2sq) > propcut) then 
!             tmp = -ci/k2sq*tmp
!          else 
!             tmp = czero 
!          endif
!
!          res = res+tmp   !  sw = 1, #1
!
!          ! ---- 1 W below gluon, 1 above, and both below
!
!          if (fl0.eq.'top') flaux = 'bot'
!          if (fl0.eq.'bot') flaux = 'top'
!          sp4 = fW_bffbf(edumm,kdumm,sp,p,fll,flaux,eW(:,2), kW(:,2),&
!               &ng1,ng2,ng3,1,giarray,qiarray,WWid(2),pol_int)    
!
!          k4 = p(:,1)+p(:,2)+p(:,3)+kW(:,2)
!          sp4 = spb2(sp4,k4) !+ mass*sp4
!
!          k4sq = sc(k4,k4)
!
!          tmp = vbqW(sp4,eW(:,1))
!
!          if (abs(k4sq) > propcut) then 
!             tmp = ci/k4sq*tmp
!          else 
!             tmp = czero
!          endif
!
!          res = res + tmp   ! sw=1,#2
!
!
!       endif

       if (sw == 2) then

          ! -- for W above gluon line
          if (fl0 == 'top') flaux = 'bot'
          if (fl0 == 'bot') flaux = 'top'
          if (case_b2 .eqv. .false.) then

             sp4 = fW(edumm,kdumm, sp(:,1),p(:,1), fl1,fl0, eW(:,2),&
                  &kW(:,2),0,giarray,qiarray(1:1),WWid(2), pol_int)
             k4 = p(:,1)+kW(:,2)          
             k4sq = sc(k4,k4)

             sp4 = spb2(sp4, k4)!+mass*sp4

             e2 = gW_fbf(edumm,kdumm, sp(:,2), p(:,2), fl2, sp(:,3), p(:,3), fl3,&
                  &eW(:,1),kW(:,1),0,0, giarray,qiarray(2:3), WWid(1),pol_int)

             k2 = p(:,2)+p(:,3)+kW(:,1)
             k2sq = sc(k2,k2)

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

             res = res + tmp  ! sw =2, #1




             ! ---- for W below gluon line          

             sp4 = fW_bffbf(edumm,kdumm,sp,p,fll,flaux,eW(:,1), kW(:,1),&
                  &ng1,ng2,ng3,3,&
                  &giarray,qiarray,WWid(1),pol_int)    

             k4 = p(:,1)+p(:,2)+p(:,3)+kW(:,1)

             sp4 = spb2(sp4,k4) !+ mass*sp4
             k4sq = sc(k4,k4)

             tmp = vbqW(sp4,eW(:,2))


             if (abs(k4sq) > propcut) then 
                tmp = ci/k4sq*tmp
             else 
                tmp = czero
             endif

             res = res + tmp   ! sw = 2, #2

          elseif (case_b2) then

             sp4 = fW_bffbf(edumm,kdumm,sp,p,fll,flaux,eW(:,1), kW(:,1),&
                  &ng1,ng2,ng3,sw,&
                  &giarray,qiarray,WWid(1),pol_int)    

             k4 = p(:,1)+p(:,2)+p(:,3)+kW(:,1)

             sp4 = spb2(sp4,k4) !+ mass*sp4
             k4sq = sc(k4,k4)

             tmp = vbqW(sp4,eW(:,2))

             if (abs(k4sq) > propcut) then 
                tmp = ci/k4sq*tmp
             else 
                tmp = czero
             endif

             res = res + tmp   ! sw = 2, case_b2, #1

          endif

       endif

!       if (sw ==3) then
!          e2 = gWW_fbf(edumm, kdumm, sp(:,2),p(:,2),fl2,sp(:,3),p(:,3),fl3,& 
!               &eW,kW,0,0,giarray,qiarray(2:3),WWid,pol_int)
!
!          k2 = p(:,2)+p(:,3)+kW(:,1)+kW(:,2)
!          k2sq = sc(k2,k2)
!
!          tmp = vqg(sp(:,1),e2)
!
!
!          if (abs(k2sq) > propcut.and.fl0==fl1) then 
!             tmp = -ci/k2sq*tmp
!          else 
!             tmp = czero 
!             !write(*,*) 'set to zero'
!             !write(*,*) fl0
!          endif
!
!          res = res + tmp  ! #1
!       endif
!
!       ! This case is used for A3.
!       if (sw==4) then 
!
!          eW(:,1) = eWin(:,2)
!          eW(:,2) = eWin(:,1)
!          kW(:,1) = kWin(:,2)
!          kW(:,2) = kWin(:,1)
!!!$          
!          spnew(:,1) = sp(:,3)
!          spnew(:,2) = sp(:,2)
!          spnew(:,3) = sp(:,1)
!          pnew(:,1) = p(:,3)
!          pnew(:,2) = p(:,2)
!          pnew(:,3)= p(:,1)
!          fllnew(1) = fll(3)
!          fllnew(2) = fll(2)
!          fllnew(3) = fll(1)
!
!
!          res = fWW_bffbf(e,k,spnew,pnew,fllnew,fl0,eW,kW,ng1,ng2,ng3,1,giarray,qiarray,WWid, pol_int)
!
!       endif
!
!       ! Also for A3...
!
!       if (sw == 5) then
!!!$
!          eW(:,1) = eWin(:,2)
!          eW(:,2) = eWin(:,1)
!          kW(:,1) = kWin(:,2)
!          kW(:,2) = kWin(:,1)
!
!
!          spnew(:,1) = sp(:,3)
!          spnew(:,2) = sp(:,2)
!          spnew(:,3) = sp(:,1)
!          pnew(:,1) = p(:,3)
!          pnew(:,2) = p(:,2)
!          pnew(:,3)= p(:,1)
!          fllnew(1) = fll(3)
!          fllnew(2) = fll(2)
!          fllnew(3) = fll(1)
!          res = fWW_bffbf(e,k,spnew,pnew,fllnew,fl0,eW,kW,ng1,ng2,ng3,3,giarray,qiarray,WWid, pol_int)
!
!       endif

    else                     ! ngluons > 0


       if (sw == 2) then

          if (fl0 == 'bot') flaux = 'top'
          if (fl0 == 'top') flaux = 'bot'

          do m = 0,ng2

             sp4 = fW(e(:,1:ng1+m),k(:,1:ng1+m),sp(:,1),p(:,1),fl1,fl0,&
                  &eW(:,2),kW(:,2),ng1,giarray(1:ng1+m),qiarray(1:1),WWid(2),&
                  & pol_int)

             k4 = sum(k(:,1:ng1+m),dim=2)+ p(:,1) + kW(:,2)

             k4sq = sc(k4,k4) !- mass2

             sp4 = spb2(sp4,k4) !+ mass*sp4

             e2 = gW_fbf(e(:,ng1+m+1:ngluons),k(:,ng1+m+1:ngluons),sp(:,2),&
                  p(:,2),fll(2),sp(:,3),p(:,3),fll(3), eW(:,1),kW(:,1), &
                  & ng2-m, ng3, giarray(ng1+m+1:ngluons),qiarray(2:3),WWid(1),&
                  & pol_int)

             k2 = sum(k(:,ng1+m+1:ngluons),dim=2)+kW(:,1)+p(:,2)+p(:,3)

             k2sq = sc(k2,k2)

             tmp = vqg(sp4,e2)

             if (abs(k4sq) > propcut) then
                tmp = ci*tmp/k4sq
             else
                tmp = czero
             endif

             if (abs(k2sq) >  propcut) then
                tmp = -ci*tmp/k2sq
             else
                tmp = czero
             endif

             res = res + tmp

          enddo


          do m = 1,ng1

             sp4 = fWW_bffbf(e(:,m+1:ngluons),k(:,m+1:ngluons),sp,p,fll,fl0, &
                  & eW,kW,ng1-m,ng2,ng3,sw,giarray(m+1:ngluons),qiarray,&
                  & WWid,pol_int)

             k4 = sum(k(:,m+1:ngluons),dim=2) + p(:,1)+p(:,2)+p(:,3) + &
                  & kW(:,1) + kW(:,2)

             k4sq = sc(k4,k4) !- mass2

             sp4 = spb2(sp4,k4) !+ mass*sp4

             e2 = vgluon(e(:,1:m),k(:,1:m),giarray(1:m),pol_int)


             k2 = sum(k(:,1:m),dim=2)
             k2sq = sc(k2,k2)

             tmp = vgq(e2,sp4)


             if (abs(k4sq) > propcut) then
                tmp = ci*tmp/k4sq
             else
                tmp = czero
             endif


             if (m>1) then
                if (abs(k2sq) >  propcut) then
                   tmp = -ci*tmp/k2sq
                else
                   tmp = czero
                endif
             endif

             res = res + tmp

          enddo


          do m = 0,ng4-1

             ngL = ng1 + ng2 + ng3 + m

             sp4 = fWW_bffbf(e(:,1:ngL),k(:,1:ngL),sp,p,fll,fl0, &
                  & eW,kW,ng1,ng2,ng3,sw,giarray(1:ngL),qiarray,&
                  & WWid,pol_int)


             k4 = sum(k(:,1:ngL),dim=2) + sum(p(:,1:3),dim=2) + &
                  & sum(kW(:,1:2),dim=2)

             sp4 = spb2(sp4,k4) !+ mass*sp4

             k4sq = sc(k4,k4) !- mass2

             e2 = vgluon(e(:,ngL+1:ngluons),k(:,ngL+1:ngluons),&
                  & giarray(ngL+1:ngluons),pol_int)

             k2 = sum(k(:,ngL+1:ngluons),dim=2)
             k2sq = sc(k2,k2)

             tmp = vqg(sp4,e2)


             if (abs(k4sq) > propcut) then
                tmp = ci*tmp/k4sq
             else
                tmp = czero
             endif


             !if (ngL+1 < ngluons) then
             if (m < ng4-1) then
                if (abs(k2sq) >  propcut) then
                   tmp = -ci*tmp/k2sq
                else
                   tmp = czero
                endif
             endif

             res = res + tmp


          enddo

          sp4 = fW_bffbf(e,k,sp,p,fll,flaux,eW(:,1),kW(:,1),ng1,ng2,ng3,3,&
               & giarray, qiarray,WWid(1),pol_int)

          k4 = sum(k(:,1:ngluons),dim=2) + kW(:,1) + sum(p(:,1:3),dim=2)

          k4sq = sc(k4,k4) !- mass2

          sp4 = spb2(sp4,k4) !+ mass*sp4


          tmp = vbqW(sp4,eW(:,2))

          if (abs(k4sq) > propcut) then
             tmp = ci*tmp/k4sq
          else
             tmp = czero
          endif

          res = res + tmp

       endif
    endif


  end function fWW_bffbf


  end module recurrence