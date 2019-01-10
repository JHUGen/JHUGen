  !------------------------------------------------------------------------!
  ! Authors: Tom Melia, Kirill Melnikov, Raoul Rontsch, Giulia Zanderighi  !
  ! Date: 25/10/2010                                                       !
  ! Used for arXiv:1007.5313 Wp Wp 2 jets                                  !
  !------------------------------------------------------------------------!

  module qqqqgampl
  use consts_dp
  use recurrence
  use spinfns
  implicit none
  private
  public :: getamplqqqqg, getamplqqqqg_stu, mtotqqqqg

contains


  !------------------------------!        
  !getamplqqqqg returns          !
  !------------------------------!
  !     q(2)|            |qb(7)  !
  !         |            |       !
  !         ^            v       !
  !    g(9) |~~~~~~~~~~~~|       !
  !         |            |       !
  !         ^            v       !
  !         |            |       !
  !    qb(1)|            |q(8)   !
  !             Amp(1)           !
  !------------------------------!
  !------------------------------!
  !     q(2)|            |qb(7)  !
  !         |   g(9)     |       !
  !         ^            v       !
  !         |~~~~~~~~~~~~|       !
  !         |            |       !
  !         ^            v       !
  !         |            |       !
  !    qb(1)|            |q(8)   !
  !             Amp(2)           !
  !------------------------------!
  !------------------------------!
  !     q(2)|            |qb(7)  !
  !         |            |       !
  !         ^            v       !
  !         |~~~~~~~~~~~~|  g(9) !
  !         |            |       !
  !         ^            v       !
  !         |            |       !
  !    qb(1)|            |q(8)   !
  !             Amp(3)           !
  !------------------------------!
  !------------------------------!
  !     q(2)|            |qb(7)  !
  !         |            |       !
  !         ^            v       !
  !         |~~~~~~~~~~~~|       !
  !         |            |       !
  !         ^    g(9)    v       !
  !         |            |       !
  !    qb(1)|            |q(8)   !
  !             Amp(4)           !
  !------------------------------!



  subroutine getamplqqqqg(pin,i1,i2,i3,i4,i5,hlg,Amp)
    double precision,intent(in)       :: pin(9,4)
    integer, intent(in)       :: i1,i2,i3,i4,i5
    integer, intent(in)       :: hlg
    double complex, intent(out)  :: Amp(4)
    !--------------------------------------
    double complex :: p(9,4)
    double complex :: u1(4), v1(4), u2(4), v2(4)
    double complex :: eW1(4), eW2(4), pW1(4), pW2(4)
    double complex :: tmp(4)
    double complex :: eg(4,1),pg(4,1)
    double complex :: eW(4,2),pW(4,2)
    double complex :: pq(4,3),sp(4,3)
    character(len=3),save :: fl0,fl1(3) 
    integer,save :: giarray(1), qiarray(3), WWid(2), pol_int  
    integer :: iswap(9),i
    logical, save :: firsttime = .true.

    Amp = czero

    if (firsttime) then
       giarray(:) = 32 
       qiarray(:) = (/1,2,4/)
       WWid = (/8,16/) 
       pol_int = 0 
       case_b1 = .true. 
       case_b2 = .false.
       case_a3=.false.
       WWqqqq = .true. 
       qbq_WW_and_gluons=.true.
       qbq_and_gluons=.false.
       fl1 = (/'top','bot','top'/)
       fl0 = 'bot'
       firsttime = .false.
    endif


    ! complexify momenta and swap spacetime indices

    iswap = (/i1,i2,3,4,5,6,i3,i4,i5/)
    do i=1,9
       p(i,1) = dcmplx(pin(iswap(i),4),0d0)
       p(i,2) = dcmplx(pin(iswap(i),1),0d0)
       p(i,3) = dcmplx(pin(iswap(i),2),0d0)
       p(i,4) = dcmplx(pin(iswap(i),3),0d0)
    enddo


    ! set Ws polarisation
    eW1 = pol_dk2mom(p(3,:),p(4,:),-1,.true.)
    eW2 = pol_dk2mom(p(5,:),p(6,:),-1,.true.)
    pW1 = p(3,:)+p(4,:)
    pW2 = p(5,:)+p(6,:)
    
    eW(:,1) = eW1
    eW(:,2) = eW2 
    pW(:,1) = pW1
    pW(:,2) = pW2 

    
    pq(:,1) = p(2,:) 
    pq(:,2) = p(7,:) 
    pq(:,3) = p(8,:) 
    sp(:,1) = ubar0(p(2,:),-1) ! u1
    sp(:,2) = v0(p(7,:),1)     ! v2
    sp(:,3) = ubar0(p(8,:),-1) ! u2
    pg(:,1) = p(9,:)
    eg(:,1) = pol_mless(p(9,:),hlg,.true.)

    !call riinitialize_Extmem
    tmp = fWW_bffbf(eg,pg,sp,pq,fl1,fl0,eW,pW,1,0,0,2,&
         &giarray,qiarray,WWid,pol_int) 
    Amp(1) = psp1(tmp,v0(p(1,:),1))


    !call riinitialize_Extmem
    tmp = fWW_bffbf(eg,pg,sp,pq,fl1,fl0,eW,pW,0,1,0,2,&
         &giarray,qiarray,WWid,pol_int) 

    Amp(2) = psp1(tmp,v0(p(1,:),1))

    !call riinitialize_Extmem
    tmp = fWW_bffbf(eg,pg,sp,pq,fl1,fl0,eW,pW,0,0,1,2,&
         &giarray,qiarray,WWid,pol_int) 

    Amp(3) = psp1(tmp,v0(p(1,:),1))

    !call riinitialize_Extmem
    tmp = fWW_bffbf(eg,pg,sp,pq,fl1,fl0,eW,pW,0,0,0,2,&
         &giarray,qiarray,WWid,pol_int) 

    Amp(4) = psp1(tmp,v0(p(1,:),1))

    
    !-- swap Ws
    eW(:,1) = eW2
    eW(:,2) = eW1 
    pW(:,1) = pW2
    pW(:,2) = pW1 

    !call riinitialize_Extmem
    tmp = fWW_bffbf(eg,pg,sp,pq,fl1,fl0,eW,pW,1,0,0,2,&
         &giarray,qiarray,WWid,pol_int) 

    Amp(1) = Amp(1) + psp1(tmp,v0(p(1,:),1))

    !call riinitialize_Extmem
    tmp = fWW_bffbf(eg,pg,sp,pq,fl1,fl0,eW,pW,0,1,0,2,&
         &giarray,qiarray,WWid,pol_int) 

    Amp(2) = Amp(2) + psp1(tmp,v0(p(1,:),1))

    !call riinitialize_Extmem
    tmp = fWW_bffbf(eg,pg,sp,pq,fl1,fl0,eW,pW,0,0,1,2,&
         &giarray,qiarray,WWid,pol_int) 

    Amp(3) = Amp(3) + psp1(tmp,v0(p(1,:),1))

    !call riinitialize_Extmem
    tmp = fWW_bffbf(eg,pg,sp,pq,fl1,fl0,eW,pW,0,0,0,2,&
         &giarray,qiarray,WWid,pol_int) 

    Amp(4) = Amp(4) + psp1(tmp,v0(p(1,:),1))

  end subroutine getamplqqqqg



  subroutine getamplqqqqg_stu(pin,chn,hlg,Amp)
    double precision,intent(in)       :: pin(9,4)
    character, intent(in)     :: chn*3
    integer, intent(in)       :: hlg
    double complex, intent(out)  :: Amp(4,2)

    Amp = czero

    ! --- now calculate the Amps
    if (chn .eq. 'qqb') then

       call getamplqqqqg(pin,1,2,7,8,9,hlg,Amp(:,1))
       call getamplqqqqg(pin,1,8,7,2,9,hlg,Amp(:,2))

    elseif (chn .eq. 'qbq') then

       call getamplqqqqg(pin,2,1,7,8,9,hlg,Amp(:,1))
       call getamplqqqqg(pin,2,8,7,1,9,hlg,Amp(:,2))

    elseif (chn .eq. 'qqq') then

       call getamplqqqqg(pin,1,7,2,8,9,hlg,Amp(:,1))
       call getamplqqqqg(pin,1,8,2,7,9,hlg,Amp(:,2))

    elseif (chn .eq. 'qbb') then

       call getamplqqqqg(pin,7,1,8,2,9,hlg,Amp(:,1))
       call getamplqqqqg(pin,7,2,8,1,9,hlg,Amp(:,2))

    elseif (chn .eq. 'qgl') then

       call getamplqqqqg(pin,1,9,7,8,2,hlg,Amp(:,1))
       call getamplqqqqg(pin,1,8,7,9,2,hlg,Amp(:,2))
       
    elseif (chn .eq. 'glq') then

       call getamplqqqqg(pin,2,9,7,8,1,hlg,Amp(:,1))
       call getamplqqqqg(pin,2,8,7,9,1,hlg,Amp(:,2))

    elseif (chn .eq. 'qbg') then

       call getamplqqqqg(pin,9,1,7,8,2,hlg,Amp(:,1))
       call getamplqqqqg(pin,9,8,7,1,2,hlg,Amp(:,2))

    elseif (chn .eq. 'gqb') then

       call getamplqqqqg(pin,9,2,7,8,1,hlg,Amp(:,1))
       call getamplqqqqg(pin,9,8,7,2,1,hlg,Amp(:,2))

    endif
   Amp(:,2) = -Amp(:,2)  ! (-1) as swap fermions
       
  end subroutine getamplqqqqg_stu


  ! Schematically: mtot(1) = |s+u|^2 , mtot(1) = |s|^2 , mtot(1) = |u|^2 
  subroutine mtotqqqqg(p,mtot,chn)
    double precision, intent(in)  :: p(9,4)
    character, intent(in) :: chn*3
    double precision, intent(out) :: mtot(3)
    !--------------------------------
!    double complex :: Amp(4,2),tmp(3) 
    double complex :: Amp(4,2)
    double precision    :: tmp(3) 
    double precision, parameter    :: nc = 3d0 , Cf = 4d0/3d0
    integer :: i, hlg  

    mtot = 0d0
    Amp = czero
   
    do hlg = -1,1,2

       call getamplqqqqg_stu(p,chn,hlg,Amp)

       do i=1,2 
          tmp(i+1) =  &
               &   nc**2 *(abs(Amp(4,i))**2+abs(Amp(2,i))**2)+&
              &        (abs(Amp(1,i))**2+abs(Amp(3,i))**2)+&
              & 2d0*real((Amp(3,i)+Amp(1,i))*dconjg(Amp(4,i)+Amp(2,i)))

       enddo
       
       tmp(1) =  +2d0*nc*(real(Amp(4,1)*dconjg(Amp(4,2)+Amp(1,2)+Amp(2,2)))&
           & + real(Amp(2,1)*dconjg(Amp(4,2)+Amp(2,2)+Amp(3,2)))&
            & + real(Amp(1,1)*dconjg(Amp(4,2))) &
            & + real(Amp(3,1)*dconjg(Amp(2,2)))) & 
            + 2d0/nc*dreal((Amp(3,1)+Amp(1,1))*dconjg(Amp(3,2)+Amp(1,2)))

       tmp(1) = tmp(1)+tmp(2)+tmp(3)
       

       mtot = mtot + tmp 
          
    enddo

    mtot = mtot*Cf
    

  end subroutine mtotqqqqg

end module qqqqgampl
