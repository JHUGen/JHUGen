  !------------------------------------------------------------------------!
  ! Authors: Tom Melia, Kirill Melnikov, Raoul Rontsch, Giulia Zanderighi  !
  ! Date: 25/10/2010                                                       !
  ! Used for arXiv:1007.5313 Wp Wp 2 jets                                  !
  !------------------------------------------------------------------------!

  module qqqqampl
  use consts_dp
  use recurrence
  use spinfns
  implicit none

  public :: getamplqqqq

  contains

  subroutine getamplqqqq(pin,i1,i2,i3,i4,Amp)
    double precision,intent(in) :: pin(8,4)
    integer, intent(in)       :: i1,i2,i3,i4
    double complex, intent(out)  :: Amp
    !--------------------------------------
    double complex :: p(8,4)
    double complex :: u1(4), v1(4), u2(4), v2(4)
    double complex :: eW1(4), eW2(4), pW1(4), pW2(4)
    double complex :: tmp(4)
    double complex :: edum(4,0),kdum(4,0)
    double complex :: eW(4,2),pW(4,2)
    double complex :: pq(4,3),sp(4,3)
    character(len=3), save :: fl0,fl1(3) 
    integer, save :: giarray(0), qiarray(3), WWid(2), pol_int  
    integer :: iswap(8),i
    logical, save :: firsttime = .true.

    Amp = czero

    if (firsttime) then
       giarray(:) = 0 
       qiarray(:) = (/1,2,4/)
       WWid = (/8,16/) 
       pol_int = 0 
       fl1 = (/'top','bot','top'/)
       fl0 = 'bot'
       firsttime = .false.
       WWqqqq = .true. 
       qbq_WW_and_gluons=.true.
       qbq_and_gluons=.false.
    endif
    ! -- virtual calculation might changes, make sure they are set here 
    case_b1 = .true. 
    case_b2 = .false.

       
    ! complexify momenta and swap spacetime indices

    iswap = (/i1,i2,3,4,5,6,i3,i4/)
    do i=1,8
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


    !call riinitialize_Extmem
    tmp = fWW_bffbf(edum,kdum,sp,pq,fl1,fl0,eW,pW,0,0,0,2,&
         &giarray,qiarray,WWid,pol_int) 


    !-- swap Ws
    eW(:,1) = eW2
    eW(:,2) = eW1 
    pW(:,1) = pW2
    pW(:,2) = pW1 

    !call riinitialize_Extmem
    tmp = tmp+fWW_bffbf(edum,kdum,sp,pq,fl1,fl0,eW,pW,0,0,0,2,&
         &giarray,qiarray,WWid,pol_int) 

    Amp = psp1(tmp,v0(p(1,:),1))  ! psp1(tmp,v1)

  end subroutine getamplqqqq


  end module qqqqampl