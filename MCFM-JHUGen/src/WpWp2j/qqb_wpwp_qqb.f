  !------------------------------------------------------------------------!
  ! Authors: Tom Melia, Kirill Melnikov, Raoul Rontsch, Giulia Zanderighi  !
  ! Date: 25/10/2010                                                       !
  ! Used for arXiv:1007.5313 Wp Wp 2 jets                                  !
  !------------------------------------------------------------------------!
      subroutine qqb_wpwp_qqb(p,msq)
      use qqqqampl
      use consts_dp
      implicit none
!      include 'constants.f'
      include 'types.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer i,j,k
      double precision msq(-5:5,-5:5),p(12,4)
      double precision mqqb(3),mqbq(3),mqqq(2),mqbb(2)
      double precision mtot(3),mtot_bits(3)
      double precision aveqqqq, fac
      double precision sw1, sw2, facprop
      double complex propw1, propw2
      double complex Amp(2)

      double precision, parameter :: c1=2._dp,c2=-2._dp/3._dp

c---set msq=0 to initialize

      msq=0._dp

      aveqqqq = 1._dp/9._dp/4._dp
      fac = (gw**8)*(gsq**2)/4._dp
      sw1 =two*(p(3,4)*p(4,4)-p(3,1)*p(4,1)-p(3,2)*p(4,2)-p(3,3)*p(4,3))
      sw2 =two*(p(5,4)*p(6,4)-p(5,1)*p(6,1)-p(5,2)*p(6,2)-p(5,3)*p(6,3))
      propw1 = sw1/(sw1-wmass**2+ci*wwidth*wmass)
      propw2 = sw2/(sw2-wmass**2+ci*wwidth*wmass)
      facprop = abs(propw1)**2*abs(propw2)**2
      fac = fac*aveqqqq*facprop


      call getamplqqqq(p(1:8,:),1,2,7,8,Amp(1))
      call getamplqqqq(p(1:8,:),1,8,7,2,Amp(2))
      Amp(2)=-Amp(2)
      mqqb(1) = fac*(c1*abs(Amp(1))**2+c1*abs(Amp(2))**2 +
     . c2*(Amp(1)*conjg(Amp(2))+conjg(Amp(1))*Amp(2)))
      mqqb(2) = fac*(c1*abs(Amp(1))**2)
      mqqb(3) = fac*(c1*abs(Amp(2))**2)

      call getamplqqqq(p(1:8,:),2,1,7,8,Amp(1))
      call getamplqqqq(p(1:8,:),2,8,7,1,Amp(2))
      Amp(2)=-Amp(2)
      mqbq(1) = fac*(c1*abs(Amp(1))**2+c1*abs(Amp(2))**2 +
     . c2*(Amp(1)*conjg(Amp(2))+conjg(Amp(1))*Amp(2)))
      mqbq(2) = fac*(c1*abs(Amp(1))**2)
      mqbq(3) = fac*(c1*abs(Amp(2))**2)

      call getamplqqqq(p(1:8,:),1,7,2,8,Amp(1))
      call getamplqqqq(p(1:8,:),1,8,2,7,Amp(2))
      Amp(2)=-Amp(2)
      mqqq(1) = fac*(c1*abs(Amp(1))**2+c1*abs(Amp(2))**2 +
     . c2*(Amp(1)*conjg(Amp(2))+conjg(Amp(1))*Amp(2)))
      mqqq(2) = fac*(c1*abs(Amp(1))**2)

      call getamplqqqq(p(1:8,:),7,1,8,2,Amp(1))
      call getamplqqqq(p(1:8,:),7,2,8,1,Amp(2))
      Amp(2)=-Amp(2)
      mqbb(1) = fac*(c1*abs(Amp(1))**2+c1*abs(Amp(2))**2 +
     . c2*(Amp(1)*conjg(Amp(2))+conjg(Amp(1))*Amp(2)))
      mqbb(2) = fac*(c1*abs(Amp(1))**2)

      !---fill msq


      msq(2,-1) = mqqb(1) + mqqb(2) ! u dbar initial state
      msq(2,-3) = mqqb(3)           ! u sbar initial state
      msq(4,-3) = mqqb(1) + mqqb(2) ! c sbar initial state
      msq(4,-1) = mqqb(3)           ! c dbar initial state

      msq(-1,2) = mqbq(1) + mqbq(2) ! dbar u initial state
      msq(-3,2) = mqbq(3)           ! sbar u initial state
      msq(-3,4) = mqbq(1) + mqbq(2) ! sbar c intital state
      msq(-1,4) = mqbq(3)           ! dbar c initial state

      msq(2,2) = mqqq(1)*(0.5_dp)  ! u u initial state
      msq(2,4) = mqqq(2)            ! u c initial state
      msq(4,2) = mqqq(2)            ! c u initial state
      msq(4,4) = mqqq(1)*(0.5_dp)  ! c c initial state

      msq(-1,-1) = mqbb(1)*(0.5_dp) ! dbar dbar initial state
      msq(-1,-3) = mqbb(2)           ! dbar sbar initial state
      msq(-3,-1) = mqbb(2)           ! sbar dbar initial state
      msq(-3,-3) = mqbb(1)*(0.5_dp) ! sbar sbar initial state

      return
      end
