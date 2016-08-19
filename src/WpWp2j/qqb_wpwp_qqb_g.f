  !------------------------------------------------------------------------!
  ! Authors: Tom Melia, Kirill Melnikov, Raoul Rontsch, Giulia Zanderighi  !
  ! Date: 25/10/2010                                                       !
  ! Used for arXiv:1007.5313 Wp Wp 2 jets                                  !
  !------------------------------------------------------------------------!
      subroutine qqb_wpwp_qqb_g(p,msq,chn)
      use qqqqgampl
      use consts_dp
      implicit none
c      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      character(len=3) :: chn 
      integer i,j,k
      double precision msq(-5:5,-5:5),p(12,4)
      double precision mqqb(3),mqbq(3),mqqq(3),mqbb(3)
      double precision mqgl(3),mglq(3),mqbg(3),mgqb(3)
      double precision mtot(3)
      double precision aveqq, aveqg,fac
      double precision sw1, sw2, facprop
      double complex propw1, propw2

c---set msq=0 to initialize

      msq=0d0

      aveqq = 1d0/9d0/4d0
      aveqg = 1d0/3d0/8d0/4d0
      fac = (gw**8)*(gsq**2)/4d0 ! Born factor 
      fac = fac*gsq*2d0 ! 2 dues to TA normalization 
      fac = fac/4d0  

      sw1 =2d0*(p(3,4)*p(4,4)-p(3,1)*p(4,1)-p(3,2)*p(4,2)-p(3,3)*p(4,3))
      sw2 =2d0*(p(5,4)*p(6,4)-p(5,1)*p(6,1)-p(5,2)*p(6,2)-p(5,3)*p(6,3))
      propw1 = sw1/(sw1-wmass**2+ci*wwidth*wmass)
      propw2 = sw2/(sw2-wmass**2+ci*wwidth*wmass)
      facprop = abs(propw1)**2*abs(propw2)**2
      fac = fac*facprop


      call mtotqqqqg(p(1:9,:),mtot,'qqb')

         mqqb(1) = fac*aveqq*mtot(1)
         mqqb(2) = fac*aveqq*mtot(2)
         mqqb(3) = fac*aveqq*mtot(3)

      call mtotqqqqg(p(1:9,:),mtot,'qbq')

         mqbq(1) = fac*aveqq*mtot(1)
         mqbq(2) = fac*aveqq*mtot(2)
         mqbq(3) = fac*aveqq*mtot(3)

      call mtotqqqqg(p(1:9,:),mtot,'qqq')

         mqqq(1) = fac*aveqq*mtot(1)
         mqqq(2) = fac*aveqq*mtot(2)
         mqqq(3) = fac*aveqq*mtot(3)

      call mtotqqqqg(p(1:9,:),mtot,'qbb')

         mqbb(1) = fac*aveqq*mtot(1)
         mqbb(2) = fac*aveqq*mtot(2)
         mqbb(3) = fac*aveqq*mtot(3)

      call mtotqqqqg(p(1:9,:),mtot,'qgl')

         mqgl(1) = fac*aveqg*mtot(1)
         mqgl(2) = fac*aveqg*mtot(2)
         mqgl(3) = fac*aveqg*mtot(3)
         
      call mtotqqqqg(p(1:9,:),mtot,'glq')

         mglq(1) = fac*aveqg*mtot(1)
         mglq(2) = fac*aveqg*mtot(2)
         mglq(3) = fac*aveqg*mtot(3)

      call mtotqqqqg(p(1:9,:),mtot,'qbg')

         mqbg(1) = fac*aveqg*mtot(1)
         mqbg(2) = fac*aveqg*mtot(2)
         mqbg(3) = fac*aveqg*mtot(3)

      call mtotqqqqg(p(1:9,:),mtot,'gqb')

         mgqb(1) = fac*aveqg*mtot(1)
         mgqb(2) = fac*aveqg*mtot(2)
         mgqb(3) = fac*aveqg*mtot(3)

c      Now fill in the matrix elements

      
         msq(2,-1) = mqqb(1) + mqqb(2)               ! u dbar initial state
         msq(2,-3) = mqqb(3)                         ! u sbar initial state
         msq(4,-3) = mqqb(1) + mqqb(2)               ! c sbar initial state
         msq(4,-1) = mqqb(3)                         ! c dbar initial state

         msq(-1,2) = mqbq(1) + mqbq(2)               ! dbar u initial state
         msq(-3,2) = mqbq(3)                         ! sbar u initial state
         msq(-3,4) = mqbq(1) + mqbq(2)               ! sbar c intital state
         msq(-1,4) = mqbq(3)                         ! dbar c initial state
      
         msq(2,2) = mqqq(1)*(1d0/2d0)                ! u u initial state
         msq(2,4) = mqqq(2)                          ! u c initial state
         msq(4,2) = mqqq(2)                          ! c u initial state
         msq(4,4) = mqqq(1)*(1d0/2d0 )               ! c c initial state

         msq(-1,-1) = mqbb(1)*(1d0/2d0)              ! dbar dbar initial state
         msq(-1,-3) = mqbb(2)                        ! dbar sbar initial state
         msq(-3,-1) = mqbb(2)                        ! sbar dbar initial state
         msq(-3,-3) = mqbb(1)*(1d0/2d0)              ! sbar sbar initial state

         msq(1,0) = 0d0                              ! d g initial state
         msq(2,0) = mqgl(1)/2d0 + mqgl(2)            ! u g initial state
         msq(3,0) = 0d0                              ! s g initial state
         msq(4,0) = mqgl(1)/2d0 + mqgl(2)            ! c g initial state

         msq(0,1) = 0d0                             ! g d initial state
         msq(0,2) = mglq(1)/2d0 + mglq(2)           ! g u initial state
         msq(0,3) = msq(0,1)                        ! g s initial state
         msq(0,4) = msq(0,2)                        ! g c initial state

         msq(-1,0) = mqbg(1)/2d0 + mqbg(2)          ! db g initial state
         msq(-2,0) = 0d0                            ! ub g initial state
         msq(-3,0) = msq(-1,0)                      ! sb g initial state
         msq(-4,0) = msq(-2,0)                      ! cb g initial state

         msq(0,-1) = mgqb(1)/2d0 + mgqb(2)           ! g db initial state
         msq(0,-2) = 0d0                             ! g ub initial state
         msq(0,-3) = msq(0,-1)                       ! g sb initial state
         msq(0,-4) = msq(0,-2)                       ! g cb initial state


      return
      end
