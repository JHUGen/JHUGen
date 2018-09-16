      subroutine qqb_ttw_gs(p,msq)
************************************************************************
*     Author: R.K. Ellis                                               *
*     April, 2012.                                                     *
************************************************************************
*     Matrix element SUBTRACTION squared and averaged over initial     *
*     colors and spins                                                 *
*                                                                      *
*     q(-p1)+qbar(-p2) --> t(p345) + t~(p678) + W +g(p11)              *
*                                               |                      *
*                                               --> nu(p9)+e^+(p10)    *
*                                                                      *
*     (positively charged W only)                                      *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'masses.f'
      include 'breit.f'
      integer j,k,nd
c --- remember: nd will count the dipoles
      
      double precision p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      double precision 
     & msq17_2(-nf:nf,-nf:nf),msq27_1(-nf:nf,-nf:nf),
     & msq57_6(-nf:nf,-nf:nf),msq67_5(-nf:nf,-nf:nf),
     & msq17_5(-nf:nf,-nf:nf),msq57_1(-nf:nf,-nf:nf),
     & msq17_6(-nf:nf,-nf:nf),msq67_1(-nf:nf,-nf:nf),
     & msq27_5(-nf:nf,-nf:nf),msq57_2(-nf:nf,-nf:nf),
     & msq27_6(-nf:nf,-nf:nf),msq67_2(-nf:nf,-nf:nf),
     & dummyv(-nf:nf,-nf:nf),
     & sub17_2(4),sub27_1(4),sub57_6(4),sub67_5(4),
     & sub17_5(4),sub57_1(4),sub27_5(4),sub57_2(4),
     & sub17_6(4),sub67_1(4),sub27_6(4),sub67_2(4),dsubv
      double precision oldmass2
      external qqb_ttw,donothing_gvec

c--- Note that, compared with the massless case, the subtractions here
c--- must separate the initial-final and final-initial dipoles
      ndmax=12
c--- note: mass of final state particle for final-initial and
c---  final-final dipoles is passed in (clumsily) via mass2

      qqproc=.true.
      qgproc=.false.
      gqproc=.false.
      ggproc=.false.

c--- save original value of mass2
      oldmass2=mass2

c--- calculate all the initial-initial dipoles
      mass2=0d0
      call      dips_mass( 1,p,1,7,2,sub17_2,dsubv,msq17_2,dummyv,
     . qqb_ttw,donothing_gvec)
      call      dips_mass( 2,p,2,7,1,sub27_1,dsubv,msq27_1,dummyv,
     . qqb_ttw,donothing_gvec)
      
c--- final-final
      mass2=mt   
      call dips_mass( 3,p,5,7,6,sub57_6,dsubv,msq57_6,dummyv,
     . qqb_ttw,donothing_gvec)
      call dips_mass( 4,p,6,7,5,sub67_5,dsubv,msq67_5,dummyv,
     . qqb_ttw,donothing_gvec)

c--- now the basic initial final and final initial
      call dips_mass( 5,p,1,7,5,sub17_5,dsubv,msq17_5,dummyv,
     . qqb_ttw,donothing_gvec)
      call dips_mass( 6,p,5,7,1,sub57_1,dsubv,msq57_1,dummyv,
     . qqb_ttw,donothing_gvec)

      call dips_mass( 7,p,1,7,6,sub17_6,dsubv,msq17_6,dummyv,
     . qqb_ttw,donothing_gvec)
      call dips_mass( 8,p,6,7,1,sub67_1,dsubv,msq67_1,dummyv,
     . qqb_ttw,donothing_gvec)

      call dips_mass( 9,p,2,7,5,sub27_5,dsubv,msq27_5,dummyv,
     . qqb_ttw,donothing_gvec)
      call dips_mass(10,p,5,7,2,sub57_2,dsubv,msq57_2,dummyv,
     . qqb_ttw,donothing_gvec)

      call dips_mass(11,p,2,7,6,sub27_6,dsubv,msq27_6,dummyv,
     . qqb_ttw,donothing_gvec)
      call dips_mass(12,p,6,7,2,sub67_2,dsubv,msq67_2,dummyv,
     . qqb_ttw,donothing_gvec)

c--- reset mass2 to original value
      mass2=oldmass2

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msq(nd,j,k)=0d0
      enddo
      enddo
      enddo

      do j=-nf,nf
      do k=-nf,nf
      
      if     ((j.gt.0) .and. (k.lt.0)) then
c--- quark-antiquark
        msq( 1,j,k)=-sub17_2(qq)*msq17_2(j,k)/xn
        msq( 2,j,k)=-sub27_1(qq)*msq27_1(j,k)/xn
        msq( 3,j,k)=-sub57_6(qq)*msq57_6(j,k)/xn
        msq( 4,j,k)=-sub67_5(qq)*msq67_5(j,k)/xn
        msq( 5,j,k)=+sub17_5(qq)*msq17_5(j,k)*(xn-two/xn)
        msq( 6,j,k)=+sub57_1(qq)*msq57_1(j,k)*(xn-two/xn)
        msq( 7,j,k)=+sub17_6(qq)*msq17_6(j,k)*two/xn
        msq( 8,j,k)=+sub67_1(qq)*msq67_1(j,k)*two/xn
        msq( 9,j,k)=+sub27_5(qq)*msq27_5(j,k)*two/xn
        msq(10,j,k)=+sub57_2(qq)*msq57_2(j,k)*two/xn
        msq(11,j,k)=+sub27_6(qq)*msq27_6(j,k)*(xn-two/xn)
        msq(12,j,k)=+sub67_2(qq)*msq67_2(j,k)*(xn-two/xn)

      elseif((j.lt.0) .and. (k.gt.0)) then
c--- antiquark-quark
        msq( 1,j,k)=-sub17_2(qq)*msq17_2(j,k)/xn
        msq( 2,j,k)=-sub27_1(qq)*msq27_1(j,k)/xn
        msq( 3,j,k)=-sub57_6(qq)*msq57_6(j,k)/xn
        msq( 4,j,k)=-sub67_5(qq)*msq67_5(j,k)/xn
        msq( 5,j,k)=+sub17_5(qq)*msq17_5(j,k)*two/xn
        msq( 6,j,k)=+sub57_1(qq)*msq57_1(j,k)*two/xn
        msq( 7,j,k)=+sub17_6(qq)*msq17_6(j,k)*(xn-two/xn)
        msq( 8,j,k)=+sub67_1(qq)*msq67_1(j,k)*(xn-two/xn)
        msq( 9,j,k)=+sub27_5(qq)*msq27_5(j,k)*(xn-two/xn)
        msq(10,j,k)=+sub57_2(qq)*msq57_2(j,k)*(xn-two/xn)
        msq(11,j,k)=+sub27_6(qq)*msq27_6(j,k)*two/xn
        msq(12,j,k)=+sub67_2(qq)*msq67_2(j,k)*two/xn

      elseif (j.eq.0) then
        if (k.gt.0) then
c--- gluon-quark
          msq(1,j,k)=sub17_2(qg)
     &    *(msq17_2(-1,k)+msq17_2(-2,k)+msq17_2(-3,k)+msq17_2(-4,k)
     .     +msq17_2(-5,k))

        elseif (k.lt.0) then
c--- gluon-antiquark
          msq(1,j,k)=sub17_2(qg)
     &    *(msq17_2(+1,k)+msq17_2(+2,k)+msq17_2(+3,k)+msq17_2(+4,k)
     .    +msq17_2(+5,k))
        endif 

      elseif (k.eq.0) then
        if (j.gt.0) then
c--- quark-gluon
          msq(2,j,k)=sub27_1(qg)
     &    *(msq27_1(j,-1)+msq27_1(j,-2)+msq27_1(j,-3)+msq27_1(j,-4)
     .    +msq27_1(j,-5))

        elseif (j.lt.0) then
c--- antiquark-gluon
          msq(2,j,k)=sub27_1(qg)
     &    *(msq27_1(j,+1)+msq27_1(j,+2)+msq27_1(j,+3)+msq27_1(j,+4)
     .     +msq27_1(j,+5))
        endif 
      endif

      enddo
      enddo

      return
      end





