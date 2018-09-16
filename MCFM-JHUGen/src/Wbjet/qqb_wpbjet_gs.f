************************************************************************
*     This is the W+ routine                                           *
************************************************************************
      subroutine qqb_wpbjet_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell                                            *
*     January, 2004.                                                   *
************************************************************************
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+b(-p2) -->  W + f(p5) + b(p6) + g(p7)
c                        |
c                        --> nu(p3) + e^+(p4)
c
c--- all momenta are incoming
c
c--- Extended to include charm quark production via the variable "flav"


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'heavyflav.f'

      integer:: n16_2,n17_2,n15_2,n26_1,
     & n27_1,n25_1,n57_6,n67_5,
     & n17_5,n57_1,n17_6,n67_1,n27_6,n67_2,n27_5,n57_2
      parameter(
     & n16_2=1,n17_2=2,n15_2=3,n26_1=4,
     & n27_1=5,n25_1=6,n57_6=7,n67_5=8,
     & n17_5=9,n57_1=9,
     & n17_6=10,n67_1=10,
     & n27_6=11,n67_2=11,
     & n27_5=12,n57_2=12)

      integer:: j,k
c --- remember: nd will count the dipoles
      integer:: nd,isub
c--- slightly obtuse notation, to simplify declaration lines
      real(dp):: p(mxpart,4),msq(maxd,fn:nf,fn:nf)
      real(dp)::
     & msq17_2(fn:nf,fn:nf),msq27_1(fn:nf,fn:nf),
     & msq15_2(fn:nf,fn:nf),msq25_1(fn:nf,fn:nf),
     & msq16_2(fn:nf,fn:nf),msq26_1(fn:nf,fn:nf),
     & msq17_5(fn:nf,fn:nf),msq17_6(fn:nf,fn:nf),
     & msq27_5(fn:nf,fn:nf),msq27_6(fn:nf,fn:nf),
     & msq57_6(fn:nf,fn:nf),msq67_5(fn:nf,fn:nf),
     & dummy(fn:nf,fn:nf),dummyv(fn:nf,fn:nf),
     & sub17_2(4),sub27_1(4),sub15_2(4),sub25_1(4),
     & sub16_2(4),sub26_1(4),
     & sub17_5(4),sub57_1(4),sub27_5(4),sub57_2(4),
     & sub17_6(4),sub67_1(4),sub27_6(4),sub67_2(4),
     & sub57_6(4),sub67_5(4),dsubv

      real(dp)::
     & m17_2(0:2,fn:nf,fn:nf),m27_1(0:2,fn:nf,fn:nf),
     & m16_2(0:2,fn:nf,fn:nf),m26_1(0:2,fn:nf,fn:nf),
     & m15_2(0:2,fn:nf,fn:nf),m25_1(0:2,fn:nf,fn:nf),
     & m57_6(0:2,fn:nf,fn:nf),m67_5(0:2,fn:nf,fn:nf),
     & m17_5(0:2,fn:nf,fn:nf),
     & m17_6(0:2,fn:nf,fn:nf),
     & m27_5(0:2,fn:nf,fn:nf),
     & m27_6(0:2,fn:nf,fn:nf),
     & mqq(0:2,fn:nf,fn:nf)

      common/isub/isub
      external qqb_wbjet,donothing_gvec

      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer,parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      if (isub == 2) then
        ndmax=6
      else
        ndmax=12
      endif

c-- initialize the matrix elements to zero
      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msq(nd,j,k)=0._dp
      enddo
      enddo
      enddo

c--- calculate all the dipoles
c--- the dipole number relates the matrix elements to the transformed
c--- momenta, used in realint to perform cuts and clustering

c--- NEW ATTEMPT - mqq plays the role of a dummy variable here
      call dips(n16_2,p,1,6,2,sub16_2,dsubv,msq16_2,dummyv,
     & qqb_wbjet,donothing_gvec)
      call storedip_mass(m16_2,mqq)
      call dips(n17_2,p,1,7,2,sub17_2,dsubv,msq17_2,dummyv,
     & qqb_wbjet,donothing_gvec)
      call storedip_mass(m17_2,mqq)
      call dips(n15_2,p,1,5,2,sub15_2,dsubv,msq15_2,dummyv,
     & qqb_wbjet,donothing_gvec)
      call storedip_mass(m15_2,mqq)
      call dips(n26_1,p,2,6,1,sub26_1,dsubv,msq26_1,dummyv,
     & qqb_wbjet,donothing_gvec)
      call storedip_mass(m26_1,mqq)
      call dips(n27_1,p,2,7,1,sub27_1,dsubv,msq27_1,dummyv,
     & qqb_wbjet,donothing_gvec)
      call storedip_mass(m27_1,mqq)
      call dips(n25_1,p,2,5,1,sub25_1,dsubv,msq25_1,dummyv,
     & qqb_wbjet,donothing_gvec)
      call storedip_mass(m25_1,mqq)

      if (isub == 1) then
      call dips(n57_6,p,5,7,6,sub57_6,dsubv,msq57_6,dummyv,
     & qqb_wbjet,donothing_gvec)
      call storedip_mass(m57_6,mqq)
      call dips(n67_5,p,6,7,5,sub67_5,dsubv,msq67_5,dummyv,
     & qqb_wbjet,donothing_gvec)
      call storedip_mass(m67_5,mqq)

      call dips(n17_5,p,1,7,5,sub17_5,dsubv,msq17_5,dummyv,
     & qqb_wbjet,donothing_gvec)
      call storedip_mass(m17_5,mqq)
      call dips(n57_1,p,5,7,1,sub57_1,dsubv,dummy,dummyv,
     & qqb_wbjet,donothing_gvec)

      call dips(n17_6,p,1,7,6,sub17_6,dsubv,msq17_6,dummyv,
     & qqb_wbjet,donothing_gvec)
      call storedip_mass(m17_6,mqq)
      call dips(n67_1,p,6,7,1,sub67_1,dsubv,dummy,dummyv,
     & qqb_wbjet,donothing_gvec)

      call dips(n27_6,p,2,7,6,sub27_6,dsubv,msq27_6,dummyv,
     & qqb_wbjet,donothing_gvec)
      call storedip_mass(m27_6,mqq)
      call dips(n67_2,p,6,7,2,sub67_2,dsubv,dummy,dummyv,
     & qqb_wbjet,donothing_gvec)

      call dips(n27_5,p,2,7,5,sub27_5,dsubv,msq27_5,dummyv,
     & qqb_wbjet,donothing_gvec)
      call storedip_mass(m27_5,mqq)
      call dips(n57_2,p,5,7,2,sub57_2,dsubv,dummy,dummyv,
     & qqb_wbjet,donothing_gvec)
      endif

C It may be helpful to have the ordering of the quarks explicit
C Lowest order case W+
c      msq(-3,-5)=Vsum(-3)*sdbbb_veepbbub(p1,p2,p3,p4,p5,p6)*gsq**2
c      msq(-1,-5)=Vsum(-1)*sdbbb_veepbbub(p1,p2,p3,p4,p5,p6)*gsq**2
c      msq(-5,-3)=Vsum(-3)*sdbbb_veepbbub(p2,p1,p3,p4,p5,p6)*gsq**2
c      msq(-5,-1)=Vsum(-1)*sdbbb_veepbbub(p2,p1,p3,p4,p5,p6)*gsq**2

c      msq(-3,+5)=Vsum(-3)*sdbb_veepbub(p1,p2,p3,p4,p5,p6)*gsq**2
c      msq(-1,+5)=Vsum(-1)*sdbb_veepbub(p1,p2,p3,p4,p5,p6)*gsq**2
c      msq(+5,-3)=Vsum(-3)*sdbb_veepbub(p2,p1,p3,p4,p5,p6)*gsq**2
c      msq(+5,-1)=Vsum(-1)*sdbb_veepbub(p2,p1,p3,p4,p5,p6)*gsq**2

c      msq(+2,-5)=Vsum(2)*subb_veepbbd(p1,p2,p3,p4,p5,p6)*gsq**2
c      msq(+4,-5)=Vsum(4)*subb_veepbbd(p1,p2,p3,p4,p5,p6)*gsq**2
c      msq(-5,+2)=Vsum(2)*subb_veepbbd(p2,p1,p3,p4,p5,p6)*gsq**2
c      msq(-5,+4)=Vsum(4)*subb_veepbbd(p2,p1,p3,p4,p5,p6)*gsq**2

c      msq(2,5)=Vsum(2)*sub_veepbd(p1,p2,p3,p4,p5,p6)*gsq**2
c      msq(4,5)=Vsum(4)*sub_veepbd(p1,p2,p3,p4,p5,p6)*gsq**2
c      msq(5,2)=Vsum(2)*sub_veepbd(p2,p1,p3,p4,p5,p6)*gsq**2
c      msq(5,4)=Vsum(2)*sub_veepbd(p2,p1,p3,p4,p5,p6)*gsq**2


C NLO order case W+
c      msq(-3,-5)=Vsum(-3)*sdbbb_veepbbubg(p1,p2,p3,p4,p5,p6,p7)*gsq**3
c      msq(-1,-5)=Vsum(-1)*sdbbb_veepbbubg(p1,p2,p3,p4,p5,p6,p7)*gsq**3
c      msq(-5,-3)=Vsum(-3)*sdbbb_veepbbubg(p2,p1,p3,p4,p5,p6,p7)*gsq**3
c      msq(-5,-1)=Vsum(-1)*sdbbb_veepbbubg(p2,p1,p3,p4,p5,p6,p7)*gsq**3

c      msq(-3,+5)=Vsum(-3)*sdbb_veepbubg(p1,p2,p3,p4,p5,p6,p7)*gsq**3
c      msq(-1,+5)=Vsum(-1)*sdbb_veepbubg(p1,p2,p3,p4,p5,p6,p7)*gsq**3
c      msq(+5,-3)=Vsum(-3)*sdbb_veepbubg(p2,p1,p3,p4,p5,p6,p7)*gsq**3
c      msq(+5,-1)=Vsum(-1)*sdbb_veepbubg(p2,p1,p3,p4,p5,p6,p7)*gsq**3

c      msq(+2,-5)=Vsum(2)*subb_veepbbdg(p1,p2,p3,p4,p5,p6,p7)*gsq**3
c      msq(+4,-5)=Vsum(4)*subb_veepbbdg(p1,p2,p3,p4,p5,p6,p7)*gsq**3
c      msq(-5,+2)=Vsum(2)*subb_veepbbdg(p2,p1,p3,p4,p5,p6,p7)*gsq**3
c      msq(-5,+4)=Vsum(4)*subb_veepbbdg(p2,p1,p3,p4,p5,p6,p7)*gsq**3

c      msq(+2,+5)=Vsum(2)*sub_veepbdg(p1,p2,p3,p4,p5,p6,p7)*gsq**3
c      msq(+4,+5)=Vsum(4)*sub_veepbdg(p1,p2,p3,p4,p5,p6,p7)*gsq**3
c      msq(+5,+2)=Vsum(2)*sub_veepbdg(p2,p1,p3,p4,p5,p6,p7)*gsq**3
c      msq(+5,+4)=Vsum(4)*sub_veepbdg(p2,p1,p3,p4,p5,p6,p7)*gsq**3

c      msq(0,+5)=(Vsum(-1)+Vsum(-3))
c     & *sgb_veepbdub(p1,p2,p3,p4,p5,p6,p7)*gsq**3
c      msq(+5,0)=(Vsum(-1)+Vsum(-3))
c     & *sgb_veepbdub(p2,p1,p3,p4,p5,p6,p7)*gsq**3
c      msq(0,-5)=(Vsum(-1)+Vsum(-3))
c     & *sgbb_veepbbubd(p1,p2,p3,p4,p5,p6,p7)*gsq**3
c      msq(-5,0)=(Vsum(-1)+Vsum(-3))
c     & *sgbb_veepbbubd(p2,p1,p3,p4,p5,p6,p7)*gsq**3
c      elseif (isub == 2) then
c      msq(0,+2)=Vsum(+2)*sgu_veepbbbd(p1,p2,p3,p4,p5,p6,p7)*gsq**3
c      msq(0,+4)=Vsum(+4)*sgu_veepbbbd(p1,p2,p3,p4,p5,p6,p7)*gsq**3
c      msq(+2,0)=Vsum(+2)*sgu_veepbbbd(p2,p1,p3,p4,p6,p5,p7)*gsq**3
c      msq(+4,0)=Vsum(+4)*sgu_veepbbbd(p2,p1,p3,p4,p6,p5,p7)*gsq**3
c      msq(0,-1)=Vsum(-1)*sgdb_veepbbbub(p1,p2,p3,p4,p6,p5,p7)*gsq**3
c      msq(0,-3)=Vsum(-3)*sgdb_veepbbbub(p1,p2,p3,p4,p6,p5,p7)*gsq**3
c      msq(-1,0)=Vsum(-1)*sgdb_veepbbbub(p2,p1,p3,p4,p5,p6,p7)*gsq**3
c      msq(-3,0)=Vsum(-3)*sgdb_veepbbbub(p2,p1,p3,p4,p5,p6,p7)*gsq**3

c--- fill the dipole contributions
      do j=-nf,nf
      do k=-nf,nf

      if ((abs(j) .ne. flav) .and. (abs(k) .ne. flav)) goto 97
      if ((abs(j) == flav) .and. (abs(k) == flav)) goto 97
c--- so that either abs(j) or abs(k) = flav (but not both).

c--- contributions with 1 b-quark
      if (isub == 1) then

c--- Note that each term must be examined individually to decide
c--- whether the leading colour term is proportional to (1) or (2)
      if ((j > 0) .and. (k>0)) then
      if (j == +flav) then
c---QQ - process bu -> bd
      msq(n57_6,j,k)=msq(n57_6,j,k)
     & +sub57_6(qq)
     & *((xn+1._dp/xn)*m57_6(0,j,k)+two*(m57_6(1,j,k)+m57_6(2,j,k))/xn)
      msq(n67_5,j,k)=msq(n67_5,j,k)
     & +sub67_5(qq)
     & *((xn+1._dp/xn)*m67_5(0,j,k)+two*(m67_5(1,j,k)+m67_5(2,j,k))/xn)
      msq(n17_5,j,k)=msq(n17_5,j,k)
     & +sub17_5(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
      msq(n17_6,j,k)=msq(n17_6,j,k)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
      msq(n27_6,j,k)=msq(n27_6,j,k)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
      msq(n27_5,j,k)=msq(n27_5,j,k)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
      msq(n17_2,j,k)=msq(n17_2,j,k)
     & +sub17_2(qq)
     & *((xn+1._dp/xn)*m17_2(0,j,k)+two*(m17_2(1,j,k)+m17_2(2,j,k))/xn)
      msq(n27_1,j,k)=msq(n27_1,j,k)
     & +sub27_1(qq)
     & *((xn+1._dp/xn)*m27_1(0,j,k)+two*(m27_1(1,j,k)+m27_1(2,j,k))/xn)
      else
c---QQ - process ub -> bd (obtained from above by switching p1 and p2)
c--- corrected this one so far (for 1 and 2 colour structure switch)
      msq(n57_6,j,k)=msq(n57_6,j,k)
     & +sub57_6(qq)
     & *((xn+1._dp/xn)*m57_6(0,j,k)+two*(m57_6(2,j,k)+m57_6(1,j,k))/xn)
      msq(n67_5,j,k)=msq(n67_5,j,k)
     & +sub67_5(qq)
     & *((xn+1._dp/xn)*m67_5(0,j,k)+two*(m67_5(2,j,k)+m67_5(1,j,k))/xn)
      msq(n27_5,j,k)=msq(n27_5,j,k)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
      msq(n27_6,j,k)=msq(n27_6,j,k)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
      msq(n17_6,j,k)=msq(n17_6,j,k)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
      msq(n17_5,j,k)=msq(n17_5,j,k)
     & +sub17_5(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
      msq(n17_2,j,k)=msq(n17_2,j,k)
     & +sub17_2(qq)
     & *((xn+1._dp/xn)*m17_2(0,j,k)+two*(m17_2(2,j,k)+m17_2(1,j,k))/xn)
      msq(n27_1,j,k)=msq(n27_1,j,k)
     & +sub27_1(qq)
     & *((xn+1._dp/xn)*m27_1(0,j,k)+two*(m27_1(2,j,k)+m27_1(1,j,k))/xn)
      endif
      elseif ((j < 0) .and. (k<0)) then
      if (j == -flav) then
c--- QbarQbar - process bd -> bu (obtained from bu -> bd by switching
c---                               p1 and p5, p2 and p6)
      msq(n57_6,j,k)=msq(n57_6,j,k)
     & +sub57_6(qq)
     & *((xn+1._dp/xn)*m57_6(0,j,k)+two*(m57_6(1,j,k)+m57_6(2,j,k))/xn)
      msq(n67_5,j,k)=msq(n67_5,j,k)
     & +sub67_5(qq)
     & *((xn+1._dp/xn)*m67_5(0,j,k)+two*(m67_5(1,j,k)+m67_5(2,j,k))/xn)
      msq(n17_5,j,k)=msq(n17_5,j,k)
     & +sub17_5(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
      msq(n17_6,j,k)=msq(n17_6,j,k)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
      msq(n27_6,j,k)=msq(n27_6,j,k)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
      msq(n27_5,j,k)=msq(n27_5,j,k)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
      msq(n17_2,j,k)=msq(n17_2,j,k)
     & +sub17_2(qq)
     & *((xn+1._dp/xn)*m17_2(0,j,k)+two*(m17_2(1,j,k)+m17_2(2,j,k))/xn)
      msq(n27_1,j,k)=msq(n27_1,j,k)
     & +sub27_1(qq)
     & *((xn+1._dp/xn)*m27_1(0,j,k)+two*(m27_1(1,j,k)+m27_1(2,j,k))/xn)
      else
c---QbarQbar - process db -> bu (obtained from above by switching p1 and p2)
      msq(n57_6,j,k)=msq(n57_6,j,k)
     & +sub57_6(qq)
     & *((xn+1._dp/xn)*m57_6(0,j,k)+two*(m57_6(2,j,k)+m57_6(1,j,k))/xn)
      msq(n67_5,j,k)=msq(n67_5,j,k)
     & +sub67_5(qq)
     & *((xn+1._dp/xn)*m67_5(0,j,k)+two*(m67_5(2,j,k)+m67_5(1,j,k))/xn)
      msq(n27_5,j,k)=msq(n27_5,j,k)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(1,j,k)-m27_5(0,j,k)/xn-m27_5(2,j,k)/xn)
      msq(n27_6,j,k)=msq(n27_6,j,k)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(2,j,k)-m27_6(0,j,k)/xn-m27_6(1,j,k)/xn)
      msq(n17_6,j,k)=msq(n17_6,j,k)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(1,j,k)-m17_6(0,j,k)/xn-m17_6(2,j,k)/xn)
      msq(n17_5,j,k)=msq(n17_5,j,k)
     & +sub17_5(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*m17_5(2,j,k)-m17_5(0,j,k)/xn-m17_5(1,j,k)/xn)
      msq(n17_2,j,k)=msq(n17_2,j,k)
     & +sub17_2(qq)
     & *((xn+1._dp/xn)*m17_2(0,j,k)+two*(m17_2(2,j,k)+m17_2(1,j,k))/xn)
      msq(n27_1,j,k)=msq(n27_1,j,k)
     & +sub27_1(qq)
     & *((xn+1._dp/xn)*m27_1(0,j,k)+two*(m27_1(2,j,k)+m27_1(1,j,k))/xn)
      endif
      elseif ((j > 0) .and. (k<0)) then
        if (j == flav) then
c--- QQbar - process bu -> bd (obtained from bu -> bd by switching
c---                               p2 and p6)
        msq(n57_6,j,k)=msq(n57_6,j,k)
     &   +sub57_6(qq)
     &   *((xn-two/xn)*m57_6(2,j,k)-m57_6(0,j,k)/xn-m57_6(1,j,k)/xn)
        msq(n67_5,j,k)=msq(n67_5,j,k)
     &   +sub67_5(qq)
     &   *((xn-two/xn)*m67_5(2,j,k)-m67_5(0,j,k)/xn-m67_5(1,j,k)/xn)
        msq(n17_5,j,k)=msq(n17_5,j,k)
     &   +sub17_5(qq)
     &   *((xn-two/xn)*m17_5(1,j,k)-m17_5(0,j,k)/xn-m17_5(2,j,k)/xn)
     &   +sub57_1(qq)
     &   *((xn-two/xn)*m17_5(1,j,k)-m17_5(0,j,k)/xn-m17_5(2,j,k)/xn)
        msq(n17_6,j,k)=msq(n17_6,j,k)
     &   +sub17_6(qq)
     &   *((xn+1._dp/xn)*m17_6(0,j,k)+two*(m17_6(2,j,k)+m17_6(1,j,k))/xn)
     &   +sub67_1(qq)
     &   *((xn+1._dp/xn)*m17_6(0,j,k)+two*(m17_6(2,j,k)+m17_6(1,j,k))/xn)
        msq(n27_6,j,k)=msq(n27_6,j,k)
     &   +sub27_6(qq)
     &   *((xn-two/xn)*m27_6(1,j,k)-m27_6(0,j,k)/xn-m27_6(2,j,k)/xn)
     &   +sub67_2(qq)
     &   *((xn-two/xn)*m27_6(1,j,k)-m27_6(0,j,k)/xn-m27_6(2,j,k)/xn)
        msq(n27_5,j,k)=msq(n27_5,j,k)
     &   +sub27_5(qq)
     &   *((xn+1._dp/xn)*m27_5(0,j,k)+two*(m27_5(2,j,k)+m27_5(1,j,k))/xn)
     &   +sub57_2(qq)
     &   *((xn+1._dp/xn)*m27_5(0,j,k)+two*(m27_5(2,j,k)+m27_5(1,j,k))/xn)
        msq(n17_2,j,k)=msq(n17_2,j,k)
     &   +sub17_2(qq)
     &   *((xn-two/xn)*m17_2(2,j,k)-m17_2(0,j,k)/xn-m17_2(1,j,k)/xn)
        msq(n27_1,j,k)=msq(n27_1,j,k)
     &   +sub27_1(qq)
     &   *((xn-two/xn)*m27_1(2,j,k)-m27_1(0,j,k)/xn-m27_1(1,j,k)/xn)
        else
c--- QQbar - process ub -> bd (obtained from ub -> bd by switching
c---                               p2 and p5)
      msq(n57_6,j,k)=msq(n57_6,j,k)
     & +sub57_6(qq)
     & *((xn-two/xn)*m57_6(1,j,k)-m57_6(0,j,k)/xn-m57_6(2,j,k)/xn)
      msq(n67_5,j,k)=msq(n67_5,j,k)
     & +sub67_5(qq)
     & *((xn-two/xn)*m67_5(1,j,k)-m67_5(0,j,k)/xn-m67_5(2,j,k)/xn)
      msq(n27_5,j,k)=msq(n27_5,j,k)
     & +sub27_5(qq)
     & *((xn-two/xn)*m27_5(2,j,k)-m27_5(0,j,k)/xn-m27_5(1,j,k)/xn)
     & +sub57_2(qq)
     & *((xn-two/xn)*m27_5(2,j,k)-m27_5(0,j,k)/xn-m27_5(1,j,k)/xn)
      msq(n27_6,j,k)=msq(n27_6,j,k)
     & +sub27_6(qq)
     & *((xn+1._dp/xn)*m27_6(0,j,k)+two*(m27_6(1,j,k)+m27_6(2,j,k))/xn)
     & +sub67_2(qq)
     & *((xn+1._dp/xn)*m27_6(0,j,k)+two*(m27_6(1,j,k)+m27_6(2,j,k))/xn)
      msq(n17_6,j,k)=msq(n17_6,j,k)
     & +sub17_6(qq)
     & *((xn-two/xn)*m17_6(2,j,k)-m17_6(0,j,k)/xn-m17_6(1,j,k)/xn)
     & +sub67_1(qq)
     & *((xn-two/xn)*m17_6(2,j,k)-m17_6(0,j,k)/xn-m17_6(1,j,k)/xn)
      msq(n17_5,j,k)=msq(n17_5,j,k)
     & +sub17_5(qq)
     & *((xn+1._dp/xn)*m17_5(0,j,k)+two*(m17_5(1,j,k)+m17_5(2,j,k))/xn)
     & +sub57_1(qq)
     & *((xn+1._dp/xn)*m17_5(0,j,k)+two*(m17_5(1,j,k)+m17_5(2,j,k))/xn)
      msq(n17_2,j,k)=msq(n17_2,j,k)
     & +sub17_2(qq)
     & *((xn-two/xn)*m17_2(1,j,k)-m17_2(0,j,k)/xn-m17_2(2,j,k)/xn)
      msq(n27_1,j,k)=msq(n27_1,j,k)
     & +sub27_1(qq)
     & *((xn-two/xn)*m27_1(1,j,k)-m27_1(0,j,k)/xn-m27_1(2,j,k)/xn)
        endif
      elseif ((j < 0) .and. (k>0)) then
        if (j == -flav) then
c--- QbarQ - process bu -> bd (obtained from QQbar process ub -> bd
c---                               by switching p1 and p2)
      msq(n57_6,j,k)=msq(n57_6,j,k)
     & +sub57_6(qq)
     & *((xn-two/xn)*m57_6(2,j,k)-m57_6(0,j,k)/xn-m57_6(1,j,k)/xn)
      msq(n67_5,j,k)=msq(n67_5,j,k)
     & +sub67_5(qq)
     & *((xn-two/xn)*m67_5(2,j,k)-m67_5(0,j,k)/xn-m67_5(1,j,k)/xn)
      msq(n17_5,j,k)=msq(n17_5,j,k)
     & +sub17_5(qq)
     & *((xn-two/xn)*m17_5(1,j,k)-m17_5(0,j,k)/xn-m17_5(2,j,k)/xn)
     & +sub57_1(qq)
     & *((xn-two/xn)*m17_5(1,j,k)-m17_5(0,j,k)/xn-m17_5(2,j,k)/xn)
      msq(n17_6,j,k)=msq(n17_6,j,k)
     & +sub17_6(qq)
     & *((xn+1._dp/xn)*m17_6(0,j,k)+two*(m17_6(2,j,k)+m17_6(1,j,k))/xn)
     & +sub67_1(qq)
     & *((xn+1._dp/xn)*m17_6(0,j,k)+two*(m17_6(2,j,k)+m17_6(1,j,k))/xn)
      msq(n27_6,j,k)=msq(n27_6,j,k)
     & +sub27_6(qq)
     & *((xn-two/xn)*m27_6(1,j,k)-m27_6(0,j,k)/xn-m27_6(2,j,k)/xn)
     & +sub67_2(qq)
     & *((xn-two/xn)*m27_6(1,j,k)-m27_6(0,j,k)/xn-m27_6(2,j,k)/xn)
      msq(n27_5,j,k)=msq(n27_5,j,k)
     & +sub27_5(qq)
     & *((xn+1._dp/xn)*m27_5(0,j,k)+two*(m27_5(2,j,k)+m27_5(1,j,k))/xn)
     & +sub57_2(qq)
     & *((xn+1._dp/xn)*m27_5(0,j,k)+two*(m27_5(2,j,k)+m27_5(1,j,k))/xn)
      msq(n17_2,j,k)=msq(n17_2,j,k)
     & +sub17_2(qq)
     & *((xn-two/xn)*m17_2(2,j,k)-m17_2(0,j,k)/xn-m17_2(1,j,k)/xn)
      msq(n27_1,j,k)=msq(n27_1,j,k)
     & +sub27_1(qq)
     & *((xn-two/xn)*m27_1(2,j,k)-m27_1(0,j,k)/xn-m27_1(1,j,k)/xn)
        else
c--- QbarQ - process ub -> bd (obtained from QQbar process bu -> bd
c---                               by switching p1 and p2)
        msq(n57_6,j,k)=msq(n57_6,j,k)
     &   +sub57_6(qq)
     &   *((xn-two/xn)*m57_6(1,j,k)-m57_6(0,j,k)/xn-m57_6(2,j,k)/xn)
        msq(n67_5,j,k)=msq(n67_5,j,k)
     &   +sub67_5(qq)
     &   *((xn-two/xn)*m67_5(1,j,k)-m67_5(0,j,k)/xn-m67_5(2,j,k)/xn)
        msq(n27_5,j,k)=msq(n27_5,j,k)
     &   +sub27_5(qq)
     &   *((xn-two/xn)*m27_5(2,j,k)-m27_5(0,j,k)/xn-m27_5(1,j,k)/xn)
     &   +sub57_2(qq)
     &   *((xn-two/xn)*m27_5(2,j,k)-m27_5(0,j,k)/xn-m27_5(1,j,k)/xn)
        msq(n27_6,j,k)=msq(n27_6,j,k)
     &   +sub27_6(qq)
     &   *((xn+1._dp/xn)*m27_6(0,j,k)+two*(m27_6(1,j,k)+m27_6(2,j,k))/xn)
     &   +sub67_2(qq)
     &   *((xn+1._dp/xn)*m27_6(0,j,k)+two*(m27_6(1,j,k)+m27_6(2,j,k))/xn)
        msq(n17_6,j,k)=msq(n17_6,j,k)
     &   +sub17_6(qq)
     &   *((xn-two/xn)*m17_6(2,j,k)-m17_6(0,j,k)/xn-m17_6(1,j,k)/xn)
     &   +sub67_1(qq)
     &   *((xn-two/xn)*m17_6(2,j,k)-m17_6(0,j,k)/xn-m17_6(1,j,k)/xn)
        msq(n17_5,j,k)=msq(n17_5,j,k)
     &   +sub17_5(qq)
     &   *((xn+1._dp/xn)*m17_5(0,j,k)+two*(m17_5(1,j,k)+m17_5(2,j,k))/xn)
     &   +sub57_1(qq)
     &   *((xn+1._dp/xn)*m17_5(0,j,k)+two*(m17_5(1,j,k)+m17_5(2,j,k))/xn)
        msq(n17_2,j,k)=msq(n17_2,j,k)
     &   +sub17_2(qq)
     &   *((xn-two/xn)*m17_2(1,j,k)-m17_2(0,j,k)/xn-m17_2(2,j,k)/xn)
        msq(n27_1,j,k)=msq(n27_1,j,k)
     &   +sub27_1(qq)
     &   *((xn-two/xn)*m27_1(1,j,k)-m27_1(0,j,k)/xn-m27_1(2,j,k)/xn)
        endif

      elseif (j==0) then
c-- G-Q
       if     (k > 0) then
       if (k == +flav) then
       msq(n16_2,j,k)=msq(n16_2,j,k)+sub16_2(qg)*(aveqg/aveqq)*(
     &   +(xn-one/xn)*(msq16_2(-1,k)+msq16_2(-2,k)
     &                +msq16_2(-3,k)+msq16_2(-4,k)))
       msq(n17_2,j,k)=msq(n17_2,j,k)+sub17_2(qg)*(aveqg/aveqq)*(
     &   +(xn-one/xn)*(msq17_2(+1,k)+msq17_2(+2,k)
     &                +msq17_2(+3,k)+msq17_2(+4,k)))
       else
       endif
C-- G-Qbar
        elseif (k < 0) then
        if (k == -flav) then
        msq(n16_2,j,k)=msq(n16_2,j,k)+sub16_2(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*(msq16_2(+1,k)+msq16_2(+2,k)
     &                 +msq16_2(+3,k)+msq16_2(+4,k)))
        msq(n17_2,j,k)=msq(n17_2,j,k)+sub17_2(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*(msq17_2(-1,k)+msq17_2(-2,k)
     &                 +msq17_2(-3,k)+msq17_2(-4,k)))
        else
        endif
        endif
      elseif (k==0) then
c-- Q-G
        if     (j > 0) then
        if ( j == +flav) then
        msq(n26_1,j,k)=msq(n26_1,j,k)+sub26_1(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*(msq26_1(j,-1)+msq26_1(j,-2)
     &                 +msq26_1(j,-3)+msq26_1(j,-4)))
        msq(n27_1,j,k)=msq(n27_1,j,k)+sub27_1(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*(msq27_1(j,+1)+msq27_1(j,+2)
     &                 +msq27_1(j,+3)+msq27_1(j,+4)))
        endif
C-- Qbar-G
        elseif (j < 0) then
        if ( j == -flav) then
        msq(n26_1,j,k)=msq(n26_1,j,k)+sub26_1(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*(msq26_1(j,+1)+msq26_1(j,+2)
     &                 +msq26_1(j,+3)+msq26_1(j,+4)))
        msq(n27_1,j,k)=msq(n27_1,j,k)+sub27_1(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*(msq27_1(j,-1)+msq27_1(j,-2)
     &                 +msq27_1(j,-3)+msq27_1(j,-4)))
        endif
        endif
      endif

      endif

 97   continue
c--- contributions with 2 b-quarks
      if (isub == 2) then
      if (j==0) then
c-- G-Q
       if (kk(k) == 2) then
        msq(n16_2,j,k)=msq(n16_2,j,k)+sub16_2(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*msq16_2(-flav,k))
        msq(n15_2,j,k)=msq(n15_2,j,k)+sub15_2(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*msq15_2(+flav,k))
C-- G-Qbar
        elseif (kk(k) == -1) then
        msq(n16_2,j,k)=msq(n16_2,j,k)+sub16_2(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*msq16_2(-flav,k))
        msq(n15_2,j,k)=msq(n15_2,j,k)+sub15_2(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*msq15_2(+flav,k))
        endif
      elseif (k==0) then
c-- Q-G
        if     (jj(j) == 2) then
        msq(n26_1,j,k)=msq(n26_1,j,k)+sub26_1(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*msq26_1(j,+flav))
        msq(n25_1,j,k)=msq(n25_1,j,k)+sub25_1(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*msq25_1(j,-flav))
C-- Qbar-G
        elseif (jj(j) == -1) then
        msq(n26_1,j,k)=msq(n26_1,j,k)+sub26_1(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*msq26_1(j,+flav))
        msq(n25_1,j,k)=msq(n25_1,j,k)+sub25_1(qg)*(aveqg/aveqq)*(
     &    +(xn-one/xn)*msq25_1(j,-flav))
        endif
      endif

      endif

      enddo
      enddo

      return
      end

