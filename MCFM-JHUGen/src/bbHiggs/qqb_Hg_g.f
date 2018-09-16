************************************************************************
*  Normalization checked 2/14/02 by JMC against those provided by      *
*  Fabio and Scott                                                     *
*  ISUB picks out the sub-process as follows:                          *
*    isub = 0 : all of the processes below                             *
*    isub = 1 : g + b -> b + g ,  b + q  -> b + q                      *
*    isub = 2 : g + g -> b~ + b , b + b  -> b + b                      *
************************************************************************
      subroutine qqb_Hg_g(p,msq)
      implicit none
      include 'types.f'

c---Matrix element squared averaged over initial colors and spins
c     parton(-p1)+parton(-p2) --> H(p)+parton(p5)+parton(p6)
c                                  |
c                                   --> b(p3)+bb(p4)
c  via coupling to b
c--all momenta incoming
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'msbarmasses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'susycoup.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'scale.f'
      include 'couple.f'
      integer:: j,k,isub
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,propsq,hdecay,coupsq_eff,ghbb_eff
      real(dp):: mb_eff,massfrun
      real(dp):: bbggh,bbaqh,bbbbh
      common/isub/isub

c--susycoup is the deviation of Higgs coupling
c-- from the standard model value

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call spinoru(6,p,za,zb)
      if (s(3,4) < 4._dp*mbsq) return

c--- run mb to appropriate scale
      mb_eff=massfrun(mb_msbar,scale,amz,2)
c      mb_eff=mb_msbar

      call hbbdecay(p,3,4,hdecay)
      hdecay=hdecay*susycoup**2
      propsq=1._dp/((s(3,4)-hmass**2)**2+(hmass*hwidth)**2)
c--- The _eff couplings include the running mass
c--- We need to separate these from the factors associated with the
c--- Higgs decay, because the Br. Ratio does not include running mb
      ghbb_eff=sqrt(esq/xw)*mb_eff/2._dp/wmass
      coupsq_eff=susycoup**2*ghbb_eff**2

      fac=gsq**2*coupsq_eff*propsq*hdecay

c     BBGGH=     0 -> bbar1 b2 g3 g4 h

      if ((isub == 2) .or. (isub == 0)) then
c --- g + g -> b + b~
      msq(0,0)=avegg*fac*BBGGH(6,5,1,2)
      endif

      if ((isub == 1) .or. (isub == 0)) then
c --- g + b -> b + g
      msq(0,+5)=aveqg*fac*BBGGH(2,5,1,6)
c --- g + b~ -> b~ + g
      msq(0,-5)=aveqg*fac*BBGGH(5,2,1,6)
c --- b + g -> b + g
      msq(+5,0)=aveqg*fac*BBGGH(1,5,2,6)
c --- b~ + g -> b~ + g
      msq(-5,0)=aveqg*fac*BBGGH(5,1,2,6)
      endif

      if ((isub == 1) .or. (isub == 0)) then
c --- q + b-> b + q
      msq(+1,+5)=aveqq*fac*BBAQH(2,5,1,6)
c --- b + q-> b + q
      msq(+5,+1)=aveqq*fac*BBAQH(1,5,2,6)
c --- b~ + q~ -> b~ + q~
      msq(-5,-1)=aveqq*fac*BBAQH(5,1,6,2)
c --- q~ + b~-> b~+ q~
      msq(-1,-5)=aveqq*fac*BBAQH(5,2,6,1)

c --- q~ + b-> b + q~
      msq(-1,+5)=aveqq*fac*BBAQH(2,5,6,1)
c --- b + q~-> b + q~
      msq(+5,-1)=aveqq*fac*BBAQH(1,5,6,2)

c --- q + b~-> b~ + q
      msq(+1,-5)=aveqq*fac*BBAQH(5,2,1,6)
c --- b~ + q -> b~ + q
      msq(-5,+1)=aveqq*fac*BBAQH(5,1,2,6)

      do j=2,4
      msq(+j,+5)=msq(+1,+5)
      msq(-j,+5)=msq(-1,+5)
      msq(+j,-5)=msq(+1,-5)
      msq(-j,-5)=msq(-1,-5)

      msq(+5,+j)=msq(+5,+1)
      msq(+5,-j)=msq(+5,-1)
      msq(-5,+j)=msq(-5,+1)
      msq(-5,-j)=msq(-5,-1)
      enddo
      endif

c 0 -> bbar1 b2 bbar3 b4 h

      if ((isub == 2) .or. (isub == 0)) then
c    BBAQH=     0 -> bbar1 b2 a3 q4 h
c --- q + q~-> b + b~
      msq(+1,-1)=aveqq*fac*BBAQH(6,5,1,2)
c --- q~ + q-> b + b~
      msq(-1,+1)=aveqq*fac*BBAQH(6,5,2,1)
      do j=2,4
      msq(j,-j)=msq(+1,-1)
      msq(-j,j)=msq(-1,+1)
      enddo

      if (s(5,6) < 4._dp*mb**2) return
c --- b + b -> b + b
      msq(+5,+5)=half*fac*aveqq*BBBBH(1,5,2,6)
c --- b~ + b~ -> b~ + b~
      msq(-5,-5)=half*fac*aveqq*BBBBH(5,1,6,2)
c --- b + b~ -> b + b~
      msq(+5,-5)=fac*aveqq*BBBBH(1,5,6,2)
c --- b~ + b -> b + b~
      msq(-5,+5)=fac*aveqq*BBBBH(6,5,2,1)
      endif

      return
      end





