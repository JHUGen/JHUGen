      subroutine qqb_wbbm_z(p,z)
      implicit none
      include 'types.f'
************************************************************************
*     Virtual ct's for Wbb with massive b quark                        *
*                                                                      *
*     q(-p1)+qbar(-p2) --> b(p5) + b~(p6) + W +g(p7)                   *
*                                           |                          *
*                                           --> nu(p3)+e^+(p4)         *
*                                                                      *
*     Author: J. Campbell, October 2010                                *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      include 'masses.f'
      include 'heavyflav.f'
      integer:: is,nu
      real(dp):: z,p(mxpart,4),dot,metric,Q56sq
      real(dp):: xl12,xl15,xl16,xl25,xl26,xl56,
     & mbar15,mbar16,mbar25,mbar26,mbar56,mq
      real(dp):: ii_qq,ii_qg,if_mqq,fi_mqq,ff_mqq,tempqg

C--- set up the correct mass, according to 'flav'
      if     (flav == 6) then
        mq=mt
      elseif (flav == 5) then
        mq=mb
      elseif (flav == 4) then
        mq=mc
      else
        write(6,*) 'Wrong flavour in qqb_wbbm_z.f: flav=',flav
        call flush(6)
        stop
      endif

      Q56sq=zero
      metric=one
      do nu=4,1,-1
      Q56sq=Q56sq+metric*(p(5,nu)+p(6,nu))*(p(5,nu)+p(6,nu))
      metric=-one
      enddo

CDTS (5.45,5.77)
      mbar15=mq/sqrt(-two*dot(p,1,5))
      mbar16=mq/sqrt(-two*dot(p,1,6))
      mbar25=mq/sqrt(-two*dot(p,2,5))
      mbar26=mq/sqrt(-two*dot(p,2,6))
CDTS (5.5)
      mbar56=mq/sqrt(Q56sq)

      xl12=log(+two*dot(p,1,2)/musq)
      xl15=log(-two*dot(p,1,5)/musq)
      xl16=log(-two*dot(p,1,6)/musq)
      xl25=log(-two*dot(p,2,5)/musq)
      xl26=log(-two*dot(p,2,6)/musq)
      xl56=log(Q56sq/musq)

      do is=1,3
      Q1(q,q,a,is)=ason4pi
     & *((xn-two/xn)*(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
     &      +two/xn *(if_mqq(z,xl16,mbar16,is)+fi_mqq(z,xl16,mbar16,is))
     &      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
      Q1(a,a,q,is)=ason4pi
     & *((xn-two/xn)*(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
     &      +two/xn *(if_mqq(z,xl26,mbar26,is)+fi_mqq(z,xl26,mbar26,is))
     &      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
      Q2(a,a,q,is)=ason4pi
     & *((xn-two/xn)*(if_mqq(z,xl26,mbar26,is)+fi_mqq(z,xl26,mbar26,is))
     &      +two/xn *(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
     &      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
      Q2(q,q,a,is)=ason4pi
     & *((xn-two/xn)*(if_mqq(z,xl16,mbar16,is)+fi_mqq(z,xl16,mbar16,is))
     &      +two/xn *(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
     &      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))

c--- next four definitions are LC only
c      Q1(q,q,a,is)=ason4pi
c     & *((xn-zip/xn)*(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
c     &      +zip/xn *(if_mqq(z,xl16,mbar16,is)+fi_mqq(z,xl16,mbar16,is))
c     &      -zip/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
c      Q1(a,a,q,is)=ason4pi
c     & *((xn-zip/xn)*(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
c     &      +zip/xn *(if_mqq(z,xl26,mbar26,is)+fi_mqq(z,xl26,mbar26,is))
c     &      -zip/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
c      Q2(a,a,q,is)=ason4pi
c     & *((xn-zip/xn)*(if_mqq(z,xl26,mbar26,is)+fi_mqq(z,xl26,mbar26,is))
c     &      +zip/xn *(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
c     &      -zip/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
c      Q2(q,q,a,is)=ason4pi
c     & *((xn-zip/xn)*(if_mqq(z,xl16,mbar16,is)+fi_mqq(z,xl16,mbar16,is))
c     &      +zip/xn *(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
c     &      -zip/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))

c--- next four definitions are SLC only
c      Q1(q,q,a,is)=ason4pi
c     & *((  -two/xn)*(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
c     &      +two/xn *(if_mqq(z,xl16,mbar16,is)+fi_mqq(z,xl16,mbar16,is))
c     &      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
c      Q1(a,a,q,is)=ason4pi
c     & *((  -two/xn)*(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
c     &      +two/xn *(if_mqq(z,xl26,mbar26,is)+fi_mqq(z,xl26,mbar26,is))
c     &      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
c      Q2(a,a,q,is)=ason4pi
c     & *((  -two/xn)*(if_mqq(z,xl26,mbar26,is)+fi_mqq(z,xl26,mbar26,is))
c     &      +two/xn *(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
c     &      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
c      Q2(q,q,a,is)=ason4pi
c     & *((  -two/xn)*(if_mqq(z,xl16,mbar16,is)+fi_mqq(z,xl16,mbar16,is))
c     &      +two/xn *(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
c     &      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))

      tempqg=ason2pi*tr*ii_qg(z,xl12,is)
      Q1(a,g,q,is)=tempqg
      Q1(q,g,a,is)=tempqg
      Q2(a,g,q,is)=tempqg
      Q2(q,g,a,is)=tempqg
      enddo

      return
      end

