      subroutine qqb_ttw_z(pin,z)
************************************************************************
*     Virtual ct's for ttw with massive t quark                        *
*                                                                      *
*     q(-p1)+qbar(-p2) --> t(p345) + t~(p678) + W                      *
*                                               |                      *
*                                               --> nu(9)+e^+(10)      *
*                                                                      *
*     Author: R.K. Ellis, March 2012                                   *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      include 'masses.f'
      integer is,nu
      double precision z,p(mxpart,4),pin(mxpart,4),dot,Q56sq
      double precision xl12,xl15,xl16,xl25,xl26,xl56,
     & mbar15,mbar16,mbar25,mbar26,mbar56
      double precision ii_qq,ii_qg,if_mqq,fi_mqq,ff_mqq,tempqg

      do nu=1,4
      p(1,nu)=pin(1,nu)
      p(2,nu)=pin(2,nu)
      p(3,nu)=pin(9,nu)
      p(4,nu)=pin(10,nu)
      p(5,nu)=pin(3,nu)+pin(4,nu)+pin(5,nu)
      p(6,nu)=pin(6,nu)+pin(7,nu)+pin(8,nu)
      enddo

      Q56sq=+(p(5,4)+p(6,4))**2
     &      -(p(5,1)+p(6,1))**2
     &      -(p(5,2)+p(6,2))**2
     &      -(p(5,3)+p(6,3))**2

CDTS (5.45,5.77)
      mbar15=mt/dSqrt(-2d0*dot(p,1,5))
      mbar16=mt/dSqrt(-2d0*dot(p,1,6))
      mbar25=mt/dSqrt(-2d0*dot(p,2,5))
      mbar26=mt/dSqrt(-2d0*dot(p,2,6))
CDTS (5.5)
      mbar56=mt/dSqrt(Q56sq)
      
      xl12=dlog(+two*dot(p,1,2)/musq)
      xl15=dlog(-two*dot(p,1,5)/musq)
      xl16=dlog(-two*dot(p,1,6)/musq)
      xl25=dlog(-two*dot(p,2,5)/musq)
      xl26=dlog(-two*dot(p,2,6)/musq)
      xl56=dlog(Q56sq/musq)

      do is=1,3
      Q1(q,q,a,is)=ason4pi
     . *((xn-two/xn)*(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
     .      +two/xn *(if_mqq(z,xl16,mbar16,is)+fi_mqq(z,xl16,mbar16,is))
     .      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
      Q1(a,a,q,is)=ason4pi
     . *((xn-two/xn)*(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
     .      +two/xn *(if_mqq(z,xl26,mbar26,is)+fi_mqq(z,xl26,mbar26,is))
     .      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
      Q2(a,a,q,is)=ason4pi
     . *((xn-two/xn)*(if_mqq(z,xl26,mbar26,is)+fi_mqq(z,xl26,mbar26,is))
     .      +two/xn *(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
     .      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
      Q2(q,q,a,is)=ason4pi
     . *((xn-two/xn)*(if_mqq(z,xl16,mbar16,is)+fi_mqq(z,xl16,mbar16,is))
     .      +two/xn *(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
     .      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))

c--- next four definitions are LC only
c      Q1(q,q,a,is)=ason4pi
c     . *((xn-zip/xn)*(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
c     .      +zip/xn *(if_mqq(z,xl16,mbar16,is)+fi_mqq(z,xl16,mbar16,is))
c     .      -zip/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
c      Q1(a,a,q,is)=ason4pi
c     . *((xn-zip/xn)*(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
c     .      +zip/xn *(if_mqq(z,xl26,mbar26,is)+fi_mqq(z,xl26,mbar26,is))
c     .      -zip/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
c      Q2(a,a,q,is)=ason4pi
c     . *((xn-zip/xn)*(if_mqq(z,xl26,mbar26,is)+fi_mqq(z,xl26,mbar26,is))
c     .      +zip/xn *(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
c     .      -zip/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
c      Q2(q,q,a,is)=ason4pi
c     . *((xn-zip/xn)*(if_mqq(z,xl16,mbar16,is)+fi_mqq(z,xl16,mbar16,is))
c     .      +zip/xn *(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
c     .      -zip/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))

c--- next four definitions are SLC only
c      Q1(q,q,a,is)=ason4pi
c     . *((  -two/xn)*(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
c     .      +two/xn *(if_mqq(z,xl16,mbar16,is)+fi_mqq(z,xl16,mbar16,is))
c     .      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
c      Q1(a,a,q,is)=ason4pi
c     . *((  -two/xn)*(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
c     .      +two/xn *(if_mqq(z,xl26,mbar26,is)+fi_mqq(z,xl26,mbar26,is))
c     .      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
c      Q2(a,a,q,is)=ason4pi
c     . *((  -two/xn)*(if_mqq(z,xl26,mbar26,is)+fi_mqq(z,xl26,mbar26,is))
c     .      +two/xn *(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))
c     .      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))
c      Q2(q,q,a,is)=ason4pi
c     . *((  -two/xn)*(if_mqq(z,xl16,mbar16,is)+fi_mqq(z,xl16,mbar16,is))
c     .      +two/xn *(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))
c     .      -one/xn *(ii_qq(z,xl12,is)+ff_mqq(z,xl56,mbar56,is)))

      tempqg=ason2pi*tr*ii_qg(z,xl12,is)
      Q1(a,g,q,is)=tempqg
      Q1(q,g,a,is)=tempqg
      Q2(a,g,q,is)=tempqg
      Q2(q,g,a,is)=tempqg
      enddo

      return
      end

