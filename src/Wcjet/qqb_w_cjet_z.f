      subroutine qqb_w_cjet_z(p,z)
************************************************************************
*     Author: J.M. Campbell                                            *
*     December, 2004.                                                  *
************************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      include 'breit.f'
      integer is
      double precision z,xl12,xl15,xl25,p(mxpart,4),dot
      double precision mbar12,mbar15,mbar25
      double precision ii_mqq,ii_mqg,ii_mgq,ii_mgg,
     .                 if_mqq,if_mgg,
     .                 fi_mqq

      xl12=dlog(+two*dot(p,1,2)/musq)
      xl15=dlog(-two*dot(p,1,5)/musq)
      xl25=dlog(-two*dot(p,2,5)/musq)
      mbar12=0d0
CDTS (5.45,5.77)
      mbar15=mass2/dsqrt(-2d0*dot(p,1,5))
      mbar25=mass2/dsqrt(-2d0*dot(p,2,5))
     
c--- sum over regular and plus terms
      do is=1,3
c--- (q,g)
      Q2(g,g,q,is)=ason4pi*xn
     & *(ii_mgg(z,xl12,mbar12,is)+if_mgg(z,xl25,mbar25,is)
     &  +fi_mqq(z,xl25,mbar25,is))
      Q1(q,q,g,is)=ason4pi*xn
     & *(ii_mqq(z,xl12,mbar12,is)
     & -(if_mqq(z,xl15,mbar15,is)+fi_mqq(z,xl15,mbar15,is))/xnsq)
c--- (qb,g)
      Q2(g,g,a,is)=Q2(g,g,q,is)
      Q1(a,a,g,is)=Q1(q,q,g,is)

c--- (g,q)
      Q1(g,g,q,is)=ason4pi*xn
     & *(ii_mgg(z,xl12,mbar12,is)+if_mgg(z,xl15,mbar15,is)
     &  +fi_mqq(z,xl15,mbar15,is))
      Q2(q,q,g,is)=ason4pi*xn
     & *(ii_mqq(z,xl12,mbar12,is)
     & -(if_mqq(z,xl25,mbar25,is)+fi_mqq(z,xl25,mbar25,is))/xnsq)
c--- (g,qb)
      Q1(g,g,a,is)=Q1(g,g,q,is)
      Q2(a,a,g,is)=Q2(q,q,g,is)

c--- (g,g)
      Q1(q,g,g,is)=ason4pi*2d0*tr*ii_mqg(z,xl12,mbar12,is)
      Q2(q,g,g,is)=Q1(q,g,g,is)
 
c--- these are terms coming from the diagrams which are
c--- a 4-quark contribution, with (q,q~)->g in the initial state
      Q2(g,q,a,is)=ason4pi*2d0*cf*ii_mgq(z,xl12,mbar12,is)
      Q2(g,a,q,is)=Q2(g,q,a,is)
      Q2(g,a,a,is)=Q2(g,q,a,is)
      Q2(g,q,q,is)=Q2(g,q,a,is)

      Q1(g,q,a,is)=Q2(g,q,a,is)
      Q1(g,a,q,is)=Q2(g,q,a,is)
      Q1(g,a,a,is)=Q2(g,q,a,is)
      Q1(g,q,q,is)=Q2(g,q,a,is)
      
      enddo
      
      return
      end
