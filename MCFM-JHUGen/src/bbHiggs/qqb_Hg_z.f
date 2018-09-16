      subroutine qqb_Hg_z(p,z)
      implicit none
      include 'types.f'
************************************************************************
*     Author: John M. Campbell                                         *
*     February, 2002                                                   *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      integer:: is
      real(dp):: z,p(mxpart,4),dot
      real(dp):: ii_qg,ii_gq,fi_qq,
     &                 ii_gg,if_gg,ii_qq,if_qq
      real(dp):: xl12,xl15,xl25

      xl12=log(+two*dot(p,1,2)/musq)
      xl15=log(-two*dot(p,1,5)/musq)
      xl25=log(-two*dot(p,2,5)/musq)

c--- 2-quark terms
c--- sum over regular and plus terms
      do is=1,3
c--- No (q,qb) terms here

c--- (q,g)
      Q2(g,g,q,is)=ason4pi*xn
     & *(ii_gg(z,xl12,is)
     &  +if_gg(z,xl25,is)+fi_qq(z,xl25,is))
      Q1(q,q,g,is)=ason4pi*xn
     & *(ii_qq(z,xl12,is)
     & -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))/xnsq)
c--- (qb,g)
      Q2(g,g,a,is)=Q2(g,g,q,is)
      Q1(a,a,g,is)=Q1(q,q,g,is)

c--- (g,q)
      Q1(g,g,q,is)=ason4pi*xn
     & *(ii_gg(z,xl12,is)
     &  +if_gg(z,xl15,is)+fi_qq(z,xl15,is))
      Q2(q,q,g,is)=ason4pi*xn
     & *(ii_qq(z,xl12,is)
     & -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))/xnsq)
c--- (g,qb)
      Q1(g,g,a,is)=Q1(g,g,q,is)
      Q2(a,a,g,is)=Q2(q,q,g,is)

c--- (g,g)
      Q1(q,g,g,is)=ason4pi*2._dp*tr*ii_qg(z,xl12,is)
      Q1(a,g,g,is)=Q1(q,g,g,is)
      Q2(q,g,g,is)=Q1(q,g,g,is)
      Q2(a,g,g,is)=Q1(q,g,g,is)

      enddo

c--- 4-quark terms
      do is=1,3
      Q1(g,q,q,is)=ason4pi*(xn-1._dp/xn)*ii_gq(z,xl12,is)
      Q2(g,q,q,is)=ason4pi*(xn-1._dp/xn)*ii_gq(z,xl12,is)
      Q1(g,a,a,is)=Q1(g,q,q,is)
      Q2(g,a,a,is)=Q2(g,q,q,is)
      Q1(g,a,q,is)=Q1(g,q,q,is)
      Q2(g,a,q,is)=Q2(g,q,q,is)
      Q1(g,q,a,is)=Q1(g,q,q,is)
      Q2(g,q,a,is)=Q2(g,q,q,is)

      enddo

      return
      end
