      subroutine gQ_zQ_z(p,z)
************************************************************************
*   Authors: R.K. Ellis and John M. Campbell                           *
*   July, 2003.                                                        *
*   Int. sub. terms for Z + heavy quark (of flavour "flav") production *
************************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      integer is
      double precision z,p(mxpart,4),dot
      double precision ii_qq,ii_gg,ii_qg,ii_gq,
     .                 if_qq,if_gg,
     .                 fi_qq
      double precision xl12,xl15,xl25

      xl12=dlog(+two*dot(p,1,2)/musq)
      xl15=dlog(-two*dot(p,1,5)/musq)
      xl25=dlog(-two*dot(p,2,5)/musq)

c--- sum over regular and plus terms
      do is=1,3
c--- (q,g)
      Q2(g,g,q,is)=ason4pi*xn
     & *(ii_gg(z,xl12,is)
     &  +if_gg(z,xl25,is)+fi_qq(z,xl25,is))
      Q1(q,q,g,is)=ason4pi*xn
     & *(ii_qq(z,xl12,is)
     & -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))/xnsq)
      Q2(a,g,q,is)=ason4pi*2d0*tr*ii_qg(z,xl12,is)
c--- (qb,g)
      Q2(g,g,a,is)=Q2(g,g,q,is)
      Q1(a,a,g,is)=Q1(q,q,g,is)
      Q2(q,g,a,is)=Q2(a,g,q,is)

c--- (g,q)
      Q1(g,g,q,is)=ason4pi*xn
     & *(ii_gg(z,xl12,is)
     &  +if_gg(z,xl15,is)+fi_qq(z,xl15,is))
      Q2(q,q,g,is)=ason4pi*xn
     & *(ii_qq(z,xl12,is)
     & -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))/xnsq)
      Q1(a,g,q,is)=ason4pi*2d0*tr*ii_qg(z,xl12,is)
c--- (g,qb)
      Q1(g,g,a,is)=Q1(g,g,q,is)
      Q2(a,a,g,is)=Q2(q,q,g,is)
      Q1(q,g,a,is)=Q1(a,g,q,is)

c--- (g,g)
      Q1(q,g,g,is)=ason4pi*2d0*tr*ii_qg(z,xl12,is)
      Q1(a,g,g,is)=Q1(q,g,g,is)
      Q2(q,g,g,is)=Q1(q,g,g,is)
      Q2(a,g,g,is)=Q1(q,g,g,is)
      
      enddo

      do is=1,3
      Q1(g,q,q,is)=ason4pi*(xn-1d0/xn)*ii_gq(z,xl12,is)
      Q2(g,q,q,is)=ason4pi*(xn-1d0/xn)*ii_gq(z,xl12,is)
      Q1(g,a,a,is)=Q1(g,q,q,is)
      Q2(g,a,a,is)=Q2(g,q,q,is)
      Q1(g,a,q,is)=Q1(g,q,q,is)
      Q2(g,a,q,is)=Q2(g,q,q,is)
      Q1(g,q,a,is)=Q1(g,q,q,is)
      Q2(g,q,a,is)=Q2(g,q,q,is)

      enddo
      
      return
      end
