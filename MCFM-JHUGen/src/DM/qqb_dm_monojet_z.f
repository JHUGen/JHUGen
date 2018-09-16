      subroutine qqb_dm_monojet_z(p,z)
      implicit none
      include 'types.f'
************************************************************************
!---- Modified for DM by CW Nov 2012 based upon Z1jet routine by 
*     Authors: R.K. Ellis and John M. Campbell                         *
*     November, 2001.                                                  *
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
      real(dp):: z,xl12,xl15,xl25,p(mxpart,4),dot
      real(dp):: ii_qq,ii_qg,ii_gq,ii_gg,
     &                 if_qq,if_gg,
     &                 fi_qq,fi_gg

      xl12=log(+two*dot(p,1,2)/musq)
      xl15=log(-two*dot(p,1,5)/musq)
      xl25=log(-two*dot(p,2,5)/musq)

c--- sum over regular and plus terms
      do is=1,3
c--- (q,qb) terms
      Q1(q,q,a,is) =ason4pi*xn*(if_qq(z,xl15,is)+0.5d0*fi_gg(z,xl15,is)
     &                         -ii_qq(z,xl12,is)/xnsq)
      Q1(a,a,q,is)=Q1(q,q,a,is)
      Q2(a,a,q,is)=ason4pi*xn*(if_qq(z,xl25,is)+0.5d0*fi_gg(z,xl25,is)
     &                         -ii_qq(z,xl12,is)/xnsq)
      Q2(q,q,a,is) =Q2(a,a,q,is)

c--- (q,g)
      Q2(g,g,q,is)=ason4pi*xn
     & *(ii_gg(z,xl12,is)+if_gg(z,xl25,is)+fi_qq(z,xl25,is))
      Q1(q,q,g,is)=ason4pi*xn*(ii_qq(z,xl12,is)
     &                       -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))/xnsq)
      Q2(a,g,q,is)=ason4pi*2d0*tr*ii_qg(z,xl12,is)
c--- (qb,g)
      Q2(g,g,a,is)=Q2(g,g,q,is)
      Q1(a,a,g,is)=Q1(q,q,g,is)
      Q2(q,g,a,is)=Q2(a,g,q,is)

c--- (g,q)
      Q1(g,g,q,is)=ason4pi*xn
     & *(ii_gg(z,xl12,is)+if_gg(z,xl15,is)+fi_qq(z,xl15,is))
      Q2(q,q,g,is)=ason4pi*xn*(ii_qq(z,xl12,is)
     &                       -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))/xnsq)
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
