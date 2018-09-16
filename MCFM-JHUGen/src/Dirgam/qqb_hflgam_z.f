      subroutine qqb_hflgam_z(p,z)
      implicit none
      include 'types.f'
************************************************************************
*     Authors: R.K. Ellis and John M. Campbell                         *
*     October, 2002.                                                   *
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
      real(dp):: z,xl12,xl14,xl24,p(mxpart,4),dot
      real(dp):: ii_qq,ii_qg,ii_gq,ii_gg,
     &                 if_qq,if_gg,
     &                 fi_qq,fi_gg

      xl12=log(+two*dot(p,1,2)/musq)
      xl14=log(-two*dot(p,1,4)/musq)
      xl24=log(-two*dot(p,2,4)/musq)

c--- sum over regular and plus terms
      do is=1,3
c--- (q,qb) terms
c      Q1(q,q,a,is) =ason4pi*xn*(if_qq(z,xl14,is)+0.5_dp*fi_gg(z,xl14,is)
c     &                            -ii_qq(z,xl12,is)/xnsq)
c      Q1(a,a,q,is)=Q1(q,q,a,is)
c      Q2(a,a,q,is)=ason4pi*xn*(if_qq(z,xl24,is)+0.5_dp*fi_gg(z,xl24,is)
c     &                            -ii_qq(z,xl12,is)/xnsq)
c      Q2(q,q,a,is) =Q2(a,a,q,is)

c--- (q,g)
      Q2(g,g,q,is)=ason4pi*xn
     & *(ii_gg(z,xl12,is)+if_gg(z,xl24,is)+fi_qq(z,xl24,is))
      Q1(q,q,g,is)=ason4pi*xn*(ii_qq(z,xl12,is)
     &                       -(if_qq(z,xl14,is)+fi_qq(z,xl14,is))/xnsq)
c      Q2(a,g,q,is)=ason4pi*2._dp*tr*ii_qg(z,xl12,is)

c--- Removed qb-g (need to do anti-heavy flavor later)
cc--- (qb,g)
c      Q2(g,g,a,is)=Q2(g,g,q,is)
c      Q1(a,a,g,is)=Q1(q,q,g,is)
cc      Q2(q,g,a,is)=Q2(a,g,q,is)

c--- (g,q)
      Q1(g,g,q,is)=ason4pi*xn
     & *(ii_gg(z,xl12,is)+if_gg(z,xl14,is)+fi_qq(z,xl14,is))
      Q2(q,q,g,is)=ason4pi*xn*(ii_qq(z,xl12,is)
     &                       -(if_qq(z,xl24,is)+fi_qq(z,xl24,is))/xnsq)
c      Q1(a,g,q,is)=ason4pi*2._dp*tr*ii_qg(z,xl12,is)

c--- Removed g-qb (need to do anti-heavy flavor later)
cc--- (g,qb)
c      Q1(g,g,a,is)=Q1(g,g,q,is)
c      Q2(a,a,g,is)=Q2(q,q,g,is)
cc      Q1(q,g,a,is)=Q1(a,g,q,is)

c--- (g,g)
      Q1(q,g,g,is)=ason4pi*2._dp*tr*ii_qg(z,xl12,is)
      Q1(a,g,g,is)=Q1(q,g,g,is)
      Q2(q,g,g,is)=Q1(q,g,g,is)
      Q2(a,g,g,is)=Q1(q,g,g,is)
      
      enddo

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
