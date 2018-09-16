      subroutine gg_hgg_z(p,z)
      implicit none
      include 'types.f'
************************************************************************
*     Author: John M. Campbell                                         *
*     August, 2005.                                                    *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_h2j.f'
      include 'agq.f'
      integer:: is
      real(dp):: z,xl12,xl15,xl25,xl16,xl26,xl56,p(mxpart,4),dot
      real(dp):: ii_qq,ii_qg,ii_gq,ii_gg,
     &                 if_qq,if_gg,
     &                 fi_qq,fi_gg,ff_qq,ff_gg
c--- pointers from msq_struc:
      integer:: igg_ab,igg_ba,igg_sym,iqq_a,iqq_b,iqq_i,
     & igggg_a,igggg_b,igggg_c,iqr
      parameter(igg_ab=4,igg_ba=5,igg_sym=6)
c--- Note that the 4-quark and 4-gluon pieces are never needed simultaneously,
c--- so we can reuse the same indices to save memory
      parameter(iqq_a=1,iqq_b=2,iqq_i=3)
      parameter(igggg_a=1,igggg_b=2,igggg_c=3)
c--- One extra parameter needed for non-identical quark pieces
      parameter(iqr=7)


      xl12=log(+two*dot(p,1,2)/musq)
      xl15=log(-two*dot(p,1,5)/musq)
      xl16=log(-two*dot(p,1,6)/musq)
      xl25=log(-two*dot(p,2,5)/musq)
      xl26=log(-two*dot(p,2,6)/musq)
      xl56=log(+two*dot(p,5,6)/musq)

c--- sum over regular and plus terms
      do is=1,3

c--- quark-quark terms (also used for antiquark-antiquark)
c--- only the iqq_a terms are used for non-identical quarks
      H1(q,q,q,iqq_a,is)=ason4pi*(
     &      +two/xn*(ii_qq(z,xl12,is)+ff_qq(z,xl56,is))
     &      -one/xn*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
     & +(xn-two/xn)*(if_qq(z,xl16,is)+fi_qq(z,xl16,is)))
      H1(q,q,q,iqq_b,is)=ason4pi*(
     &      +two/xn*(ii_qq(z,xl12,is)+ff_qq(z,xl56,is))
     & +(xn-two/xn)*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
     &      -one/xn*(if_qq(z,xl16,is)+fi_qq(z,xl16,is)))
      H1(q,q,q,iqq_i,is)=ason4pi*(
     & +(xn+one/xn)*(ii_qq(z,xl12,is)+ff_qq(z,xl56,is))
     &      -one/xn*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
     &      -one/xn*(if_qq(z,xl16,is)+fi_qq(z,xl16,is)))

      H2(q,q,q,iqq_a,is)=ason4pi*(
     &      +two/xn*(ii_qq(z,xl12,is)+ff_qq(z,xl56,is))
     & +(xn-two/xn)*(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
     &      -one/xn*(if_qq(z,xl26,is)+fi_qq(z,xl26,is)))
      H2(q,q,q,iqq_b,is)=ason4pi*(
     &      +two/xn*(ii_qq(z,xl12,is)+ff_qq(z,xl56,is))
     &      -one/xn*(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
     & +(xn-two/xn)*(if_qq(z,xl26,is)+fi_qq(z,xl26,is)))
      H2(q,q,q,iqq_i,is)=ason4pi*(
     & +(xn+one/xn)*(ii_qq(z,xl12,is)+ff_qq(z,xl56,is))
     &      -one/xn*(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
     &      -one/xn*(if_qq(z,xl26,is)+fi_qq(z,xl26,is)))

      H1(g,q,q,igg_ab,is)=ason4pi*two*cf*ii_gq(z,xl12,is)
      H1(g,q,q,igg_ba,is)=H1(g,q,q,igg_ab,is)
      H1(g,q,q,igg_sym,is)=H1(g,q,q,igg_ab,is)

      H2(g,q,q,igg_ab,is)=ason4pi*two*cf*ii_gq(z,xl12,is)
      H2(g,q,q,igg_ba,is)=H2(g,q,q,igg_ab,is)
      H2(g,q,q,igg_sym,is)=H2(g,q,q,igg_ab,is)

c--- quark-antiquark terms (also used for antiquark-quark)
c---    4-quark terms
c--- only the iqq_a terms are used for non-identical quarks
      H1(q,q,a,iqq_a,is)=ason4pi*(
     & +(xn-two/xn)*(ii_qq(z,xl12,is)+ff_qq(z,xl56,is))
     &      -one/xn*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
     &      +two/xn*(if_qq(z,xl16,is)+fi_qq(z,xl16,is)))

      H1(q,q,a,iqq_b,is)=ason4pi*(
     &      -one/xn*(ii_qq(z,xl12,is)+ff_qq(z,xl56,is))
     & +(xn-two/xn)*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
     &      +two/xn*(if_qq(z,xl16,is)+fi_qq(z,xl16,is)))

      H1(q,q,a,iqq_i,is)=ason4pi*(
     &      -one/xn*(ii_qq(z,xl12,is)+ff_qq(z,xl56,is))
     &      -one/xn*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
     & +(xn+one/xn)*(if_qq(z,xl16,is)+fi_qq(z,xl16,is)))

      H2(a,a,q,iqq_a,is)=ason4pi*(
     & +(xn-two/xn)*(ii_qq(z,xl12,is)+ff_qq(z,xl56,is))
     &      -one/xn*(if_qq(z,xl26,is)+fi_qq(z,xl26,is))
     &      +two/xn*(if_qq(z,xl25,is)+fi_qq(z,xl25,is)))

      H2(a,a,q,iqq_b,is)=ason4pi*(
     &      -one/xn*(ii_qq(z,xl12,is)+ff_qq(z,xl56,is))
     & +(xn-two/xn)*(if_qq(z,xl26,is)+fi_qq(z,xl26,is))
     &      +two/xn*(if_qq(z,xl25,is)+fi_qq(z,xl25,is)))

      H2(a,a,q,iqq_i,is)=ason4pi*(
     &      -one/xn*(ii_qq(z,xl12,is)+ff_qq(z,xl56,is))
     &      -one/xn*(if_qq(z,xl26,is)+fi_qq(z,xl26,is))
     & +(xn+one/xn)*(if_qq(z,xl25,is)+fi_qq(z,xl25,is)))

      H1(q,q,a,iqr,is)=H1(q,q,a,iqq_b,is)
      H2(a,a,q,iqr,is)=H2(a,a,q,iqq_b,is)

      H1(g,q,a,igg_ab,is)=ason4pi*two*cf*ii_gq(z,xl12,is)
      H1(g,q,a,igg_ba,is)=H1(g,q,a,igg_ab,is)
      H1(g,q,a,igg_sym,is)=H1(g,q,a,igg_ab,is)

      H2(g,a,q,igg_ab,is)=ason4pi*two*cf*ii_gq(z,xl12,is)
      H2(g,a,q,igg_ba,is)=H2(g,a,q,igg_ab,is)
      H2(g,a,q,igg_sym,is)=H2(g,a,q,igg_ab,is)

c---    2-quark terms
      H1(q,q,a,igg_ab,is)=ason4pi*xn*(
     & +ff_gg(z,xl56,is)/2._dp
     & +if_qq(z,xl15,is)+fi_gg(z,xl15,is)/2._dp
     & -ii_qq(z,xl12,is)/xn**2)

      H1(q,q,a,igg_ba,is)=ason4pi*xn*(
     & +ff_gg(z,xl56,is)/2._dp
     & +if_qq(z,xl25,is)+fi_gg(z,xl25,is)/2._dp
     & -ii_qq(z,xl12,is)/xn**2)

      H1(q,q,a,igg_sym,is)=ason4pi*xn*(
     & +if_qq(z,xl25,is)+fi_gg(z,xl25,is)/2._dp
     & +if_qq(z,xl15,is)+fi_gg(z,xl15,is)/2._dp
     & -ii_qq(z,xl12,is)*(1._dp+1._dp/xn**2))

      H2(a,a,q,igg_ab,is)=ason4pi*xn*(
     & +ff_gg(z,xl56,is)/2._dp
     & +if_qq(z,xl26,is)+fi_gg(z,xl26,is)/2._dp
     & -ii_qq(z,xl12,is)/xn**2)

      H2(a,a,q,igg_ba,is)=ason4pi*xn*(
     & +ff_gg(z,xl56,is)/2._dp
     & +if_qq(z,xl16,is)+fi_gg(z,xl16,is)/2._dp
     & -ii_qq(z,xl12,is)/xn**2)

      H2(a,a,q,igg_sym,is)=ason4pi*xn*(
     & +if_qq(z,xl26,is)+fi_gg(z,xl26,is)/2._dp
     & +if_qq(z,xl16,is)+fi_gg(z,xl16,is)/2._dp
     & -ii_qq(z,xl12,is)*(1._dp+1._dp/xn**2))


c--- gluon-gluon terms
c---    0-quark terms
      H1(g,g,g,igggg_a,is)=ason4pi*xn*(
     & +(ii_gg(z,xl12,is)+ff_gg(z,xl56,is)/2._dp)
     & +(if_gg(z,xl16,is)+fi_gg(z,xl16,is)/2._dp))
      H1(g,g,g,igggg_b,is)=ason4pi*xn*(
     & +(if_gg(z,xl16,is)+fi_gg(z,xl16,is)/2._dp)
     & +(if_gg(z,xl15,is)+fi_gg(z,xl15,is)/2._dp))
      H1(g,g,g,igggg_c,is)=ason4pi*xn*(
     & +(ii_gg(z,xl12,is)+ff_gg(z,xl56,is)/2._dp)
     & +(if_gg(z,xl15,is)+fi_gg(z,xl15,is)/2._dp))

      H2(g,g,g,igggg_a,is)=ason4pi*xn*(
     & +(ii_gg(z,xl12,is)+ff_gg(z,xl56,is)/2._dp)
     & +(if_gg(z,xl25,is)+fi_gg(z,xl25,is)/2._dp))
      H2(g,g,g,igggg_b,is)=ason4pi*xn*(
     & +(if_gg(z,xl25,is)+fi_gg(z,xl25,is)/2._dp)
     & +(if_gg(z,xl26,is)+fi_gg(z,xl26,is)/2._dp))
      H2(g,g,g,igggg_c,is)=ason4pi*xn*(
     & +(ii_gg(z,xl12,is)+ff_gg(z,xl56,is)/2._dp)
     & +(if_gg(z,xl26,is)+fi_gg(z,xl26,is)/2._dp))

c---    2-quark terms
      H1(g,g,g,igg_ab,is)=ason4pi*xn*(
     & +ii_gg(z,xl12,is)
     & +if_gg(z,xl16,is)+fi_qq(z,xl16,is)
     & -ff_qq(z,xl56,is)/xn**2)

      H1(g,g,g,igg_ba,is)=ason4pi*xn*(
     & +ii_gg(z,xl12,is)
     & +if_gg(z,xl15,is)+fi_qq(z,xl15,is)
     & -ff_qq(z,xl56,is)/xn**2)

      H1(g,g,g,igg_sym,is)=ason4pi*xn*(
     & +if_gg(z,xl15,is)+fi_qq(z,xl15,is)
     & +if_gg(z,xl16,is)+fi_qq(z,xl16,is)
     & -ff_qq(z,xl56,is)*(1._dp+1._dp/xn**2))

      H2(g,g,g,igg_ab,is)=ason4pi*xn*(
     & +ii_gg(z,xl12,is)
     & +if_gg(z,xl25,is)+fi_qq(z,xl25,is)
     & -ff_qq(z,xl56,is)/xn**2)

      H2(g,g,g,igg_ba,is)=ason4pi*xn*(
     & +ii_gg(z,xl12,is)
     & +if_gg(z,xl26,is)+fi_qq(z,xl26,is)
     & -ff_qq(z,xl56,is)/xn**2)

      H2(g,g,g,igg_sym,is)=ason4pi*xn*(
     & +if_gg(z,xl25,is)+fi_qq(z,xl25,is)
     & +if_gg(z,xl26,is)+fi_qq(z,xl26,is)
     & -ff_qq(z,xl56,is)*(1._dp+1._dp/xn**2))

      H1(q,g,g,igg_ab,is)=ason4pi*ii_qg(z,xl12,is)
      H1(q,g,g,igg_ba,is)=H1(q,g,g,igg_ab,is)
      H1(q,g,g,igg_sym,is)=H1(q,g,g,igg_ab,is)

      H2(q,g,g,igg_ab,is)=ason4pi*ii_qg(z,xl12,is)
      H2(q,g,g,igg_ba,is)=H2(q,g,g,igg_ab,is)
      H2(q,g,g,igg_sym,is)=H2(q,g,g,igg_ab,is)

c--- quark-gluon terms (also used for antiquark-gluon)
      H1(q,q,g,igg_ba,is)=ason4pi*xn*(
     & +if_qq(z,xl16,is)+fi_gg(z,xl16,is)/2._dp
     & -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))/xn**2)

      H1(q,q,g,igg_ab,is)=ason4pi*xn*(
     & +ii_qq(z,xl12,is)+ff_qq(z,xl56,is)
     & -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))/xn**2)

      H1(q,q,g,igg_sym,is)=ason4pi*xn*(
     & +ii_qq(z,xl12,is)+ff_qq(z,xl56,is)
     & +if_qq(z,xl16,is)+fi_gg(z,xl16,is)/2._dp
     & -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*(1._dp+1._dp/xn**2))

      H2(g,g,q,igg_ba,is)=ason4pi*xn*(
     & +if_gg(z,xl26,is)+fi_gg(z,xl26,is)/2._dp
     & +if_gg(z,xl25,is)+fi_qq(z,xl25,is))

      H2(g,g,q,igg_ab,is)=ason4pi*xn*(
     & +if_gg(z,xl26,is)+fi_gg(z,xl26,is)/2._dp
     & +ii_gg(z,xl12,is)+ff_gg(z,xl56,is)/2._dp)

      H2(g,g,q,igg_sym,is)=ason4pi*xn*(
     & +if_gg(z,xl25,is)+fi_qq(z,xl25,is)
     & +ii_gg(z,xl12,is)+ff_gg(z,xl56,is)/2._dp)

      H1(g,q,g,igg_ab,is)=ason4pi*(aveqg/avegg)*ii_gq(z,xl12,is)
      H1(g,q,g,igg_ba,is)=H1(g,q,g,igg_ab,is)
      H1(g,q,g,igg_sym,is)=H1(g,q,g,igg_ab,is)

      H1(g,q,g,igggg_a,is)=ason4pi*(aveqg/avegg)*ii_gq(z,xl12,is)
      H1(g,q,g,igggg_b,is)=H1(g,q,g,igggg_a,is)
      H1(g,q,g,igggg_c,is)=H1(g,q,g,igggg_a,is)

      H2(a,g,q,igg_ab,is)=ason4pi*ii_qg(z,xl12,is)
      H2(a,g,q,igg_ba,is)=H2(a,g,q,igg_ab,is)
      H2(a,g,q,igg_sym,is)=H2(a,g,q,igg_ab,is)

      H2(q,g,q,iqq_a,is)=ason4pi*ii_qg(z,xl12,is)
      H2(q,g,q,iqq_b,is)=H2(q,g,q,iqq_a,is)
      H2(q,g,q,iqq_i,is)=H2(q,g,q,iqq_a,is)

      H2(a,g,q,iqq_a,is)=ason4pi*ii_qg(z,xl12,is)
      H2(a,g,q,iqq_b,is)=H2(a,g,q,iqq_a,is)
      H2(a,g,q,iqq_i,is)=H2(a,g,q,iqq_a,is)

      H2(a,g,q,iqr,is)=ason4pi*ii_qg(z,xl12,is)

c--- gluon-quark terms (also used for gluon-antiquark)
      H1(g,g,q,igg_ba,is)=ason4pi*xn*(
     & +if_gg(z,xl16,is)+fi_gg(z,xl16,is)/2._dp
     & +if_gg(z,xl15,is)+fi_qq(z,xl15,is))

      H1(g,g,q,igg_ab,is)=ason4pi*xn*(
     & +if_gg(z,xl16,is)+fi_gg(z,xl16,is)/2._dp
     & +ii_gg(z,xl12,is)+ff_gg(z,xl56,is)/2._dp)

      H1(g,g,q,igg_sym,is)=ason4pi*xn*(
     & +if_gg(z,xl15,is)+fi_qq(z,xl15,is)
     & +ii_gg(z,xl12,is)+ff_gg(z,xl56,is)/2._dp)

      H2(q,q,g,igg_ba,is)=ason4pi*xn*(
     & +if_qq(z,xl26,is)+fi_gg(z,xl26,is)/2._dp
     & -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))/xn**2)

      H2(q,q,g,igg_ab,is)=ason4pi*xn*(
     & +ii_qq(z,xl12,is)+ff_qq(z,xl56,is)
     & -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))/xn**2)

      H2(q,q,g,igg_sym,is)=ason4pi*xn*(
     & +ii_qq(z,xl12,is)+ff_qq(z,xl56,is)
     & +if_qq(z,xl26,is)+fi_gg(z,xl26,is)/2._dp
     & -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*(1._dp+1._dp/xn**2))

      H2(g,q,g,igggg_a,is)=ason4pi*(aveqg/avegg)*ii_gq(z,xl12,is)
      H2(g,q,g,igggg_b,is)=H2(g,q,g,igggg_a,is)
      H2(g,q,g,igggg_c,is)=H2(g,q,g,igggg_a,is)

      H2(g,q,g,igg_ab,is)=ason4pi*(aveqg/avegg)*ii_gq(z,xl12,is)
      H2(g,q,g,igg_ba,is)=H2(g,q,g,igg_ab,is)
      H2(g,q,g,igg_sym,is)=H2(g,q,g,igg_ab,is)

      H1(a,g,q,igg_ab,is)=ason4pi*ii_qg(z,xl12,is)
      H1(a,g,q,igg_ba,is)=H1(a,g,q,igg_ab,is)
      H1(a,g,q,igg_sym,is)=H1(a,g,q,igg_ab,is)

      H1(q,g,q,iqq_a,is)=ason4pi*ii_qg(z,xl12,is)
      H1(q,g,q,iqq_b,is)=H1(q,g,q,iqq_a,is)
      H1(q,g,q,iqq_i,is)=H1(q,g,q,iqq_a,is)

      H1(a,g,q,iqq_a,is)=ason4pi*ii_qg(z,xl12,is)
      H1(a,g,q,iqq_b,is)=H1(a,g,q,iqq_a,is)
      H1(a,g,q,iqq_i,is)=H1(a,g,q,iqq_a,is)

      H1(a,g,q,iqr,is)=ason4pi*ii_qg(z,xl12,is)

      enddo

      return
      end
