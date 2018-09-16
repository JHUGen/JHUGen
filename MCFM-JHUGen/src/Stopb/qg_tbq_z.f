      subroutine qg_tbq_z(p,z)
      implicit none
      include 'types.f'
************************************************************************
*     Virtual ct's for t-channel single top, with explicit b-quark     *
*                                                                      *
*     q(p1) + g(p2) -> t(p3) + b(p4) + q'(p5)                          *
*                                                                      *
*     Author: J. Campbell, February 27, 2008                           *
*                                                                      *
************************************************************************
*                                                                      *
*    Modified so that off-diagonal subtractions are initial-final      *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'agq.f'
      include 'PR_new.f'
      include 'noglue.f'
      include 'stopscales.f'
      include 'breit.f'
      real(dp):: z,p(mxpart,4),dot,metric,Q34sq,
     & xl13,xl14,xl23,xl24,xl34,xl15,xl25,xl15h,xl25h,
     & mbar13,mbar14,mbar23,mbar24,mbar34a,mbar34b,tempgq1,tempgq2,
     & if_mgg,fi_mqq,ff_2mqq,if_qg,if_gq,if_qq,fi_qq,ason4pi_H,ason4pi_L
      integer:: is,nu


      Q34sq=0._dp
      metric=1._dp
      do nu=4,1,-1
      Q34sq=Q34sq+metric*(p(3,nu)+p(4,nu))*(p(3,nu)+p(4,nu))
      metric=-1._dp
      enddo

CDTS (5.45,5.77)
      mbar13=mass2/sqrt(-two*dot(p,1,3))
      mbar14=mass3/sqrt(-two*dot(p,1,4))
      mbar23=mass2/sqrt(-two*dot(p,2,3))
      mbar24=mass3/sqrt(-two*dot(p,2,4))
CDTS (5.5)
      mbar34a=mass2/sqrt(Q34sq)
      mbar34b=mass3/sqrt(Q34sq)

      xl15=log(-two*dot(p,1,5)/renscale_L**2)
      xl25=log(-two*dot(p,2,5)/renscale_L**2)
      xl13=log(-two*dot(p,1,3)/renscale_H**2)
      xl14=log(-two*dot(p,1,4)/renscale_H**2)
      xl23=log(-two*dot(p,2,3)/renscale_H**2)
      xl24=log(-two*dot(p,2,4)/renscale_H**2)
      xl15h=log(-two*dot(p,1,5)/renscale_H**2)
      xl25h=log(-two*dot(p,2,5)/renscale_H**2)
      xl34=log(Q34sq/renscale_H**2)

      ason4pi_H=as_H/fourpi
      ason4pi_L=as_L/fourpi

c--- The variables R and P provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (Q1(a,b,c,is) and leg 2 (Q2(a,b,c,is)
c--- In each case the parton labelling is Using the normal QM notation of putting
c--- everything backward
c---       emitted line after emission =   a
c---       emitter before emission     =   b
c---       spectator                   =   c
c--- There is no label for he or she who is emitted.
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

      do is=1,3

c--- implement shortcut for checking
      if ((ggonly .eqv. .false.) .and. (noglue .eqv. .false.)) then
c--- counterterms for extra gluon emission
      Q1(q,q,g,is)=ason4pi_L*two*cf*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
      Q2(g,g,q,is)=ason4pi_H*(
     &         xn *(if_mgg(z,xl23,mbar23,is)+fi_mqq(z,xl23,mbar23,is))
     &        +xn *(if_mgg(z,xl24,mbar24,is)+fi_mqq(z,xl24,mbar24,is))
     &    -one/xn *(ff_2mqq(z,xl34,mbar34a,mbar34b,is)
     &             +ff_2mqq(z,xl34,mbar34b,mbar34a,is))
     &                      )
      Q1(a,a,g,is)=Q1(q,q,g,is)
      Q2(g,g,a,is)=Q2(g,g,q,is)

      Q2(q,q,g,is)=ason4pi_L*two*cf*(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
      Q1(g,g,q,is)=ason4pi_H*(
     &         xn *(if_mgg(z,xl13,mbar13,is)+fi_mqq(z,xl13,mbar13,is))
     &        +xn *(if_mgg(z,xl14,mbar14,is)+fi_mqq(z,xl14,mbar14,is))
     &    -one/xn *(ff_2mqq(z,xl34,mbar34a,mbar34b,is)
     &             +ff_2mqq(z,xl34,mbar34b,mbar34a,is))
     &                      )
      Q2(a,a,g,is)=Q2(q,q,g,is)
      Q1(g,g,a,is)=Q1(g,g,q,is)
      endif

c--- counterterms for the splitting g->qq~ on the light quark line
      Q1(q,g,g,is)=ason4pi_L*two*tr*if_qg(z,xl15,is)
      Q2(q,g,g,is)=ason4pi_L*two*tr*if_qg(z,xl25,is)

      tempgq1=ason4pi_H*two*cf*if_gq(z,xl15h,is)
      tempgq2=ason4pi_H*two*cf*if_gq(z,xl25h,is)
c--- counterterms for the splitting q->gq on the heavy quark line
      Q2(g,q,q,is)=tempgq2
      Q2(g,a,q,is)=tempgq2
      Q2(g,q,a,is)=tempgq2
      Q2(g,a,a,is)=tempgq2
      Q1(g,q,q,is)=tempgq1
      Q1(g,a,q,is)=tempgq1
      Q1(g,q,a,is)=tempgq1
      Q1(g,a,a,is)=tempgq1

      enddo

      return
      end

