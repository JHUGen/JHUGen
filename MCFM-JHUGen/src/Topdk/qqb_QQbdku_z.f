      subroutine qqb_QQbdku_z(p,z)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis and John Campbell                             *
*     July, 2008.                                                      *
*     calculate the subtraction terms for the process                  *
*                                                                      *
*     q(-p1) +qbar(-p2)=bbar(p6)+e-(p7)+nubar(p8)+nu(p3)+e+(p4)+b(p5)  *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*                                                                      *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'agq.f'
      include 'PR_cs_new.f'
      include 'breit.f'
      real(dp):: z,p(mxpart,4),pt(mxpart,4),dot,metric,Q34sq,
     & xl12,xl13,xl14,xl23,xl24,xl34,
     & mbar12,mbar13,mbar14,mbar23,mbar24,mbar34,tempgq,tempqg,
     &                 ii_mqq,ii_mgg,
     &                 if_mqq,if_mgg,
     &                 fi_mqq,
     &                 ff_mqq,
     &                 ii_qg,ii_gq
      integer:: is,nu,icol

c--- reduce momenta to top and anti-top sums
      do nu=1,4
      pt(1,nu)=p(1,nu)
      pt(2,nu)=p(2,nu)
      pt(3,nu)=p(3,nu)+p(4,nu)+p(5,nu)
      pt(4,nu)=p(6,nu)+p(7,nu)+p(8,nu)
      enddo

      Q34sq=0._dp
      metric=1._dp
      do nu=4,1,-1
      Q34sq=Q34sq+metric*(pt(3,nu)+pt(4,nu))*(pt(3,nu)+pt(4,nu))
      metric=-1._dp
      enddo

      mbar12=0._dp
CDTS (5.45,5.77)
      mbar13=mass2/sqrt(-two*Dot(pt,1,3))
      mbar23=mass2/sqrt(-two*Dot(pt,2,3))
      mbar14=mass2/sqrt(-two*Dot(pt,1,4))
      mbar24=mass2/sqrt(-two*Dot(pt,2,4))
CDTS (5.5)
      mbar34=mass2/sqrt(Q34sq)

      xl12=log(+two*Dot(pt,1,2)/musq)
      xl13=log(-two*Dot(pt,1,3)/musq)
      xl14=log(-two*Dot(pt,1,4)/musq)
      xl23=log(-two*Dot(pt,2,3)/musq)
      xl24=log(-two*Dot(pt,2,4)/musq)
      xl34=log(Q34sq/musq)

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

c-- in the q-qb piece there is no colour structure
      R1(q,q,a,0,is)=ason4pi*(
     & (xn-two/xn)*(if_mqq(z,xl13,mbar13,is)+fi_mqq(z,xl13,mbar13,is))
     &    +two/xn *(if_mqq(z,xl14,mbar14,is)+fi_mqq(z,xl14,mbar14,is))
     &    -one/xn *(ii_mqq(z,xl12,mbar12,is)+ff_mqq(z,xl34,mbar34,is)))
      R2(a,a,q,0,is)=ason4pi*(
     & (xn-two/xn)*(if_mqq(z,xl24,mbar24,is)+fi_mqq(z,xl24,mbar24,is))
     &    +two/xn *(if_mqq(z,xl23,mbar23,is)+fi_mqq(z,xl23,mbar23,is))
     &    -one/xn *(ii_mqq(z,xl12,mbar12,is)+ff_mqq(z,xl34,mbar34,is)))
      R1(a,a,q,0,is)=ason4pi*(
     & (xn-two/xn)*(if_mqq(z,xl14,mbar14,is)+fi_mqq(z,xl14,mbar14,is))
     &    +two/xn *(if_mqq(z,xl13,mbar13,is)+fi_mqq(z,xl13,mbar13,is))
     &    -one/xn *(ii_mqq(z,xl12,mbar12,is)+ff_mqq(z,xl34,mbar34,is)))
      R2(q,q,a,0,is)=ason4pi*(
     & (xn-two/xn)*(if_mqq(z,xl23,mbar23,is)+fi_mqq(z,xl23,mbar23,is))
     &    +two/xn *(if_mqq(z,xl24,mbar24,is)+fi_mqq(z,xl24,mbar24,is))
     &    -one/xn *(ii_mqq(z,xl12,mbar12,is)+ff_mqq(z,xl34,mbar34,is)))

      do icol=1,2
      R1(q,q,a,icol,is)=R1(q,q,a,0,is)
      R2(a,a,q,icol,is)=R2(a,a,q,0,is)
      R1(a,a,q,icol,is)=R1(a,a,q,0,is)
      R2(q,q,a,icol,is)=R2(q,q,a,0,is)
      enddo

c--- no colour structure for gq either
      tempqg=ason2pi*tr*ii_qg(z,xl12,is)
      tempgq=ason4pi*two*cf*ii_gq(z,xl12,is)
      do icol=0,2
      R1(q,g,q,icol,is)=tempqg
      R1(a,g,q,icol,is)=tempqg
      R1(q,g,a,icol,is)=tempqg
      R1(a,g,a,icol,is)=tempqg
      R2(q,g,q,icol,is)=tempqg
      R2(a,g,q,icol,is)=tempqg
      R2(q,g,a,icol,is)=tempqg
      R2(a,g,a,icol,is)=tempqg

      R1(g,q,g,icol,is)=tempgq
      R1(g,a,g,icol,is)=tempgq
      R2(g,q,g,icol,is)=tempgq
      R2(g,a,g,icol,is)=tempgq
      enddo

c-- in the g-g piece, this separation is required
      R1(g,g,g,0,is)=ason4pi*xn*(if_mgg(z,xl14,mbar14,is)
     &                          +if_mgg(z,xl13,mbar13,is)
     &                          +fi_mqq(z,xl14,mbar14,is)
     &                          +fi_mqq(z,xl13,mbar13,is)
     &                          -ff_mqq(z,xl34,mbar34,is))
     &              -ason4pi/xn*(ff_mqq(z,xl34,mbar34,is))
      R1(g,g,g,1,is)=ason4pi*xn*(if_mgg(z,xl14,mbar14,is)
     &                          +fi_mqq(z,xl14,mbar14,is)
     &                          +ii_mgg(z,xl12,mbar12,is))
     &              -ason4pi/xn*(ff_mqq(z,xl34,mbar34,is))
      R1(g,g,g,2,is)=ason4pi*xn*(if_mgg(z,xl13,mbar13,is)
     &                          +fi_mqq(z,xl13,mbar13,is)
     &                          +ii_mgg(z,xl12,mbar12,is))
     &              -ason4pi/xn*(ff_mqq(z,xl34,mbar34,is))
      R2(g,g,g,0,is)=ason4pi*xn*(if_mgg(z,xl24,mbar24,is)
     &                          +if_mgg(z,xl23,mbar23,is)
     &                          +fi_mqq(z,xl24,mbar24,is)
     &                          +fi_mqq(z,xl23,mbar23,is)
     &                          -ff_mqq(z,xl34,mbar34,is))
     &              -ason4pi/xn*(ff_mqq(z,xl34,mbar34,is))
      R2(g,g,g,1,is)=ason4pi*xn*(if_mgg(z,xl23,mbar23,is)
     &                          +fi_mqq(z,xl23,mbar23,is)
     &                          +ii_mgg(z,xl12,mbar12,is))
     &              -ason4pi/xn*(ff_mqq(z,xl34,mbar34,is))
      R2(g,g,g,2,is)=ason4pi*xn*(if_mgg(z,xl24,mbar24,is)
     &                          +fi_mqq(z,xl24,mbar24,is)
     &                          +ii_mgg(z,xl12,mbar12,is))
     &              -ason4pi/xn*(ff_mqq(z,xl34,mbar34,is))

      enddo

      return
      end

