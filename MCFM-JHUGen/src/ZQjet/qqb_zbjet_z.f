      subroutine qqb_zbjet_z(p,z)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell                                            *
*     August, 2005.                                                    *
************************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'agq.f'
      include 'PR_cs_new.f'
      real(dp):: z,p(mxpart,4),dot
      real(dp):: xl12,xl15,xl16,xl25,xl26,xl56
      real(dp)::
     &                 ii_qq,ii_qg,ii_gq,ii_gg,
     &                 if_qq,if_gg,
     &                 fi_qq,fi_gg,
     &                 ff_qq,ff_gg
      real(dp):: tempgq,tempqg
      integer:: is

      xl12=log(+two*dot(p,1,2)/musq)
      xl15=log(-two*dot(p,1,5)/musq)
      xl16=log(-two*dot(p,1,6)/musq)
      xl25=log(-two*dot(p,2,5)/musq)
      xl26=log(-two*dot(p,2,6)/musq)
      xl56=log(+two*dot(p,5,6)/musq)

************************************************************************
*     Contributions from QQGG matrix elements                          *
************************************************************************

c--- QUARK-GLUON contributions
c--- additional final-final pieces that are 0:
c--- ason4pi*xn*(ff_qq(z,xl56,2) + ff_gg(z,xl56,2)/2._dp) * (1)
      do is=1,3
      R1(q,q,g,1,is)=ason4pi*half*xn*(ii_qq(z,xl12,is)+ii_qq(z,xl12,is)
     &              +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2._dp)
      R1(q,q,g,2,is)=ason4pi*xn*(if_qq(z,xl16,is)+fi_gg(z,xl16,is)/2._dp)
      R2(g,g,q,1,is)=ason4pi*half*xn
     & *(2._dp*if_gg(z,xl26,is)+fi_gg(z,xl26,is)
     &      +ii_gg(z,xl12,is)+ii_gg(z,xl12,is)
     &      +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2._dp)
      R2(g,g,q,2,is)=ason4pi*xn*(if_gg(z,xl26,is)+fi_gg(z,xl26,is)/2._dp
     &                          +if_gg(z,xl25,is)+fi_qq(z,xl25,is))
      R2(a,g,q,1,is)=ason4pi*xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      R2(a,g,q,2,is)=ason4pi*xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      enddo

      do is=1,3
      R1(q,q,g,1,is)=R1(q,q,g,1,is)
     &           -ason4pi/xn*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
      R1(q,q,g,2,is)=R1(q,q,g,2,is)
     &           -ason4pi/xn*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
      R1(q,q,g,0,is)=R1(q,q,g,0,is)+ason4pi*half*xn*
     &            (2._dp*if_qq(z,xl16,is)+fi_gg(z,xl16,is)
     &            +ii_qq(z,xl12,is)+ii_qq(z,xl12,is)
     &            +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2._dp)
      R2(g,g,q,0,is)=R2(g,g,q,0,is)+ason4pi*half*xn*
     &            (2._dp*if_gg(z,xl25,is)+2._dp*fi_qq(z,xl25,is)
     &            +ii_gg(z,xl12,is)+ii_gg(z,xl12,is)
     &            +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2._dp)
      R2(a,g,q,1,is)=R2(a,g,q,1,is)-
     & ason4pi/xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      R2(a,g,q,2,is)=R2(a,g,q,2,is)-
     & ason4pi/xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      R2(a,g,q,0,is)=R2(a,g,q,0,is)+
     & ason4pi*2._dp*xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      enddo

      do is=1,3
      R1(q,q,g,0,is)=R1(q,q,g,0,is)-ason4pi*(xn+1._dp/xn)
     &                       *(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
      R2(a,g,q,0,is)=R2(a,g,q,0,is)-
     & ason4pi*(xn+1._dp/xn)*(aveqg/aveqq)*ii_qg(z,xl12,is)
      enddo

      do is=1,3
      R2(q,g,q,0,is)=R2(a,g,q,0,is)
      R2(q,g,q,1,is)=R2(a,g,q,1,is)
      R2(q,g,q,2,is)=R2(a,g,q,2,is)
      enddo

c--- GLUON-QUARK contributions
c--- additional final-final pieces that are 0:
c--- ason4pi*xn*(ff_qq(z,xl56,2) + ff_gg(z,xl56,2)/2._dp) * (1)
      do is=1,3
      R1(g,g,q,1,is)=ason4pi*half*xn*(ii_gg(z,xl12,is)+ii_gg(z,xl12,is)
     &                           +2._dp*if_gg(z,xl16,is)+fi_gg(z,xl16,is)
     &              +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2._dp)
      R1(g,g,q,2,is)=ason4pi*xn*(if_gg(z,xl16,is)+fi_gg(z,xl16,is)/2._dp
     &                           +if_gg(z,xl15,is)+fi_qq(z,xl15,is))
      R2(q,q,g,1,is)=ason4pi*half*xn*(ii_qq(z,xl12,is)+ii_qq(z,xl12,is)
     &              +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2._dp)
      R2(q,q,g,2,is)=ason4pi*xn*(if_qq(z,xl26,is)+fi_gg(z,xl26,is)/2._dp)
      R1(a,g,q,1,is)=ason4pi*xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      R1(a,g,q,2,is)=ason4pi*xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      enddo

      do is=1,3
      R2(q,q,g,1,is)=R2(q,q,g,1,is)
     &           -ason4pi/xn*(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
      R2(q,q,g,2,is)=R2(q,q,g,2,is)
     &           -ason4pi/xn*(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
      R1(g,g,q,0,is)=R1(g,g,q,0,is)+ason4pi*half*xn*
     &            (2._dp*if_gg(z,xl15,is)+2._dp*fi_qq(z,xl15,is)
     &            +ii_gg(z,xl12,is)+ii_gg(z,xl12,is)
     &            +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2._dp)
      R2(q,q,g,0,is)=R2(q,q,g,0,is)+ason4pi*half*xn*
     &            (2._dp*if_qq(z,xl26,is)+fi_gg(z,xl26,is)
     &            +ii_qq(z,xl12,is)+ii_qq(z,xl12,is)
     &            +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2._dp)
      R1(a,g,q,1,is)=R1(a,g,q,1,is)-
     & ason4pi/xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      R1(a,g,q,2,is)=R1(a,g,q,2,is)-
     & ason4pi/xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      R1(a,g,q,0,is)=R1(a,g,q,0,is)+
     & ason4pi*2._dp*xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      enddo

      do is=1,3
      R2(q,q,g,0,is)=R2(q,q,g,0,is)-ason4pi*(xn+1._dp/xn)
     &                       *(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
      R1(a,g,q,0,is)=R1(a,g,q,0,is)-
     & ason4pi*(xn+1._dp/xn)*(aveqg/aveqq)*ii_qg(z,xl12,is)
      enddo

      do is=1,3
      R1(q,g,q,0,is)=R1(a,g,q,0,is)
      R1(q,g,q,1,is)=R1(a,g,q,1,is)
      R1(q,g,q,2,is)=R1(a,g,q,2,is)
      enddo

c--- GLUON-ANTIQUARK contributions
c--- additional final-final pieces that are 0:
c--- ason4pi*xn*(ff_qq(z,xl56,2) + ff_gg(z,xl56,2)/2._dp) * (1)
      do is=1,3
      R1(g,g,a,1,is)= R1(g,g,q,2,is)
      R1(g,g,a,2,is)= R1(g,g,q,1,is)
      R1(g,g,a,0,is)= R1(g,g,q,0,is)
      R2(a,a,g,1,is)= R2(q,q,g,2,is)
      R2(a,a,g,2,is)= R2(q,q,g,1,is)
      R2(a,a,g,0,is)= R2(q,q,g,0,is)
      R1(q,g,a,0,is)= R1(a,g,q,0,is)
      R1(q,g,a,1,is)= R1(a,g,q,1,is)
      R1(q,g,a,2,is)= R1(a,g,q,2,is)
      R1(a,g,a,0,is)= R1(q,g,q,0,is)
      R1(a,g,a,1,is)= R1(q,g,q,1,is)
      R1(a,g,a,2,is)= R1(q,g,q,2,is)
      enddo

c--- ANTIQUARK-GLUON contributions
c--- additional final-final pieces that are 0:
c--- ason4pi*xn*(ff_qq(z,xl56,2) + ff_gg(z,xl56,2)/2._dp) * (1)
      do is=1,3
      R1(a,a,g,1,is)= R1(q,q,g,2,is)
      R1(a,a,g,2,is)= R1(q,q,g,1,is)
      R1(a,a,g,0,is)= R1(q,q,g,0,is)
      R2(g,g,a,1,is)= R2(g,g,q,2,is)
      R2(g,g,a,2,is)= R2(g,g,q,1,is)
      R2(g,g,a,0,is)= R2(g,g,q,0,is)
      R2(q,g,a,0,is)= R2(a,g,q,0,is)
      R2(q,g,a,1,is)= R2(a,g,q,1,is)
      R2(q,g,a,2,is)= R2(a,g,q,2,is)
      R2(a,g,a,0,is)= R2(q,g,q,0,is)
      R2(a,g,a,1,is)= R2(q,g,q,1,is)
      R2(a,g,a,2,is)= R2(q,g,q,2,is)
      enddo

c--- GLUON-GLUON (isub=2) contribution
c--- originating from diagrams such as:
c---
c---                   Q
c---                 /
c---                /
c---               /
c---         oooooo------ooooooooo
c---                    |
c---                    |
c---                    |
c---         ooooooooooo---------- Qbar
c---
      do is=1,3
      R1(q,g,g,0,is)=ason4pi*ii_qg(z,xl12,is)
      R1(q,g,g,1,is)=R1(q,g,g,0,is)
      R1(q,g,g,2,is)=R1(q,g,g,0,is)
      R2(q,g,g,0,is)=R1(q,g,g,0,is)
      R2(q,g,g,1,is)=R1(q,g,g,0,is)
      R2(q,g,g,2,is)=R1(q,g,g,0,is)
      enddo

C--- off-diagonal Qflag->Gflag pieces, for example:
c---
c---                   q
c---                 /
c---                /
c---               /
c---       q ------ooooooooooooooo
c---                    o
c---                    o
c---                    o
c---       Q --------------------- Q
c---
      do is=1,3
      tempgq=ason4pi*two*cf*ii_gq(z,xl12,is)
c--- Q-Qbar diagrams
      R1(g,q,a,0,is)=tempgq
      R1(g,q,a,1,is)=tempgq
      R1(g,q,a,2,is)=tempgq
      R2(g,a,q,0,is)=tempgq
      R2(g,a,q,1,is)=tempgq
      R2(g,a,q,2,is)=tempgq

c--- Q-Q diagrams
      R1(g,q,q,0,is)=tempgq
      R1(g,q,q,1,is)=tempgq
      R1(g,q,q,2,is)=tempgq
      R2(g,q,q,0,is)=tempgq
      R2(g,q,q,1,is)=tempgq
      R2(g,q,q,2,is)=tempgq

c--- Qbar-Qbar diagrams
      R1(g,a,a,0,is)=tempgq
      R1(g,a,a,1,is)=tempgq
      R1(g,a,a,2,is)=tempgq
      R2(g,a,a,0,is)=tempgq
      R2(g,a,a,1,is)=tempgq
      R2(g,a,a,2,is)=tempgq

c--- Qbar-Q diagrams
      R1(g,a,q,0,is)=tempgq
      R1(g,a,q,1,is)=tempgq
      R1(g,a,q,2,is)=tempgq
      R2(g,q,a,0,is)=tempgq
      R2(g,q,a,1,is)=tempgq
      R2(g,q,a,2,is)=tempgq
      enddo

************************************************************************
*     Contributions from QQQQ matrix elements                          *
************************************************************************
c--- QUARK-QUARK contributions
      do is=1,3
      R1(q,q,q,0,is)=ason4pi*(
     & -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))/xn
     & -(if_qq(z,xl16,is)+fi_qq(z,xl16,is))/xn
     & +ii_qq(z,xl12,is)*(xn+one/xn)
     & +ff_qq(z,xl56,is)*(xn+one/xn))
      R1(q,q,q,1,is)=ason4pi*(
     & -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))/xn
     & +(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*(xn-two/xn)
     & +ii_qq(z,xl12,is)*two/xn
     & +ff_qq(z,xl56,is)*two/xn)
      R1(q,q,q,2,is)=ason4pi*(
     & +(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*(xn-two/xn)
     & -(if_qq(z,xl16,is)+fi_qq(z,xl16,is))/xn
     & +ii_qq(z,xl12,is)*two/xn
     & +ff_qq(z,xl56,is)*two/xn)
      R2(q,q,q,0,is)=ason4pi*(
     & -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))/xn
     & -(if_qq(z,xl26,is)+fi_qq(z,xl26,is))/xn
     & +ii_qq(z,xl12,is)*(xn+one/xn)
     & +ff_qq(z,xl56,is)*(xn+one/xn))
      R2(q,q,q,1,is)=ason4pi*(
     & +(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*(xn-two/xn)
     & -(if_qq(z,xl26,is)+fi_qq(z,xl26,is))/xn
     & +ii_qq(z,xl12,is)*two/xn
     & +ff_qq(z,xl56,is)*two/xn)
      R2(q,q,q,2,is)=ason4pi*(
     & -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))/xn
     & +(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*(xn-two/xn)
     & +ii_qq(z,xl12,is)*two/xn
     & +ff_qq(z,xl56,is)*two/xn)

      enddo


c--- ANTIQUARK-ANTIQUARK contributions
      do is=1,3
      R1(a,a,a,0,is)=R1(q,q,q,0,is)
      R1(a,a,a,1,is)=R1(q,q,q,1,is)
      R1(a,a,a,2,is)=R1(q,q,q,2,is)
      R2(a,a,a,0,is)=R2(q,q,q,0,is)
      R2(a,a,a,1,is)=R2(q,q,q,1,is)
      R2(a,a,a,2,is)=R2(q,q,q,2,is)

      enddo

c--- QUARK-ANTIQUARK contributions
      do is=1,3
      R1(q,q,a,0,is)=ason4pi*(
     & -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))/xn
     & +(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*(xn+one/xn)
     & -ii_qq(z,xl12,is)/xn
     & -ff_qq(z,xl56,is)/xn)
      R1(q,q,a,1,is)=ason4pi*(
     & +(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*(xn-two/xn)
     & +(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*two/xn
     & -ii_qq(z,xl12,is)/xn
     & -ff_qq(z,xl56,is)/xn)
      R1(q,q,a,2,is)=ason4pi*(
     & -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))/xn
     & +(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*two/xn
     & +ii_qq(z,xl12,is)*(xn-two/xn)
     & +ff_qq(z,xl56,is)*(xn-two/xn))
      R2(a,a,q,0,is)=ason4pi*(
     & +(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*(xn+one/xn)
     & -(if_qq(z,xl26,is)+fi_qq(z,xl26,is))/xn
     & -ii_qq(z,xl12,is)/xn
     & -ff_qq(z,xl56,is)/xn)
      R2(a,a,q,1,is)=ason4pi*(
     & +(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*two/xn
     & +(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*(xn-two/xn)
     & -ii_qq(z,xl12,is)/xn
     & -ff_qq(z,xl56,is)/xn)
      R2(a,a,q,2,is)=ason4pi*(
     & +(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*two/xn
     & -(if_qq(z,xl26,is)+fi_qq(z,xl26,is))/xn
     & +ii_qq(z,xl12,is)*(xn-two/xn)
     & +ff_qq(z,xl56,is)*(xn-two/xn))

      enddo

c--- ANTIQUARK-QUARK contributions
      do is=1,3
      R1(a,a,q,0,is)=ason4pi*(
     & +(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*(xn+one/xn)
     & -(if_qq(z,xl16,is)+fi_qq(z,xl16,is))/xn
     & -ii_qq(z,xl12,is)/xn
     & -ff_qq(z,xl56,is)/xn)
      R1(a,a,q,1,is)=ason4pi*(
     & +(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*two/xn
     & +(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*(xn-two/xn)
     & -ii_qq(z,xl12,is)/xn
     & -ff_qq(z,xl56,is)/xn)
      R1(a,a,q,2,is)=ason4pi*(
     & +(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*two/xn
     & -(if_qq(z,xl16,is)+fi_qq(z,xl16,is))/xn
     & +ii_qq(z,xl12,is)*(xn-two/xn)
     & +ff_qq(z,xl56,is)*(xn-two/xn))
      R2(q,q,a,0,is)=ason4pi*(
     & -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))/xn
     & +(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*(xn+one/xn)
     & -ii_qq(z,xl12,is)/xn
     & -ff_qq(z,xl56,is)/xn)
      R2(q,q,a,1,is)=ason4pi*(
     & +(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*(xn-two/xn)
     & +(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*two/xn
     & -ii_qq(z,xl12,is)/xn
     & -ff_qq(z,xl56,is)/xn)
      R2(q,q,a,2,is)=ason4pi*(
     & -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))/xn
     & +(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*two/xn
     & +ii_qq(z,xl12,is)*(xn-two/xn)
     & +ff_qq(z,xl56,is)*(xn-two/xn))

      enddo

      do is=1,3

c--- These contributions arise from, e.g.
c---
c---                   q
c---                 /
c---                /
c---               /
c---         oooooo--------------- qbar
c---                    o
c---                    o
c---                    o
c---       Q --------------------- Q
c---

c--- GLUON-QUARK contributions
      tempqg=ason4pi*ii_qg(z,xl12,is)
      R1(q,g,q,0,is)=tempqg
      R1(q,g,q,1,is)=tempqg
      R1(q,g,q,2,is)=tempqg

      R1(a,g,q,0,is)=tempqg
      R1(a,g,q,1,is)=tempqg
      R1(a,g,q,2,is)=tempqg

c--- GLUON-ANTIQUARK contributions
      R1(q,g,a,0,is)=tempqg
      R1(q,g,a,1,is)=tempqg
      R1(q,g,a,2,is)=tempqg

      R1(a,g,a,0,is)=tempqg
      R1(a,g,a,1,is)=tempqg
      R1(a,g,a,2,is)=tempqg


c--- QUARK-GLUON contributions
      R2(q,g,q,0,is)=tempqg
      R2(q,g,q,1,is)=tempqg
      R2(q,g,q,2,is)=tempqg

      R2(a,g,q,0,is)=tempqg
      R2(a,g,q,1,is)=tempqg
      R2(a,g,q,2,is)=tempqg

c--- ANTIQUARK-GLUON contributions
      R2(q,g,a,0,is)=tempqg
      R2(q,g,a,1,is)=tempqg
      R2(q,g,a,2,is)=tempqg

      R2(a,g,a,0,is)=tempqg
      R2(a,g,a,1,is)=tempqg
      R2(a,g,a,2,is)=tempqg

      enddo

      return
      end

