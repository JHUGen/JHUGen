      subroutine qqb_wbjet_z(p,z)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell                                            *
*     January 2004.                                                    *
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
     &                 ii_qq,ii_qg,
     &                 if_qq,
     &                 fi_qq,
     &                 ff_qq
      real(dp):: tempqg
      integer:: is

      xl12=log(+two*dot(p,1,2)/musq)
      xl15=log(-two*dot(p,1,5)/musq)
      xl16=log(-two*dot(p,1,6)/musq)
      xl25=log(-two*dot(p,2,5)/musq)
      xl26=log(-two*dot(p,2,6)/musq)
      xl56=log(+two*dot(p,5,6)/musq)

************************************************************************
*     Contributions from QQQQ matrix elements                          *
************************************************************************
c--- Note that, in all cases, the (1) contribution is equal to the
c--- (2) contribution, after the interchange of 5 and 6
c--- In this way, the 2 pieces provide all the necessary terms to deal
c--- with the fact that particle 5 is the b-quark, regardless of
c--- whether the initial state is, e.g. (bq) or (qb)

c--- QUARK-QUARK contributions
      do is=1,3
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
      R1(a,a,a,1,is)=R1(q,q,q,1,is)
      R1(a,a,a,2,is)=R1(q,q,q,2,is)
      R2(a,a,a,1,is)=R2(q,q,q,1,is)
      R2(a,a,a,2,is)=R2(q,q,q,2,is)

      enddo

c--- QUARK-ANTIQUARK contributions
      do is=1,3
      R1(q,q,a,1,is)=ason4pi*(
     & -(if_qq(z,xl16,is)+fi_qq(z,xl16,is))/xn
     & +(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*two/xn
     & +ii_qq(z,xl12,is)*(xn-two/xn)
     & +ff_qq(z,xl56,is)*(xn-two/xn))
      R1(q,q,a,2,is)=ason4pi*(
     & -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))/xn
     & +(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*two/xn
     & +ii_qq(z,xl12,is)*(xn-two/xn)
     & +ff_qq(z,xl56,is)*(xn-two/xn))
      R2(a,a,q,1,is)=ason4pi*(
     & +(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*two/xn
     & -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))/xn
     & +ii_qq(z,xl12,is)*(xn-two/xn)
     & +ff_qq(z,xl56,is)*(xn-two/xn))
      R2(a,a,q,2,is)=ason4pi*(
     & +(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*two/xn
     & -(if_qq(z,xl26,is)+fi_qq(z,xl26,is))/xn
     & +ii_qq(z,xl12,is)*(xn-two/xn)
     & +ff_qq(z,xl56,is)*(xn-two/xn))

      enddo

c--- ANTIQUARK-QUARK contributions
c      do is=1,3
c      R1(a,a,q,1,is)=ason4pi*(
c     & +(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*two/xn
c     & -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))/xn
c     & +ii_qq(z,xl12,is)*(xn-two/xn)
c     & +ff_qq(z,xl56,is)*(xn-two/xn))
c      R1(a,a,q,2,is)=ason4pi*(
c     & +(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*two/xn
c     & -(if_qq(z,xl16,is)+fi_qq(z,xl16,is))/xn
c     & +ii_qq(z,xl12,is)*(xn-two/xn)
c     & +ff_qq(z,xl56,is)*(xn-two/xn))
c      R2(q,q,a,1,is)=ason4pi*(
c     & -(if_qq(z,xl26,is)+fi_qq(z,xl26,is))/xn
c     & +(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*two/xn
c     & +ii_qq(z,xl12,is)*(xn-two/xn)
c     & +ff_qq(z,xl56,is)*(xn-two/xn))
c      R2(q,q,a,2,is)=ason4pi*(
c     & -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))/xn
c     & +(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*two/xn
c     & +ii_qq(z,xl12,is)*(xn-two/xn)
c     & +ff_qq(z,xl56,is)*(xn-two/xn))
c      enddo

C---RKE modification
c--- ANTIQUARK-QUARK contributions
      do is=1,3
      R1(a,a,q,1,is)=R1(q,q,a,1,is)
      R1(a,a,q,2,is)=R1(q,q,a,2,is)
      R2(q,q,a,1,is)=R2(a,a,q,1,is)
      R2(q,q,a,2,is)=R2(a,a,q,2,is)
      enddo


      do is=1,3
c--- GLUON-QUARK contributions
      tempqg=ason4pi*ii_qg(z,xl12,is)
      R1(q,g,q,1,is)=tempqg
      R1(q,g,q,2,is)=tempqg

      R1(a,g,q,1,is)=tempqg
      R1(a,g,q,2,is)=tempqg

c--- GLUON-ANTIQUARK contributions
      R1(q,g,a,1,is)=tempqg
      R1(q,g,a,2,is)=tempqg

      R1(a,g,a,1,is)=tempqg
      R1(a,g,a,2,is)=tempqg


c--- QUARK-GLUON contributions
      R2(q,g,q,1,is)=tempqg
      R2(q,g,q,2,is)=tempqg

      R2(a,g,q,1,is)=tempqg
      R2(a,g,q,2,is)=tempqg

c--- ANTIQUARK-GLUON contributions
      R2(q,g,a,1,is)=tempqg
      R2(q,g,a,2,is)=tempqg

      R2(a,g,a,1,is)=tempqg
      R2(a,g,a,2,is)=tempqg

      enddo

      return
      end

