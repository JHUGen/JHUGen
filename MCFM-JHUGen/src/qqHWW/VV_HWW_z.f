      subroutine VV_HWW_z(p,z)
      implicit none
      include 'types.f'
*     Weak Boson Fusion by V-V exchange                                *

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      integer:: is
      real(dp):: z,xl12,xl17,xl28,
     & tempqq1,tempqq2,tempgq,tempqg,
     & p(mxpart,4),dot,if_qq,fi_qq,ii_qg

      xl12=log(+two*dot(p,1,2)/musq)
      xl17=log(-two*dot(p,1,7)/musq)
      xl28=log(-two*dot(p,2,8)/musq)

      do is=1,3
      tempqq1=+ason2pi*cf*(if_qq(z,xl17,is)+fi_qq(z,xl17,is))
      tempqq2=+ason2pi*cf*(if_qq(z,xl28,is)+fi_qq(z,xl28,is))
      tempgq =+ason2pi*tr*ii_qg(z,xl12,is)
      tempqg =+ason2pi*tr*ii_qg(z,xl12,is)

      Q1(q,q,q,is)=tempqq1
      Q2(q,q,q,is)=tempqq2
      Q1(a,a,a,is)=tempqq1
      Q2(a,a,a,is)=tempqq2

      Q1(q,q,a,is)=tempqq1
      Q2(q,q,a,is)=tempqq2
      Q1(a,a,q,is)=tempqq1
      Q2(a,a,q,is)=tempqq2

      Q1(q,g,q,is)=tempgq
      Q1(a,g,a,is)=tempgq
      Q2(q,g,q,is)=tempqg
      Q2(a,g,a,is)=tempqg

      Q1(q,g,a,is)=tempgq
      Q1(a,g,q,is)=tempgq
      Q2(q,g,a,is)=tempqg
      Q2(a,g,q,is)=tempqg

      enddo

      return
      end
