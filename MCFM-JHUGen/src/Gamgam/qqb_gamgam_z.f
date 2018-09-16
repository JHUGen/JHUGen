      subroutine qqb_gamgam_z(p,z)
      implicit none
      include 'types.f'
************************************************************************
*     Authors: R.K. Ellis and John M. Campbell                         *
*     December, 2010.                                                  *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'qcdcouple.f'
      include 'facscale.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      integer:: is
      real(dp):: z,xl12,p(mxpart,4),dot,ii_qq,ii_qg,
     & tempqq,tempqg,tempgg

      xl12=log(two*dot(p,1,2)/musq)
c----contributions for one leg

      tempgg=0._dp
c---- NOTE: this contribution removes the AP subtraction for the qa->g
c----       splitting that is normally automatically included in virtint
c      do is=1,3
c      Q1(g,q,g,is)=0._dp
c      Q1(g,a,g,is)=0._dp
c      Q2(g,q,g,is)=0._dp
c      Q2(g,a,g,is)=0._dp
c      enddo
c      Q1(g,q,g,2)=-ason2pi*Cf*(1._dp+(1._dp-z)**2)/z
c     &            *(epinv+log(scale/facscale))
c      Q1(g,a,g,2)=Q1(g,q,g,2)
c      Q2(g,q,g,2)=Q1(g,q,g,2)
c      Q2(g,a,g,2)=Q1(g,q,g,2)
      
      do is=1,3
      tempqq=+ason2pi*cf*ii_qq(z,xl12,is)
      tempqg=+ason2pi*tr*ii_qg(z,xl12,is)
! JC removing gg contribution
c--- Comment out the following line to remove gg contribution
!      tempgg=+ason2pi*xn*ii_gg(z,xl12,is)
      
      Q1(q,q,a,is)=tempqq
      Q2(a,a,q,is)=tempqq
      Q1(a,a,q,is)=tempqq
      Q2(q,q,a,is)=tempqq

      Q2(a,g,q,is)=tempqg
      Q2(q,g,a,is)=tempqg
      Q1(a,g,q,is)=tempqg
      Q1(q,g,a,is)=tempqg

      Q1(g,g,g,is)=tempgg
      Q2(g,g,g,is)=tempgg
      enddo

      return
      end
