      subroutine qqb_gamgam_z(p,z)
************************************************************************
*     Authors: R.K. Ellis and John M. Campbell                         *
*     December, 2010.                                                  *
************************************************************************
      implicit none
      include 'constants.f'
      include 'epinv.f'
      include 'qcdcouple.f'
      include 'facscale.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      integer is
      double precision z,xl12,p(mxpart,4),dot,ii_qq,ii_qg,ii_gg,
     & tempqq,tempqg,tempgg

      xl12=dlog(two*dot(p,1,2)/musq)
c----contributions for one leg

      tempgg=0d0
c---- NOTE: this contribution removes the AP subtraction for the qa->g
c----       splitting that is normally automatically included in virtint
      do is=1,3
      Q1(g,q,g,is)=0d0
      Q1(g,a,g,is)=0d0
      Q2(g,q,g,is)=0d0
      Q2(g,a,g,is)=0d0
      enddo
      Q1(g,q,g,2)=-ason2pi*Cf*(1d0+(1d0-z)**2)/z
     &            *(epinv+dlog(scale/facscale))
      Q1(g,a,g,2)=Q1(g,q,g,2)
      Q2(g,q,g,2)=Q1(g,q,g,2)
      Q2(g,a,g,2)=Q1(g,q,g,2)
      
      do is=1,3
      tempqq=+ason2pi*cf*ii_qq(z,xl12,is)
      tempqg=+ason2pi*tr*ii_qg(z,xl12,is)
c--- Comment out the following line to remove gg contribution
      tempgg=+ason2pi*xn*ii_gg(z,xl12,is)
      
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
