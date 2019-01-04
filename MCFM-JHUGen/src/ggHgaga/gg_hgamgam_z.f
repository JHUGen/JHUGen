      subroutine gg_hgamgam_z(p,z)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      integer is
      double precision z,p(mxpart,4),xl12,dot,ii_gg,ii_gq,tempgg,tempgq

      xl12=log(two*dot(p,1,2)/musq)

      do is=1,3
      tempgg=ason4pi*2d0*xn*ii_gg(z,xl12,is)
      tempgq=ason4pi*2d0*cf*ii_gq(z,xl12,is)
      Q1(g,g,g,is)=tempgg
      Q2(g,g,g,is)=tempgg

      Q1(g,q,g,is)=tempgq
      Q1(g,a,g,is)=tempgq

      Q2(g,q,g,is)=tempgq
      Q2(g,a,g,is)=tempgq
      enddo
      return
      end


