      subroutine tree(mq,ig,is,ie,in,it,amp)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_com.f'
      integer:: is,ig,ie,in,it,i,j
      complex(dp):: amp(2,2)
      real(dp):: mq,prop

      prop=sqrt((real(za(ie,in)*zb(in,ie))-wmass**2)**2
     &          +(wmass*wwidth)**2)

c---- helitities: amp(ht,hg)
c--- labels on amplitudes represent helicities for the heavy quark
c--- and the gluon respectively
c---   1 = negative helitity, 2 = positive helitity
c--- heavy quark momentum is made massless (it) with the gluon momentum ig
      amp(1,2)=za(ie,it)
     & /za(ig,is)/za(ig,it)*(za(is,it)*zb(is,in)+za(ig,it)*zb(ig,in))
      amp(2,2)=-mq/za(ig,it)*za(ig,ie)
     & /za(ig,is)/za(ig,it)*(za(is,it)*zb(is,in)+za(ig,it)*zb(ig,in))
      amp(1,1)=-(za(ig,ie)*zb(ig,is)+za(ie,it)*zb(is,it))*zb(is,in)
     & /zb(ig,is)/zb(ig,it)
      amp(2,1)=-mq/za(ig,it)/zb(it,ig)*za(ig,ie)*zb(is,in)*zb(is,it)
     & /zb(ig,is)

      do i=1,2
      do j=1,2
      amp(i,j)=amp(i,j)/prop
      enddo
      enddo

      return
      end
