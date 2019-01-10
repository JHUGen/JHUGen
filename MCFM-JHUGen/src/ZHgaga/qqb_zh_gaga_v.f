      subroutine qqb_zh_gaga_v(p,msqv)
      implicit none
      integer j,k
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      double precision p(mxpart,4),
     . msqv(-nf:nf,-nf:nf),msq(-nf:nf,-nf:nf),
     . dot,virt,xl12

      xl12=log(two*dot(p,1,2)/musq)
      scheme='dred'
c---  calculate lowest order matrix element
      call qqb_zh_gaga(p,msq)
c---calculate the multiple of the lowest order
      virt=ason2pi*cf*(-2d0*(epinv*epinv2-epinv*xl12+half*xl12**2)
     .                 -3d0*(epinv-xl12)
     .                 +pisq-7d0)
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=virt*msq(j,k)
      enddo
      enddo
      end
