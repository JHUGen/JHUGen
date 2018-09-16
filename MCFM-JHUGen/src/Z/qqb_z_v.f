      subroutine qqb_z_v(p,msqv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      integer:: j,k
      real(dp):: p(mxpart,4),
     & msqv(-nf:nf,-nf:nf),msq(-nf:nf,-nf:nf),
     & dot,virt,xl12

      scheme='dred'
      xl12=log(two*dot(p,1,2)/musq)

c---  calculate lowest order matrix element
      call qqb_z(p,msq)
c---calculate the multiple of the lowest order
      virt=ason2pi*cf*(-2._dp*(epinv*epinv2-epinv*xl12+half*xl12**2)
     &                 -3._dp*(epinv-xl12)
     &                 +pisq-7._dp)
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=virt*msq(j,k)
      enddo
      enddo
      end

