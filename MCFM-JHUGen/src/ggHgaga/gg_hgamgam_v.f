      subroutine gg_hgamgam_v(p,msqv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'scale.f'
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     & p(mxpart,4),dot,xl12
      integer:: j,k
      
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      enddo
      enddo

      call gg_hgamgam(p,msq)
      xl12=log(two*dot(p,1,2)/musq)

c--sum of virtual diagram in DRED and UV counterterm including 
C--term required to bring into the MSbar scheme
      scheme='dred'

      msqv(0,0)=ason2pi*xn*2._dp*(
     &  -epinv*(epinv2-xl12)-0.5_dp*xl12**2+11._dp/6._dp+0.5_dp*pisq
     &  -((11._dp-two*real(nf,dp)/xn)*epinv-1._dp)/6._dp)*msq(0,0)

      return
      end
     
