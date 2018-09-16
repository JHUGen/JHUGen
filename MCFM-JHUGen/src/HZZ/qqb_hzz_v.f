      subroutine qqb_hzz_v(p,msqv)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'scale.f'
      double precision msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     . p(mxpart,4),dot,xl12
      integer j,k
      
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0d0
      enddo
      enddo

      call qqb_hzz(p,msq)

      xl12=log(two*dot(p,1,2)/musq)
      scheme='dred'
c--sum of virtual diagram in DRED and UV counterterm including 
C--term required to bring into the MSbar scheme

      msqv(0,0)=ason2pi*xn*2d0*(
     .  -epinv*(epinv2-xl12)-0.5d0*xl12**2+11d0/6d0+0.5d0*pisq
     .  -((11d0-two*dble(nf)/xn)*epinv-1d0)/6d0)*msq(0,0)

      return
      end
     
