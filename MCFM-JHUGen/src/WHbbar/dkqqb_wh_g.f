      subroutine dkqqb_wh_g(p,msq)
c---Matrix element squared averaged over initial colors and spins
c for nwz=1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+b(p6)+g(p7)
c for nwz=-1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> e^-(p3)+nubar(p4)
c                           |
c                           ---> b(p5)+b(p6)+g(p7)
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ckm.f'
      integer j,k
      double precision p(mxpart,4)
      double precision s,prop,fac,qqbWH,qbqWH,s567
      double precision msq(-nf:nf,-nf:nf),hdecay

      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

      msq(:,:)=0d0
      s567=s(5,6)+s(5,7)+s(6,7)+2d0*mb**2

c---calculate the 2 W propagators
      prop=     ((s(1,2)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
      
      fac=xn*gwsq**3*wmass**2/prop
      call hbbdecay_g(p,5,6,7,hdecay)
      hdecay=hdecay/((s567-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay

      qqbWH=aveqq*fac*s(1,4)*s(2,3)
      qbqWH=aveqq*fac*s(2,4)*s(1,3)

      do j=-nf,nf
      do k=-nf,nf
      if ((j .gt. 0) .and. (k .lt. 0)) msq(j,k)=Vsq(j,k)*qqbWH
      if ((j .lt. 0) .and. (k .gt. 0)) msq(j,k)=Vsq(j,k)*qbqWH
      enddo
      enddo
      
      return
      end

