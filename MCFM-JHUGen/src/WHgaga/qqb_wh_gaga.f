      subroutine qqb_wh_gaga(p,msq)
c---Matrix element squared averaged over initial colors and spins
c for nwz=1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> ga(p5)+ga(p6)
c for nwz=-1
c     q(-p1)+qbar(-p2) -->  H  + W
c                           |    |
c                           |    --> e^-(p3)+nubar(p4)
c                           |
c                           ---> ga(p5)+ga(p6)
c---- Extension to photon decay contributed by Fabian Stoeckli
      implicit none 
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ckm.f'
      integer j,k
      double precision p(mxpart,4)
      double precision s,prop,fac,qqbWH,qbqWH,s56
      double precision msq(-nf:nf,-nf:nf),hdecay      
      double precision msqgamgam

      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      s56=s(5,6)+2*mb**2
c---cut to ensure hard process
      if (
     .      (s(5,6) .lt. four*mbsq) 
     . .or. (s(1,5)*s(2,5)/s(1,2) .lt. mbsq) 
     . .or. (s(1,6)*s(2,6)/s(1,2) .lt. mbsq) ) return

c---calculate the 2 W propagators
      prop=     ((s(1,2)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
      
      fac=xn*gwsq**3*wmass**2/prop
      hdecay=msqgamgam(hmass)/((s56-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay

c-- Old form of this matrix element (modified to facilitate extension
c--- to H->WW decay)
c      spinave only (color factors cancel)
c       fac=spinave*gw**8*0.5d0*mbsq*(s56-4d0*mb**2)/prop
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

