      subroutine qqb_wh_ww_g(P,msq)
c---Matrix element squared averaged over initial colors and spins
c---for nwz=1
c     q(-p1)+qbar(-p2) -->  H  + W +g(p7)
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> W^+(nu(p5),e^+(p6))W^-(e^-(p7),nub(p8))
c---for nwz=-1
c     q(-p1)+qbar(-p2) -->  H  + W +g(p7)
c                           |    |
c                           |    --> e^-(p3)+nubar(p4)
c                           |
c                           ---> W^+(nu(p5),e^+(p6))W^-(e^-(p7),nub(p8))
c   for the moment --- radiation only from initial line
      implicit none 
      include 'constants.f'
      include 'ckm.f'
      include 'sprods_com.f'
      integer j,k
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf)
      double precision radi_ww,hdecay
      double precision qqbWHg,qbqWHg,qgWHq,gqWHq,gqbWHqb,qbgWHqb

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call dotem(9,p,s)
      call hwwdecay(p,5,6,7,8,hdecay)
      qqbWHg=aveqq*radi_ww(1,2,9,5,6,7,8,3,4)*hdecay
      qbqWHg=aveqq*radi_ww(2,1,9,5,6,7,8,3,4)*hdecay
      qgWHq=-radi_ww(1,9,2,5,6,7,8,3,4)*aveqg*hdecay
      gqWHq=-radi_ww(2,9,1,5,6,7,8,3,4)*aveqg*hdecay

      gqbWHqb=-radi_ww(9,2,1,5,6,7,8,3,4)*aveqg*hdecay
      qbgWHqb=-radi_ww(9,1,2,5,6,7,8,3,4)*aveqg*hdecay

c      write(6,*) 'qqbWHg',qqbWHg
c      write(6,*) 'qbqWHg',qbqWHg
c      write(6,*) 'qbgWHqb',qbgWHqb
c      write(6,*) 'gqbWHqb',gqbWHqb
c      write(6,*) 'qgWHq',qgWHq
c      write(6,*) 'gqWHq',gqWHq


      do j=-nf,nf
      do k=-nf,nf

      if     ((j .gt. 0) .and. (k .lt. 0)) then
          msq(j,k)=Vsq(j,k)*qqbWHg
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=Vsq(j,k)*qbqWHg
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &   (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWHq
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWHqb
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWHq
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWHqb
      endif

      enddo
      enddo
      return
      
      end


      double precision function radi_ww(j1,j2,j3,j4,j5,j6,j7,j8,j9)
      implicit none 
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7,j8,j9
      double precision s12,s13,s23,s123,s4567,prop
      double precision fac

      s12=s(j1,j2)
      s13=s(j1,j3)
      s23=s(j2,j3)
      s123=s12+s13+s23
      s4567=s(j4,j5)+s(j4,j6)+s(j4,j7)+s(j5,j6)+s(j5,j7)+s(j6,j7)
c---calculate the 2 W propagators
      prop=       ((s123-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(j8,j9)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s4567-hmass**2)**2+(hmass*hwidth)**2)
      
      fac=2d0*cf*xn*gsq*gwsq**3*wmass**2/prop

C---old
c      radi_ww=s12/s13/s23
c     & *(2d0*s(j1,j9)*s(j2,j8)+s(j1,j9)*s(j3,j8)+s(j2,j8)*s(j3,j9))
c     & +(s(j1,j9)*s(j2,j8)+s(j2,j8)*s(j3,j9)-s(j1,j8)*s(j1,j9))/s13
c     & +(s(j1,j9)*s(j2,j8)+s(j1,j9)*s(j3,j8)-s(j2,j8)*s(j2,j9))/s23
c---
      radi_ww=
     & (s(j2,j8)*((s12+s23)*(s(j1,j9)+s(j3,j9))-s13*s(j2,j9))
     & +s(j1,j9)*((s12+s13)*(s(j2,j8)+s(j3,j8))-s23*s(j1,j8)))/(s13*s23)
      radi_ww=fac*radi_ww
      return
      end

