      subroutine qqb_wh_g(P,msq)
c---Matrix element squared averaged over initial colors and spins
c---for nwz=1
c     q(-p1)+qbar(-p2) -->  H  + W +g(p7)
c                           |    |
c                           |    --> nu(p3)+e^+(p4)
c                           |
c                           ---> b(p5)+b(p6)
c---for nwz=-1
c     q(-p1)+qbar(-p2) -->  H  + W +g(p7)
c                           |    |
c                           |    --> e^-(p3)+nubar(p4)
c                           |
c                           ---> b(p5)+b(p6)
c   for the moment --- radiation only from initial line
      implicit none 
      include 'constants.f'
      include 'ckm.f'
      include 'sprods_com.f'
      integer j,k
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf)
      double precision radi
c      double precision radf
      double precision qqbWHg,qbqWHg,qgWHq,gqWHq,gqbWHqb,qbgWHqb

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call dotem(7,p,s)

      qqbWHg=aveqq*radi(1,2,7,5,6,3,4)
c    &  +zip*aveqq*radf(1,2,7,5,6,3,4)
      qbqWHg=aveqq*radi(2,1,7,5,6,3,4)
c    &   +zip*aveqq*radf(2,1,7,5,6,3,4)
c---turn off radiation from final legs
      qgWHq=-radi(1,7,2,5,6,3,4)*aveqg
      gqWHq=-radi(2,7,1,5,6,3,4)*aveqg

      gqbWHqb=-radi(7,2,1,5,6,3,4)*aveqg
      qbgWHqb=-radi(7,1,2,5,6,3,4)*aveqg

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


      double precision function radi(j1,j2,j3,j4,j5,j6,j7)
      implicit none 
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      integer j1,j2,j3,j4,j5,j6,j7
      double precision s45,s12,s13,s23,s123,prop
      double precision fac,hdecay,msqhbb


      s45=s(j4,j5)+2d0*mb**2
      s12=s(j1,j2)
      s13=s(j1,j3)
      s23=s(j2,j3)
      s123=s12+s13+s23
c---calculate the 2 W propagators
      prop=       ((s123-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(j6,j7)-wmass**2)**2+(wmass*wwidth)**2)
      
      fac=2d0*cf*xn*gsq*gwsq**3*wmass**2/prop
c      hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s45-4d0*mb**2)
      hdecay=msqhbb(s45)
      hdecay=hdecay/((s45-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay
c-- Old form of this matrix element (modified to facilitate extension
c--- to H->WW decay)
c      fac=CF*xnsq*gsq*gw**8*mbsq*(s45-4d0*mb**2)/prop
      radi=s12/s13/s23
     & *(2d0*s(j1,j7)*s(j2,j6)+s(j1,j7)*s(j3,j6)+s(j2,j6)*s(j3,j7))
     & +(s(j1,j7)*s(j2,j6)+s(j2,j6)*s(j3,j7)-s(j1,j6)*s(j1,j7))/s13
     & +(s(j1,j7)*s(j2,j6)+s(j1,j7)*s(j3,j6)-s(j2,j6)*s(j2,j7))/s23
      radi=fac*radi
      return
      end

c      double precision function radf(j1,j2,j3,j4,j5,j6,j7)
c      implicit none 
c      include 'constants.f'
c      include 'ewcouple.f'
c      include 'qcdcouple.f'
c      include 'masses.f'
c      include 'sprods_com.f'
c      integer j1,j2,j3,j4,j5,j6,j7
c      double precision s34,s35,s45,s12,s67,s345,prop
c      double precision fac


c      s12=s(j1,j2)
c      s67=s(j6,j7)
c      s34=s(j3,j4)
c      s35=s(j3,j5)
c      s45=s(j4,j5)
c      s345=s34+s35+s45+2d0*mb**2

c---calculate the 3 propagators
c      prop=     ((s12-wmass**2)**2+(wmass*wwidth)**2)
c      prop=prop*((s67-wmass**2)**2+(wmass*wwidth)**2)
c      prop=prop*((s345-hmass**2)**2+(hmass*hwidth)**2)
c      
c      fac=two*CF*xnsq*gsq*gw**8*mbsq*s(j1,j7)*s(j2,j6)/prop
c      radf=(s45/(s34*s35)-mb**2*(1d0/s34**2+1d0/s35**2))*(s45-2d0*mb**2)
c     & +0.5d0*((s34+s35)**2+2*(s45-mb**2)*(s34+s35))/s34/s35
c     & -mb**2*(s35/s34**2+s34/s35**2) 
c      radf=fac*radf
c      return
c      end


