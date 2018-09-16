      subroutine qqb_wh_gaga_g(P,msq)
      implicit none
      include 'types.f'
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
c---- Extension to photon decay contributed by Fabian Stoeckli
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ckm.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: radi_gaga
c      real(dp):: radf
      real(dp):: qqbWHg,qbqWHg,qgWHq,gqWHq,gqbWHqb,qbgWHqb

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(7,p,s)
c      if (
c     &      (s(5,6) < four*mbsq) 
c     & .or. (s(1,5)*s(2,5)/s(1,2) < mbsq) 
c     & .or. (s(1,6)*s(2,6)/s(1,2) < mbsq) ) return

      qqbWHg=aveqq*radi_gaga(1,2,7,5,6,3,4)
c    &  +zip*aveqq*radf(1,2,7,5,6,3,4)
      qbqWHg=aveqq*radi_gaga(2,1,7,5,6,3,4)
c    &   +zip*aveqq*radf(2,1,7,5,6,3,4)
c---turn off radiation from final legs
      qgWHq=-radi_gaga(1,7,2,5,6,3,4)*aveqg
      gqWHq=-radi_gaga(2,7,1,5,6,3,4)*aveqg

      gqbWHqb=-radi_gaga(7,2,1,5,6,3,4)*aveqg
      qbgWHqb=-radi_gaga(7,1,2,5,6,3,4)*aveqg

c      write(6,*) 'qqbWHg',qqbWHg
c      write(6,*) 'qbqWHg',qbqWHg
c      write(6,*) 'qbgWHqb',qbgWHqb
c      write(6,*) 'gqbWHqb',gqbWHqb
c      write(6,*) 'qgWHq',qgWHq
c      write(6,*) 'gqWHq',gqWHq


      do j=-nf,nf
      do k=-nf,nf

      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=Vsq(j,k)*qqbWHg
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=Vsq(j,k)*qbqWHg
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=
     &   (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWHq
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=
     &    (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWHqb
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=
     &    (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWHq
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=
     &    (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWHqb
      endif

      enddo
      enddo
      return
      
      end


      function radi_gaga(j1,j2,j3,j4,j5,j6,j7)
      implicit none
      include 'types.f'
      real(dp):: radi_gaga
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5,j6,j7
      real(dp):: s45,s12,s13,s23,s123,prop
      real(dp):: fac,hdecay
      real(dp):: msqhgamgam


      s45=s(j4,j5)+2._dp*mb**2
      s12=s(j1,j2)
      s13=s(j1,j3)
      s23=s(j2,j3)
      s123=s12+s13+s23
c---calculate the 2 W propagators
      prop=       ((s123-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((s(j6,j7)-wmass**2)**2+(wmass*wwidth)**2)
      
      fac=2._dp*cf*xn*gsq*gwsq**3*wmass**2/prop

      hdecay=msqhgamgam(s45)/((s45-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay
      
c-- Old form of this matrix element (modified to facilitate extension
c--- to H->WW decay)
c      fac=CF*xnsq*gsq*gw**8*mbsq*(s45-4._dp*mb**2)/prop
      radi_gaga=s12/s13/s23
     & *(2._dp*s(j1,j7)*s(j2,j6)+s(j1,j7)*s(j3,j6)+s(j2,j6)*s(j3,j7))
     & +(s(j1,j7)*s(j2,j6)+s(j2,j6)*s(j3,j7)-s(j1,j6)*s(j1,j7))/s13
     & +(s(j1,j7)*s(j2,j6)+s(j1,j7)*s(j3,j6)-s(j2,j6)*s(j2,j7))/s23
      radi_gaga=fac*radi_gaga
      return
      end

c      function radf(j1,j2,j3,j4,j5,j6,j7)
c      implicit none
c      include 'types.f'
c      real(dp):: radf
c       
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'ewcouple.f'
c      include 'qcdcouple.f'
c      include 'masses.f'
c      include 'sprods_com.f'
c      integer:: j1,j2,j3,j4,j5,j6,j7
c      real(dp):: s34,s35,s45,s12,s67,s345,prop
c      real(dp):: fac


c      s12=s(j1,j2)
c      s67=s(j6,j7)
c      s34=s(j3,j4)
c      s35=s(j3,j5)
c      s45=s(j4,j5)
c      s345=s34+s35+s45+2._dp*mb**2

c---calculate the 3 propagators
c      prop=     ((s12-wmass**2)**2+(wmass*wwidth)**2)
c      prop=prop*((s67-wmass**2)**2+(wmass*wwidth)**2)
c      prop=prop*((s345-hmass**2)**2+(hmass*hwidth)**2)
c      
c      fac=two*CF*xnsq*gsq*gw**8*mbsq*s(j1,j7)*s(j2,j6)/prop
c      radf=(s45/(s34*s35)-mb**2*(1._dp/s34**2+1._dp/s35**2))*(s45-2._dp*mb**2)
c     & +0.5_dp*((s34+s35)**2+2*(s45-mb**2)*(s34+s35))/s34/s35
c     & -mb**2*(s35/s34**2+s34/s35**2) 
c      radf=fac*radf
c      return
c      end


