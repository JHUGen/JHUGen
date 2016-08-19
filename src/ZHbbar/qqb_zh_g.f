      subroutine qqb_zh_g(P,msq)
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  H  + Z +g(p7)
c                           |    |
c                           |    --> fermion(p3)+antifermion(p4)
c                           |
c                           ---> b(p5)+b(p6)
c   for the moment --- radiation only from initial line
      implicit none 
      include 'constants.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      integer j,k
      double precision P(mxpart,4),msq(-nf:nf,-nf:nf)
      double precision radiLL
      double precision qqbZHgL,qbqZHgL,qgZHqL,gqZHqL,gqbZHqbL,qbgZHqbL
      double precision qqbZHgR,qbqZHgR,qgZHqR,gqZHqR,gqbZHqbR,qbgZHqbR


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call dotem(7,p,s)

      qqbZHgL=aveqq*radiLL(1,2,7,5,6,3,4)
      qqbZHgR=aveqq*radiLL(1,2,7,5,6,4,3)
      qbqZHgL=aveqq*radiLL(2,1,7,5,6,3,4)
      qbqZHgR=aveqq*radiLL(2,1,7,5,6,4,3)

      qgZHqL=-radiLL(1,7,2,5,6,3,4)*aveqg
      qgZHqR=-radiLL(1,7,2,5,6,4,3)*aveqg
      gqZHqL=-radiLL(2,7,1,5,6,3,4)*aveqg
      gqZHqR=-radiLL(2,7,1,5,6,4,3)*aveqg

      gqbZHqbL=-radiLL(7,2,1,5,6,3,4)*aveqg
      gqbZHqbR=-radiLL(7,2,1,5,6,4,3)*aveqg

      qbgZHqbL=-radiLL(7,1,2,5,6,3,4)*aveqg
      qbgZHqbR=-radiLL(7,1,2,5,6,4,3)*aveqg

      do j=-nf,nf
      do k=-nf,nf

      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 40
      if ((j .eq. 0) .and. (k .eq. 0)) goto 40

      if     ((j .gt. 0) .and. (k .lt. 0)) then
          msq(j,k)=((L(j)*l1)**2+(R(j)*r1)**2)*qqbZHgL
     .            +((R(j)*l1)**2+(L(j)*r1)**2)*qqbZHgR
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=((L(k)*l1)**2+(R(k)*r1)**2)*qbqZHgL
     .            +((R(k)*l1)**2+(L(k)*r1)**2)*qbqZHgR
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=((L(j)*l1)**2+(R(j)*r1)**2)*qgZHqL
     .            +((R(j)*l1)**2+(L(j)*r1)**2)*qgZHqR
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=((L(-j)*l1)**2+(R(-j)*r1)**2)*qbgZHqbL
     .            +((R(-j)*l1)**2+(L(-j)*r1)**2)*qbgZHqbR
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=((L(k)*l1)**2+(R(k)*r1)**2)*gqZHqL
     .            +((R(k)*l1)**2+(L(k)*r1)**2)*gqZHqR
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=((L(-k)*l1)**2+(R(-k)*r1)**2)*gqbZHqbL
     .            +((R(-k)*l1)**2+(L(-k)*r1)**2)*gqbZHqbR
      endif
 40   continue
      enddo
      enddo
      return
      end


      double precision function radiLL(j1,j2,j3,j4,j5,j6,j7)
      implicit none 
      include 'constants.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      integer j1,j2,j3,j4,j5,j6,j7
      double precision s45,s12,s13,s23,s123,prop
      double precision fac,hdecay,msqhbb


      s45=s(j4,j5)+2d0*mb**2
      s12=s(j1,j2)
      s13=s(j1,j3)
      s23=s(j2,j3)
      s123=s12+s13+s23
c---calculate the 2 Z propagators
      prop=         ((s123-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(j6,j7)-zmass**2)**2+(zmass*zwidth)**2)

      fac=8d0*cf*xn*(xw/(1d0-xw))**2*gsq*gwsq**3*wmass**2/prop
C      hdecay=xn*gwsq*mbsq/(4d0*wmass**2)*2d0*(s45-4d0*mb**2)
      hdecay=msqhbb(s45)

      hdecay=hdecay/((s45-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay

      radiLL=s12/s13/s23
     & *(2d0*s(j1,j7)*s(j2,j6)+s(j1,j7)*s(j3,j6)+s(j2,j6)*s(j3,j7))
     & +(s(j1,j7)*s(j2,j6)+s(j2,j6)*s(j3,j7)-s(j1,j6)*s(j1,j7))/s13
     & +(s(j1,j7)*s(j2,j6)+s(j1,j7)*s(j3,j6)-s(j2,j6)*s(j2,j7))/s23

      radiLL=fac*radiLL
      return
      end
