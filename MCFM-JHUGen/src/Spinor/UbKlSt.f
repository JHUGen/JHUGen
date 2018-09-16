      subroutine UbKlSt(q,m,aux,i,f)
      implicit none
      include 'types.f'
      include 'constants.f'
c     subroutine for massive ubar spinor alla Kleiss and Stirling
C     q complex vector q(4),q(1),q(2),q(3)=E_q,qx,qy,qz
C     'aux' massless auxiliary complex vector
c     i polarization 1 or -1
C     f output spinor

      integer:: i,j
      real(dp):: m
      complex(dp):: q(4),qt(4),aux(4),f(4),fqt(4),eta(4),den,
     & qDaux,qsq

      if (m == zip) then
        call ubarspinor0(q,i,f)
        return
      else
        qDaux=q(4)*aux(4)-q(1)*aux(1)-q(2)*aux(2)-q(3)*aux(3)
        qsq=q(4)**2-q(1)**2-q(2)**2-q(3)**2
        qt(:)=q(:)-half*qsq/qDaux*aux(:)
        call ubarspinor0(aux,-i,eta)
        call uspinor0(qt,i,fqt)
        den=eta(1)*fqt(1)+eta(2)*fqt(2)+eta(3)*fqt(3)+eta(4)*fqt(4)
      endif

      call Ubkslash(eta,q,f)
      do j=1,4
      f(j)=(f(j)+eta(j)*m)/den
      enddo
      return
      end


