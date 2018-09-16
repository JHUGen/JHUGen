      subroutine bq_tpq_gsdk(p,msqc)
      implicit none
      include 'types.f'


c     Matrix element for t-bbar production with radiation in decay
c      b(-p1)+u(-p2)-->t(n(p3)+e^+(p4)+b(p5)+g(p7))+d(p6)
C     averaged(summed) over initial(final) colours and spins
c--- g(p7) represents a gluon


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ptilde.f'
      include 'qcdcouple.f'
      include 'alfacut.f'
      include 'incldip.f'
      real(dp):: msq(-nf:nf,-nf:nf),msqc(maxd,-nf:nf,-nf:nf),
     & p(mxpart,4),q(mxpart,4),omz,z,fac,ptDpg,pbDpg,ptDpb,pwsq,xr,
     & y,ymax
      integer:: j,k

      do j=-nf,nf
      do k=-nf,nf
      msqc(1,j,k)=0._dp
      enddo
      enddo

      ndmax=1
      incldip(1)=.true.

      call wtransform(p,q,pbDpg,ptDpg,ptDpb)
      omz=ptDpg/(ptDpb+ptDpg-pbDpg)
      z=1._dp-omz
      pwsq=2._dp*(q(3,4)*q(4,4)-q(3,1)*q(4,1)-q(3,2)*q(4,2)-q(3,3)*q(4,3))
      xr=sqrt(pwsq/mt**2)
      ymax=(1._dp+xr)**2*z*omz/(z+xr**2*omz)
      y=2._dp*pbDpg/mt**2/(1._dp-xr)**2
      if ((z < 1._dp-aff) .and. (y > aff*ymax)) then
        incldip(1)=.false.
        return
      endif

      call bq_tpq(q,msq)
      fac=gsq*cf*(1._dp/pbDpg*(2._dp/omz-1._dp-z)-(mt/ptDpg)**2)

      do j=-nf,nf
      do k=-nf,nf
      msqc(1,j,k)=fac*msq(j,k)
      enddo
      enddo
      return
      end

