      subroutine dkqqb_w_twdk_gs(p,msqc)
      implicit none
      include 'types.f'


c--- January, 2005.
c--- Subtraction terms for radiation in top decay
c--- matrix element squared and averaged over initial colours and spins
c     q(-p1) + qbar(-p2) --> W + t(p5678)
c                            |   |
c                            |   --> nu(p5) + e^+(p6) + b(p7) + f(p8)
c                            |
c                            --> e^-(p3) + nubar(p4)


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
      msqc(1,j,k)=zero
      enddo
      enddo

      ndmax=1
      incldip(1)=.true.

      call wtransform_wt(p,q,pbDpg,ptDpg,ptDpb)
      omz=ptDpg/(ptDpb+ptDpg-pbDpg)
      z=one-omz
      pwsq=two*(q(3,4)*q(4,4)-q(3,1)*q(4,1)-q(3,2)*q(4,2)-q(3,3)*q(4,3))
      xr=sqrt(pwsq/mt**2)
      ymax=(one+xr)**2*z*omz/(z+xr**2*omz)
      y=two*pbDpg/mt**2/(one-xr)**2
      if ((z < one-aff) .and. (y > aff*ymax)) then
        incldip(1)=.false.
        return
      endif

      call qqb_w_twdk(q,msq)
      fac=gsq*cf*(one/pbDpg*(two/omz-one-z)-(mt/ptDpg)**2)


      do j=-nf,nf,nf
      do k=-nf,nf,nf
      msqc(1,j,k)=fac*msq(j,k)
      enddo
      enddo

      return
      end

