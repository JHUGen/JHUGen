      subroutine dkqqb_w_twdk_gs(p,msqc)
      implicit none

c--- January, 2005.
c--- Subtraction terms for radiation in top decay
c--- matrix element squared and averaged over initial colours and spins
c     q(-p1) + qbar(-p2) --> W + t(p5678)
c                            |   |
c                            |   --> nu(p5) + e^+(p6) + b(p7) + f(p8)
c                            |
c                            --> e^-(p3) + nubar(p4)


      include 'constants.f'
      include 'masses.f'
      include 'ptilde.f'
      include 'qcdcouple.f'
      include 'alfacut.f'
      include 'incldip.f'
      double precision msq(-nf:nf,-nf:nf),msqc(maxd,-nf:nf,-nf:nf),
     . p(mxpart,4),q(mxpart,4),omz,z,fac,ptDpg,pbDpg,ptDpb,pwsq,xr,
     . y,ymax
      integer j,k

      do j=-nf,nf
      do k=-nf,nf
      msqc(1,j,k)=0d0
      enddo
      enddo

      ndmax=1
      incldip(1)=.true.

      call wtransform_wt(p,q,pbDpg,ptDpg,ptDpb)
      omz=ptDpg/(ptDpb+ptDpg-pbDpg)
      z=1d0-omz
      pwsq=2d0*(q(3,4)*q(4,4)-q(3,1)*q(4,1)-q(3,2)*q(4,2)-q(3,3)*q(4,3))
      xr=dsqrt(pwsq/mt**2)
      ymax=(1d0+xr)**2*z*omz/(z+xr**2*omz)
      y=2d0*pbDpg/mt**2/(1d0-xr)**2
      if ((z .lt. 1d0-aff) .and. (y .gt. aff*ymax)) then
        incldip(1)=.false.
        return
      endif

      call qqb_w_twdk(q,msq)
      fac=gsq*cf*(1d0/pbDpg*(2d0/omz-1d0-z)-(mt/ptDpg)**2)


      do j=-nf,nf,nf
      do k=-nf,nf,nf
      msqc(1,j,k)=fac*msq(j,k)
      enddo
      enddo

      return
      end

