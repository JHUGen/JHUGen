      subroutine dkqg_tbqdk_gs(p,msqc)
      implicit none
      include 'types.f'

************************************************************************
*     Author: J. Campbell, January 2012                                *
*                                                                      *
*     Calculates the subtraction term for radiation in the             *
*     top quark decay for the process                                  *
*                                                                      *
*     q(p1) + g(p2) -> t(e(p3)+nu(p4)+b(p5)+g(p8)) + b(p6) + q'(p7)    *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained.                                                    *
*                                                                      *
************************************************************************
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
     & y,ymax,dot
      integer:: j,k

      do j=-nf,nf
      do k=-nf,nf
      msqc(1,j,k)=0._dp
      enddo
      enddo

      ndmax=1
      incldip(1)=.true.

      call wtransform_generic(p,3,4,5,8,q,pbDpg,ptDpg,ptDpb)

      pwsq=2._dp*(q(3,4)*q(4,4)-q(3,1)*q(4,1)-q(3,2)*q(4,2)-q(3,3)*q(4,3))

c--- form of subtraction depends on whether b-quark in decay is massless or not
      if (abs(dot(p,5,5)) < 1.e-6_dp) then
c----- massless case
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)
        z=1._dp-omz
        xr=sqrt(pwsq/mt**2)
        ymax=(1._dp+xr)**2*z*omz/(z+xr**2*omz)
        y=2._dp*pbDpg/mt**2/(1._dp-xr)**2
        if ((z < 1._dp-aff) .and. (y > aff*ymax)) then
          incldip(1)=.false.
          return
        endif
        fac=gsq*cf*(1._dp/pbDpg*(2._dp/omz-1._dp-z)-(mt/ptDpg)**2)
      else
c----- massive case
c-----  (no alpha-dependence at present)
        fac=gsq*cf*((mt**2+mb**2-pwsq)/(ptDpg*pbDpg)
     &             -(mt/ptDpg)**2-(mb/pbDpg)**2)
      endif

      call qg_tbqdk(q,msq)

      do j=-nf,nf
      do k=-nf,nf
      msqc(1,j,k)=fac*msq(j,k)
      enddo
      enddo

      return
      end

