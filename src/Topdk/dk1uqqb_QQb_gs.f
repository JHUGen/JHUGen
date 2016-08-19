      subroutine dk1uqqb_QQb_gs(p,msqc)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     November, 2011.                                                  *
*     calculate the subtraction term for radiation in the              *
*     top quark decay for the process                                  *
*                                                                      *
*     q(-p1)+qbar(-p2) = nu(p3)+e+(p4)+b(p5)                           *
*                        +bbar(p6)+e-(p7)+nubar(p8)+g(p9)              *
*                                                                      *
*     Top is kept strictly on-shell although all spin correlations     *
*     are retained. B-quark is taken to be either massive or massless. *
*                                                                      *
************************************************************************
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

      ndmax=1
      
      do j=-nf,nf
      do k=-nf,nf
        msqc(1,j,k)=0d0
        incldip(1)=.true.
      enddo
      enddo

c--- special dipole for radiation in top decay
      call wtransform_generic(p,3,4,5,9,q,pbDpg,ptDpg,ptDpb)

      pwsq=2d0*(q(3,4)*q(4,4)-q(3,1)*q(4,1)-q(3,2)*q(4,2)-q(3,3)*q(4,3))

c--- form of subtraction depends on whether b-quark in decay is massless or not
      if (mb .lt. 1d-6) then
c----- massless case      
        omz=ptDpg/(ptDpb+ptDpg-pbDpg)
        z=1d0-omz
        xr=dsqrt(pwsq/mt**2)
        ymax=(1d0+xr)**2*z*omz/(z+xr**2*omz)
        y=2d0*pbDpg/mt**2/(1d0-xr)**2
        if ((z .lt. 1d0-aff) .and. (y .gt. aff*ymax)) then
          incldip(1)=.false.
          return
        endif
        fac=gsq*cf*(1d0/pbDpg*(2d0/omz-1d0-z)-(mt/ptDpg)**2)
      else
c----- massive case    
c-----  (no alpha-dependence at present)  
        fac=gsq*cf*((mt**2+mb**2-pwsq)/(ptDpg*pbDpg)
     &             -(mt/ptDpg)**2-(mb/pbDpg)**2)
      endif      
      
      call qqb_QQbdku(q,msq) 

      do j=-nf,nf
      do k=-nf,nf
      msqc(1,j,k)=fac*msq(j,k)
      enddo
      enddo
                        
      return
      end

