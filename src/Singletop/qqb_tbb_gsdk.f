      subroutine qqb_tbb_gsdk(p,msqc)
      implicit none
c     Subtraction Matrix element for real corrections to decay in 
C     single top production
C     (nwz=+1)
c      u(-p1)+dbar(-p2)-->t(=> n(p3)+e^+(p4)+b(p5)+g(p7))+F(p6)
C     or for
C     (nwz=-1)
c      ubar(-p1)+d(-p2)-->t~(=> e^-(p3)+n(p4)+bbar(p5)+g(p7))+F(p6)
C     averaged(summed) over initial(final) colours and spins

c--- g(p7) represents a gluon 


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

      call wtransform(p,q,pbDpg,ptDpg,ptDpb)
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

      call qqb_tbb(q,msq) 
      fac=gsq*cf*(1d0/pbDpg*(2d0/omz-1d0-z)-(mt/ptDpg)**2)

      do j=-nf,nf
      do k=-nf,nf
         msqc(1,j,k)=fac*msq(j,k)
      enddo
      enddo
      return
      end

      subroutine wtransform(p,q,pbDpg,ptDpg,ptDpb)
      implicit none
      include 'constants.f'
      double precision p(mxpart,4),pw(4),pt(4),lDt(3:4),lDw(3:4)
      double precision ptDpt,pwDpw,ptDpw,q(mxpart,4),root,hsin,hcos,a,b
      double precision ptDpg,pbDpg,ptDpb

      integer j,nu
      do nu=1,4
         pw(nu)=p(3,nu)+p(4,nu)
         pt(nu)=pw(nu)+p(5,nu)+p(7,nu)
         do j=1,6
            q(j,nu)=p(j,nu)
         enddo
      enddo
      pbDpg=p(5,4)*p(7,4)-p(5,1)*p(7,1)-p(5,2)*p(7,2)-p(5,3)*p(7,3)
      ptDpg=pt(4)*p(7,4)-pt(1)*p(7,1)-pt(2)*p(7,2)-pt(3)*p(7,3)
      ptDpb=pt(4)*p(5,4)-pt(1)*p(5,1)-pt(2)*p(5,2)-pt(3)*p(5,3)
      ptDpw=pt(4)*pw(4)-pt(1)*pw(1)-pt(2)*pw(2)-pt(3)*pw(3)
      ptDpt=pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2
      pwDpw=pw(4)**2-pw(1)**2-pw(2)**2-pw(3)**2
      root=sqrt(ptDpw**2-ptDpt*pwDpw)
      hsin=0.5d0/(ptDpt*pwDpw)*(-(ptDpt-pwDpw)*ptDpw+(ptDpt+pwDpw)*root)
      hcos=0.5d0/(ptDpt*pwDpw)*(+(ptDpt+pwDpw)*ptDpw-(ptDpt-pwDpw)*root)
C---calculate coefficients of lorentz transformation
      a=hsin/root
      b=(hcos-1d0)/root**2
c---dot t and w into decay products of w
      do j=3,4
         lDt(j)=p(j,4)*pt(4)-p(j,1)*pt(1)-p(j,2)*pt(2)-p(j,3)*pt(3)
         lDw(j)=p(j,4)*pw(4)-p(j,1)*pw(1)-p(j,2)*pw(2)-p(j,3)*pw(3)
      enddo
      do nu=1,4
      do j=3,4
      q(j,nu)=p(j,nu)+a*(pt(nu)*lDw(j)-pw(nu)*lDt(j))
     . +b*(ptDpw*(pt(nu)*ldw(j)+pw(nu)*ldt(j))
     .  -pwDpw*lDt(j)*pt(nu)-ptDpt*lDw(j)*pw(nu))
      enddo
      q(5,nu)=-q(1,nu)-q(2,nu)-q(3,nu)-q(4,nu)-q(6,nu)
      enddo
      call storeptilde(1,q)
      return
      end

