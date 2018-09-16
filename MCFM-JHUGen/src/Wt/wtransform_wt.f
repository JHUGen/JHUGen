      subroutine wtransform_wt(p,q,pbDpg,ptDpg,ptDpb)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),pw(4),pt(4),lDt(5:6),lDw(5:6)
      real(dp):: ptDpt,pwDpw,ptDpw,q(mxpart,4),root,hsin,hcos,a,b
      real(dp):: ptDpg,pbDpg,ptDpb

      integer:: j,nu

      do nu=1,4
         pw(nu)=p(5,nu)+p(6,nu)
         pt(nu)=pw(nu)+p(7,nu)+p(8,nu)
         do j=1,4
            q(j,nu)=p(j,nu)
         enddo
      enddo

      pbDpg=p(7,4)*p(8,4)-p(7,1)*p(8,1)-p(7,2)*p(8,2)-p(7,3)*p(8,3)
      ptDpg=pt(4)*p(8,4)-pt(1)*p(8,1)-pt(2)*p(8,2)-pt(3)*p(8,3)
      ptDpb=pt(4)*p(7,4)-pt(1)*p(7,1)-pt(2)*p(7,2)-pt(3)*p(7,3)
      ptDpw=pt(4)*pw(4)-pt(1)*pw(1)-pt(2)*pw(2)-pt(3)*pw(3)
      ptDpt=pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2
      pwDpw=pw(4)**2-pw(1)**2-pw(2)**2-pw(3)**2

      root=sqrt(ptDpw**2-ptDpt*pwDpw)
      hsin=half/(ptDpt*pwDpw)*(-(ptDpt-pwDpw)*ptDpw+(ptDpt+pwDpw)*root)
      hcos=half/(ptDpt*pwDpw)*(+(ptDpt+pwDpw)*ptDpw-(ptDpt-pwDpw)*root)

C---calculate coefficients of lorentz transformation
      a=hsin/root
      b=(hcos-1._dp)/root**2

c---dot t and w into decay products of w
      do j=5,6
         lDt(j)=p(j,4)*pt(4)-p(j,1)*pt(1)-p(j,2)*pt(2)-p(j,3)*pt(3)
         lDw(j)=p(j,4)*pw(4)-p(j,1)*pw(1)-p(j,2)*pw(2)-p(j,3)*pw(3)
      enddo

      do nu=1,4
      do j=5,6
      q(j,nu)=p(j,nu)+a*(pt(nu)*lDw(j)-pw(nu)*lDt(j))
     & +b*(ptDpw*(pt(nu)*ldw(j)+pw(nu)*ldt(j))
     &  -pwDpw*lDt(j)*pt(nu)-ptDpt*lDw(j)*pw(nu))
      enddo
      q(7,nu)=-q(1,nu)-q(2,nu)-q(3,nu)-q(4,nu)-q(5,nu)-q(6,nu)
      q(8,nu)=0._dp
      enddo

      call storeptilde(1,q)

      return
      end

