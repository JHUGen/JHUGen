      subroutine wtransform_generic(p,n3,n4,n5,nmax,q,pbDpg,ptDpg,ptDpb)
      implicit none
      include 'types.f'
c---  Author: J. Campbell, Nov. 2011
c---  input momenta p, top decay products W(->n3+n4) + b(n5)
c---  nmax is the gluon momentum label and the last momentum considered
c---  output momenta q and dot products pbDpg,ptDpg,ptDpb
c---  q == p for all momenta except n3,n4,n5 and q(nmax) absent
c---- form of Lorentz transformation given in Section IV of arxiv:hep-ph/0408158
c---- Updated: February 24, 2012 to allow for a non-zero b-quark mass
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),pw(4),pt(4),lDt(12),lDw(12)
      real(dp):: ptDpt,pwDpw,ptDpw,q(mxpart,4),root,hsin,hcos,a,b
      real(dp):: ptDpg,pbDpg,ptDpb,pbDpb,rtlambda
      integer:: j,nu,n3,n4,n5,nmax
      
      if (nmax > 12) then
        write(6,*) 'Routine wtransform_generic fails, limit nmax=12'
        stop
      endif
      
      if (n3> n4) then
        write(6,*) 'Routine wtransform_generic fails, n3>n4',n3,n4
        stop
      endif
      
      do nu=1,4
        pw(nu)=p(n3,nu)+p(n4,nu)
        pt(nu)=pw(nu)+p(n5,nu)+p(nmax,nu)
        do j=1,nmax-1
          q(j,nu)=p(j,nu)
        enddo
      enddo
      
      pbDpg=p(n5,4)*p(nmax,4)-p(n5,1)*p(nmax,1)
     &     -p(n5,2)*p(nmax,2)-p(n5,3)*p(nmax,3)
      ptDpg=pt(4)*p(nmax,4)-pt(1)*p(nmax,1)
     &     -pt(2)*p(nmax,2)-pt(3)*p(nmax,3)
      ptDpb=pt(4)*p(n5,4)-pt(1)*p(n5,1)-pt(2)*p(n5,2)-pt(3)*p(n5,3)
      ptDpw=pt(4)*pw(4)-pt(1)*pw(1)-pt(2)*pw(2)-pt(3)*pw(3)
      ptDpt=pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2
      pwDpw=pw(4)**2-pw(1)**2-pw(2)**2-pw(3)**2
      pbDpb=p(n5,4)**2-p(n5,1)**2-p(n5,2)**2-p(n5,3)**2
      root=sqrt(ptDpw**2-ptDpt*pwDpw)
      rtlambda=sqrt((ptDpt-pwDpw-pbDpb)**2-4._dp*pwDpw*pbDpb)
      hsin=0.5_dp/(ptDpt*pwDpw)*(-rtlambda*ptDpw
     &                          +(ptDpt+pwDpw-pbDpb)*root)
      hcos=0.5_dp/(ptDpt*pwDpw)*(+(ptDpt+pwDpw-pbDpb)*ptDpw
     &                          -rtlambda*root)
C---calculate coefficients of lorentz transformation
      a=hsin/root
      b=(hcos-1._dp)/root**2
c---dot t and w into decay products of w
      do j=n3,n4
         lDt(j)=p(j,4)*pt(4)-p(j,1)*pt(1)-p(j,2)*pt(2)-p(j,3)*pt(3)
         lDw(j)=p(j,4)*pw(4)-p(j,1)*pw(1)-p(j,2)*pw(2)-p(j,3)*pw(3)
      enddo

      do nu=1,4
        do j=n3,n4
          q(j,nu)=p(j,nu)+a*(pt(nu)*lDw(j)-pw(nu)*lDt(j))
     &     +b*(ptDpw*(pt(nu)*ldw(j)+pw(nu)*ldt(j))
     &     -pwDpw*lDt(j)*pt(nu)-ptDpt*lDw(j)*pw(nu))
        enddo
        q(n5,nu)=-q(1,nu)-q(2,nu)
        do j=3,nmax-1
          if (j .ne. n5) q(n5,nu)=q(n5,nu)-q(j,nu)
        enddo      
      enddo
      
      call storeptilde(1,q)
      
      return
      end

