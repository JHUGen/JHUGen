      subroutine qqb_ttz(pin,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     May, 2011.                                                       *
*     Calculate the element squared for the process                    *
*                                                                      *
*     q(-p1) +qbar(-p2) -> tbar(bbar(p6)+e-(p7)+nubar(p8))             *
*                         +t(b(p5)+nu(p3)+e+(p4))                      *
*                         +Z(l(p9)+a(p10))                             *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'zprods_com.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'topzlabels.f'
      include 'kprocess.f'

      integer:: j,k,nu
      real(dp):: msq(-nf:nf,-nf:nf),pin(mxpart,4),p(mxpart,4)
      real(dp):: pw1(4),pw2(4),q(4),a(4),u(4),b(4),z(4),
     & q1(4),q2(4),a1(4),a2(4),z1(4),z2(4)
      real(dp):: sz,sw1,sw2,qDq,aDa,uDu,bDb,densq,
     & p3Dp5,p6Dp8,q1Dq1,q2Dq2,a1Da1,a2Da2
      real(dp):: wtqqb(2),wtqbq(2),wtgg
      real(dp):: gamu1,gamu2,gamb1,gamb2,gamq4,gama7,dot
      real(dp):: facqqb,facgg,denu,denb,denq1,denq2,dena1,dena2
      real(dp):: p1Du,p2Du,p1Db,p2Db,p4Dq,p7Da,denz1,denz2
      complex(dp):: propz
      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

c      write(6,*) 'mt',mt,twidth,wmass,wwidth,zmass,zwidth
c      write(6,*) 'Enter scalefac**2'
c      read(5,*)  scalefac
c      scalefac=sqrt(scalefac)

c      mt=scalefac*mt
c      twidth=scalefac*twidth
c      wmass=scalefac*wmass
c      wwidth=scalefac*wwidth
c      zmass=scalefac*zmass
c      zwidth=scalefac*zwidth

      do nu=1,4
      do j=1,mxpart
c      p(j,nu)=scalefac*pin(j,nu)
      p(j,nu)=pin(j,nu)
      enddo

      z(nu)=p(9,nu)+p(10,nu)
      pw1(nu)=p(3,nu)+p(4,nu)
      pw2(nu)=p(7,nu)+p(8,nu)

      q(nu)=+p(3,nu)+p(4,nu)+p(5,nu)
      u(nu)=q(nu)+z(nu)
      q1(nu)=q(nu)+p(1,nu)
      q2(nu)=q(nu)+p(2,nu)

      a(nu)=-p(6,nu)-p(7,nu)-p(8,nu)
      b(nu)=a(nu)-z(nu)
      a1(nu)=a(nu)-p(1,nu)
      a2(nu)=a(nu)-p(2,nu)
      z1(nu)=z(nu)+p(1,nu)
      z2(nu)=-z(nu)-p(2,nu)
      enddo


      sz=(z(4)**2-z(1)**2-z(2)**2-z(3)**2)
      sw1=(pw1(4)**2-pw1(1)**2-pw1(2)**2-pw1(3)**2)
      sw2=(pw2(4)**2-pw2(1)**2-pw2(2)**2-pw2(3)**2)
      qDq=(q(4)**2-q(1)**2-q(2)**2-q(3)**2)
      q1Dq1=(q1(4)**2-q1(1)**2-q1(2)**2-q1(3)**2)
      q2Dq2=(q2(4)**2-q2(1)**2-q2(2)**2-q2(3)**2)
      aDa=(a(4)**2-a(1)**2-a(2)**2-a(3)**2)
      a1Da1=(a1(4)**2-a1(1)**2-a1(2)**2-a1(3)**2)
      a2Da2=(a2(4)**2-a2(1)**2-a2(2)**2-a2(3)**2)
      uDu=(u(4)**2-u(1)**2-u(2)**2-u(3)**2)
      bDb=(b(4)**2-b(1)**2-b(2)**2-b(3)**2)

      p3Dp5=dot(p,3,5)
      p6Dp8=dot(p,6,8)

      densq=((qDq-mt**2)**2+mt**2*twidth**2)
     &     *((aDa-mt**2)**2+mt**2*twidth**2)
     &     *((sw1-wmass**2)**2+wmass**2*wwidth**2)
     &     *((sw2-wmass**2)**2+wmass**2*wwidth**2)
     &     *sz**2

      denu=uDu-mt**2
      denb=bDb-mt**2

      propz=sz/cplx2((sz-zmass**2),zmass*zwidth)
      denq1=q1Dq1-mt**2
      denq2=q2Dq2-mt**2
      dena1=a1Da1-mt**2
      dena2=a2Da2-mt**2
      denz1=z1(4)**2-z1(1)**2-z1(2)**2-z1(3)**2
      denz2=z2(4)**2-z2(1)**2-z2(2)**2-z2(3)**2

      facqqb=V/4d0*(4d0*p3Dp5*p6Dp8)
     & *16d0*gsq**2*esq**2*gwsq**4/densq
      facgg=xn*facqqb
      if (kcase==kqqtthz) then
      facqqb=2d0*xn*facqqb
      facgg=2d0*xn*facgg
      endif

      p4Dq=p(4,4)*q(4)-p(4,1)*q(1)-p(4,2)*q(2)-p(4,3)*q(3)
      p7Da=p(7,4)*a(4)-p(7,1)*a(1)-p(7,2)*a(2)-p(7,3)*a(3)

      p1Du=p(1,4)*u(4)-p(1,1)*u(1)-p(1,2)*u(2)-p(1,3)*u(3)
      p2Du=p(2,4)*u(4)-p(2,1)*u(1)-p(2,2)*u(2)-p(2,3)*u(3)

      p1Db=p(1,4)*b(4)-p(1,1)*b(1)-p(1,2)*b(2)-p(1,3)*b(3)
      p2Db=p(2,4)*b(4)-p(2,1)*b(1)-p(2,2)*b(2)-p(2,3)*b(3)

      gamq4=qDq/(2d0*p4Dq)
      gama7=aDa/(2d0*p7Da)

      gamb1=bDb/(2d0*p1Db)
      gamb2=bDb/(2d0*p2Db)
      gamu1=uDu/(2d0*p1Du)
      gamu2=uDu/(2d0*p2Du)

C     now the momenta 3,5,6,8 are no longer needed
C     so set
      do nu=1,4
      p(q4,nu)=q(nu)-gamq4*p(4,nu)
      p(a7,nu)=a(nu)-gama7*p(7,nu)
      p(u1,nu)=u(nu)-gamu1*p(1,nu)
      p(u2,nu)=u(nu)-gamu2*p(2,nu)
      p(b1,nu)=b(nu)-gamb1*p(1,nu)
      p(b2,nu)=b(nu)-gamb2*p(2,nu)
      enddo

      call spinoru(12,p,za,zb)

      call qqbttz(denu,denb,denz1,denz2,propz,wtqqb,wtqbq)
      call ggttz(denu,denb,denq1,denq2,dena1,dena2,propz,wtgg)

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if (j < 0) then
          msq(j,-j)=aveqq*facqqb*wtqbq(-jj(j))
      elseif (j == 0) then
          msq(j,j)=avegg*facgg*wtgg
      elseif (j > 0) then
          msq(j,-j)=aveqq*facqqb*wtqqb(jj(j))
      endif
      enddo

      return
      end
