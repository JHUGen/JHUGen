      subroutine qqbZtt(pin,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     Jan, 2011.                                                       *
*     Calculate the lowest order element squared for the process       *
*                                                                      *
*     q(-p1) +qbar(-p2) -> tbar(bbar(p6)+e-(p7)+nubar(p8))             *
*                         +t(b(p5)+nu(p3)+e+(p4))                      *
*     mediated by s-channel Z-boson exchange                           *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ewcouple.f'
      include 'zprods_com.f'
      include 'masses.f'
      include 'topzlabels.f'
      include 'cplx.h'

      integer:: j,k,nu
      real(dp):: msq(-nf:nf,-nf:nf),pin(mxpart,4),p(mxpart,4)
      real(dp):: pw1(4),pw2(4),p12(4),q(4),a(4)
      real(dp):: s12,sw1,sw2,qDq,aDa,densq,p3Dp5,p6Dp8
      real(dp):: wtqqb(2),wtqbq(2)
      real(dp):: gamq4,gama7,dot
      real(dp):: facqqb,p4Dq,p7Da
      complex(dp):: propz
      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      do nu=1,4
      do j=1,8
      p(j,nu)=pin(j,nu)
      enddo

      p12(nu)=p(1,nu)+p(2,nu)
      pw1(nu)=p(3,nu)+p(4,nu)
      pw2(nu)=p(7,nu)+p(8,nu)

      q(nu)=+p(3,nu)+p(4,nu)+p(5,nu)
      a(nu)=-p(6,nu)-p(7,nu)-p(8,nu)
      enddo


      s12=(p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2)

      if (s12 < 4._dp*mt**2) then
      write(6,*) 'qqbZtt: s12',s12
      return
      endif

      sw1=(pw1(4)**2-pw1(1)**2-pw1(2)**2-pw1(3)**2)
      sw2=(pw2(4)**2-pw2(1)**2-pw2(2)**2-pw2(3)**2)
      qDq=(q(4)**2-q(1)**2-q(2)**2-q(3)**2)
      aDa=(a(4)**2-a(1)**2-a(2)**2-a(3)**2)

      p3Dp5=dot(p,3,5)
      p6Dp8=dot(p,6,8)

      densq=((qDq-mt**2)**2+mt**2*twidth**2)
     &     *((aDa-mt**2)**2+mt**2*twidth**2)
     &     *((sw1-wmass**2)**2+wmass**2*wwidth**2)
     &     *((sw2-wmass**2)**2+wmass**2*wwidth**2)
     &     *s12**2

      propz=s12/(cplx2((s12-zmass**2),zmass*zwidth))

      facqqb=xnsq*(4._dp*p3Dp5*p6Dp8)
     & *4._dp*esq**2*gwsq**4/densq

      p4Dq=p(4,4)*q(4)-p(4,1)*q(1)-p(4,2)*q(2)-p(4,3)*q(3)
      p7Da=p(7,4)*a(4)-p(7,1)*a(1)-p(7,2)*a(2)-p(7,3)*a(3)

      gamq4=qDq/(2._dp*p4Dq)
      gama7=aDa/(2._dp*p7Da)


      do nu=1,4
      p(q4,nu)=q(nu)-gamq4*p(4,nu)
      p(a7,nu)=a(nu)-gama7*p(7,nu)
      enddo

      call spinoru(7,p,za,zb)
      call qqbZtt1(propz,wtqqb,wtqbq)


C---fill qb-q and q-qb elements
      do j=-nf,nf
      if (j < 0) then
          msq(j,-j)=aveqq*facqqb*wtqbq(-jj(j))
      elseif (j > 0) then
          msq(j,-j)=aveqq*facqqb*wtqqb(jj(j))
      endif
      enddo

      return
      end
