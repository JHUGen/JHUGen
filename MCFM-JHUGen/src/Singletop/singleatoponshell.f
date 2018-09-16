      subroutine singleatoponshell(q1,q2,q7,p,iswitch,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*                                                                      *
*     q(-p1)+g(-p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+b(p6)+q'(p7)     *
*                                                                      *
*     keeping polarization information for t~,pb and gluon             *
*     iswitch= 0 for no gluon emission                                 *
*     iswitch=-1 for gluon emission in atop decay                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      real(dp):: p(mxpart,4),q(mxpart,4),
     & ala,alb,dot,s17,mt2,twopaDp2,twopbDp2
      complex(dp):: m(2,2,2),iza,izb,cprop
      integer:: p1,p2,a,p7,b,e,si,aa,bb,iswitch,q1,q2,q7
      parameter(p1=1,p2=2,a=3,e=4,b=5,p7=6)
C-----matrix element for p1+p2(g) -> t+c where both t and c are on shell
C-----t rendered massless wrt eb, and c rendered massless wrt p2
c--- statement functions
      iza(aa,bb)=cone/za(aa,bb)
      izb(aa,bb)=cone/zb(aa,bb)
c--- end statement functions
c      write(6,*) 'iswitch',iswitch
c      call writeout(p)
c      pause
C---zero all arrays
      m(:,:,:)=czip
      do si=1,4
      q(p1,si)=p(q1,si)
      q(p2,si)=p(q2,si)
      if (iswitch == 0) then
      q(a,si)=p(3,si)+p(4,si)+p(5,si)
      elseif (iswitch == -1) then
      q(a,si)=p(3,si)+p(4,si)+p(5,si)+p(8,si)
      endif
      q(e,si)=p(3,si)
      q(b,si)=p(6,si)
      q(p7,si)=p(q7,si)
      enddo
      mt2=mt**2
C---- now render "a" massless wrt to vector e
C---- now render "b" massless wrt to vector p2
      ala=mt2/(2._dp*dot(q,a,e))
      alb=mb**2/(2._dp*dot(q,b,p2))
      s17=2._dp*dot(q,p1,p7)
      cprop=cplx2(s17-wmass**2,zip)
      twopaDp2=2._dp*dot(q,a,p2)
      twopbDp2=2._dp*dot(q,b,p2)
      do si=1,4
      q(a,si)=q(a,si)-ala*q(e,si)
      q(b,si)=q(b,si)-alb*q(p2,si)
      enddo
      call spinoru(6,q,za,zb)
      
C----order of indices is polb,polg,polatop
      m(1,1,2)= + cprop**(-1)*twopbDp2**(-1) * (  - za(b,p2)*za(p2,p7)*
     &    zb(a,p1) - za(b,p2)*za(p2,p7)*zb(a,p1)*alb )
      m(1,1,2) = m(1,1,2) + cprop**(-1)*twopbDp2**(-1)*mb**2 * (  - za(
     &    p2,p7)*zb(a,p1)*izb(b,p2) )
      m(1,1,2) = m(1,1,2) + cprop**(-1)*twopaDp2**(-1) * ( za(e,p2)*za(
     &    b,p7)*zb(e,p1)*zb(b,a)*izb(b,p2)*ala + za(b,p7)*za(a,p2)*zb(b
     &    ,a)*zb(a,p1)*izb(b,p2) )
      m(1,1,2) = m(1,1,2) + cprop**(-1)*mt2*twopaDp2**(-1) * ( za(e,p2)
     &    *za(b,p7)*zb(b,p1)*iza(e,a)*izb(b,p2) )

      m(1,2,2)= + cprop**(-1)*twopaDp2**(-1) * ( za(e,b)*za(b,p7)*zb(e,
     &    p1)*zb(a,p2)*iza(b,p2)*ala - za(b,a)*za(b,p7)*zb(a,p1)*zb(a,
     &    p2)*iza(b,p2) + za(b,p7)*zb(a,p2)*zb(p1,p2) )
      m(1,2,2) = m(1,2,2) + cprop**(-1)*mt2*twopaDp2**(-1) * ( za(e,b)*
     &    za(b,p7)*zb(p1,p2)*iza(e,a)*iza(b,p2) )

      m(1,1,1)= + cprop**(-1)*mt*twopbDp2**(-1) * (  - za(b,p2)*za(p2,
     &    p7)*zb(e,p1)*izb(e,a) - za(b,p2)*za(p2,p7)*zb(e,p1)*izb(e,a)*
     &    alb )
      m(1,1,1) = m(1,1,1) + cprop**(-1)*mt*twopbDp2**(-1)*mb**2 * (  - 
     &    za(p2,p7)*zb(e,p1)*izb(e,a)*izb(b,p2) )
      m(1,1,1) = m(1,1,1) + cprop**(-1)*mt*twopaDp2**(-1) * (  - za(e,
     &    p2)*za(b,p7)*zb(e,b)*zb(e,p1)*izb(e,a)*izb(b,p2)*ala - za(b,
     &    p7)*za(a,p2)*zb(e,b)*zb(a,p1)*izb(e,a)*izb(b,p2) + za(b,p7)*
     &    za(a,p2)*zb(b,p1)*izb(b,p2) )

      m(1,2,1)= + cprop**(-1)*mt*twopaDp2**(-1) * ( za(e,b)*za(b,p7)*
     &    zb(e,p1)*zb(e,p2)*iza(b,p2)*izb(e,a)*ala - za(b,a)*za(b,p7)*
     &    zb(e,p2)*zb(a,p1)*iza(b,p2)*izb(e,a) - za(b,a)*za(b,p7)*zb(p1
     &    ,p2)*iza(b,p2) + za(b,p7)*zb(e,p2)*zb(p1,p2)*izb(e,a) )

      
      m(2,1,2)= + cprop**(-1)*twopaDp2**(-1)*mb * (  - za(e,p2)*za(p2,
     &    p7)*zb(e,p1)*zb(b,a)*iza(b,p2)*izb(b,p2)*ala - za(a,p2)*za(p2
     &    ,p7)*zb(b,a)*zb(a,p1)*iza(b,p2)*izb(b,p2) )
      m(2,1,2) = m(2,1,2) + cprop**(-1)*mt2*twopaDp2**(-1)*mb * (  - 
     &    za(e,p2)*za(p2,p7)*zb(b,p1)*iza(e,a)*iza(b,p2)*izb(b,p2) )

      m(2,2,2)= + cprop**(-1)*twopaDp2**(-1)*mb * (  - za(e,b)*za(p2,p7
     &    )*zb(e,p1)*zb(a,p2)*iza(b,p2)**2*ala + za(b,a)*za(p2,p7)*zb(a
     &    ,p1)*zb(a,p2)*iza(b,p2)**2 - za(p2,p7)*zb(a,p2)*zb(p1,p2)*
     &    iza(b,p2) )
      m(2,2,2) = m(2,2,2) + cprop**(-1)*mt2*twopaDp2**(-1)*mb * (  - 
     &    za(e,b)*za(p2,p7)*zb(p1,p2)*iza(e,a)*iza(b,p2)**2 )

      m(2,1,1)= + cprop**(-1)*mt*twopaDp2**(-1)*mb * ( za(e,p2)*za(p2,
     &    p7)*zb(e,b)*zb(e,p1)*iza(b,p2)*izb(e,a)*izb(b,p2)*ala + za(a,
     &    p2)*za(p2,p7)*zb(e,b)*zb(a,p1)*iza(b,p2)*izb(e,a)*izb(b,p2)
     &     - za(a,p2)*za(p2,p7)*zb(b,p1)*iza(b,p2)*izb(b,p2) )

      m(2,2,1)= + cprop**(-1)*mt*twopaDp2**(-1)*mb * (  - za(e,b)*za(p2
     &    ,p7)*zb(e,p1)*zb(e,p2)*iza(b,p2)**2*izb(e,a)*ala + za(b,a)*
     &    za(p2,p7)*zb(e,p2)*zb(a,p1)*iza(b,p2)**2*izb(e,a) + za(b,a)*
     &    za(p2,p7)*zb(p1,p2)*iza(b,p2)**2 - za(p2,p7)*zb(e,p2)*zb(p1,
     &    p2)*iza(b,p2)*izb(e,a) )

      return
      end
