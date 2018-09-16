      subroutine Watoponshell(q1,q2,p,iswitch,m)
      implicit none
      include 'types.f'
      
************************************************************************
*     Author: R.K. Ellis, January 2012                                 *
*                                                                      *
*     q(-p1)+g(-p2) --> t~(-->e^-(p3)+nu~(p4)+b~(p5))+b(p6)+q'(p3)     *
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
     & ala,dot,s34,mt2,twopaDp2
      complex(dp):: m(2,2,2),iza,izb,cprop
      integer:: p2,a,p3,p4,b,e,si,aa,bb,iswitch,q1,q2
      parameter(b=1,p2=2,p3=3,p4=4,a=5,e=6)
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
      q(b,si)=p(q1,si)
      q(p2,si)=p(q2,si)
      q(p3,si)=p(3,si)
      q(p4,si)=p(4,si)
      if (iswitch == 0) then
      q(a,si)=p(5,si)+p(6,si)+p(7,si)
      elseif (iswitch == -1) then
      q(a,si)=p(5,si)+p(6,si)+p(7,si)+p(8,si)
      endif
      q(e,si)=p(5,si)
      enddo
      mt2=mt**2
C---- now render "a" massless wrt to vector e
      ala=mt2/(two*dot(q,a,e))
      s34=two*dot(q,p3,p4)
      cprop=cplx2(s34-wmass**2,wmass*wwidth)
      twopaDp2=two*dot(q,a,p2)
      do si=1,4
      q(a,si)=q(a,si)-ala*q(e,si)
      enddo
      call spinoru(6,q,za,zb)
      
C----order of indices is polb,polg,polatop
      m(1,1,2)= + cprop**(-1)*twopaDp2**(-1) * ( za(e,p2)*za(b,p3)*zb(e
     &    ,p4)*zb(b,a)*izb(b,p2)*ala + za(b,p3)*za(a,p2)*zb(b,a)*zb(a,
     &    p4)*izb(b,p2) )
      m(1,1,2) = m(1,1,2) + cprop**(-1) * ( za(p2,p3)*zb(a,p4)*izb(b,p2
     &    ) )
      m(1,1,2) = m(1,1,2) + cprop**(-1)*mt2*twopaDp2**(-1) * ( za(e,p2)
     &    *za(b,p3)*zb(b,p4)*iza(e,a)*izb(b,p2) )

      m(1,2,2)= + cprop**(-1)*twopaDp2**(-1) * ( za(e,b)*za(b,p3)*zb(e,
     &    p4)*zb(a,p2)*iza(b,p2)*ala - za(b,a)*za(b,p3)*zb(a,p2)*zb(a,
     &    p4)*iza(b,p2) - za(b,p3)*zb(a,p2)*zb(p2,p4) )
      m(1,2,2) = m(1,2,2) + cprop**(-1)*mt2*twopaDp2**(-1) * (  - za(e,
     &    b)*za(b,p3)*zb(p2,p4)*iza(e,a)*iza(b,p2) )

      m(1,1,1)= + cprop**(-1)*mt*twopaDp2**(-1) * (  - za(e,p2)*za(b,p3
     &    )*zb(e,b)*zb(e,p4)*izb(e,a)*izb(b,p2)*ala - za(b,p3)*za(a,p2)
     &    *zb(e,b)*zb(a,p4)*izb(e,a)*izb(b,p2) + za(b,p3)*za(a,p2)*zb(b
     &    ,p4)*izb(b,p2) )
      m(1,1,1) = m(1,1,1) + cprop**(-1)*mt * ( za(p2,p3)*zb(e,p4)*izb(e
     &    ,a)*izb(b,p2) )

      m(1,2,1)= + cprop**(-1)*mt*twopaDp2**(-1) * ( za(e,b)*za(b,p3)*
     &    zb(e,p2)*zb(e,p4)*iza(b,p2)*izb(e,a)*ala - za(b,a)*za(b,p3)*
     &    zb(e,p2)*zb(a,p4)*iza(b,p2)*izb(e,a) + za(b,a)*za(b,p3)*zb(p2
     &    ,p4)*iza(b,p2) - za(b,p3)*zb(e,p2)*zb(p2,p4)*izb(e,a) )

      
      m(2,1,2)=czip

      m(2,2,2)=czip

      m(2,1,1)=czip

      m(2,2,1)=czip

      return
      end
