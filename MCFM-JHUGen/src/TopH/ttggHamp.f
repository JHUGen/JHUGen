      subroutine ttggHamp(q,q1,q2,q3,q4,eta1,eta2,ampAB)
      implicit none
      include 'constants.f'
      include 'masses.f'
      double precision q(mxpart,4),s134mmtsq,s234mmtsq,s134,s234,s34
      double precision s14,s23,s14mmtsq,s23mmtsq
      double complex spstrng1,spstrng2,
     & spstrng3,spstrng4,cdot,
     & e1(4),e2(4),e3(4),e4(4),p1(4),p2(4),p3(4),p4(4),et1(4),et2(4),
     & e1m(4),e2m(4),e3m(4),e4m(4),e1p(4),e2p(4),e3p(4),e4p(4),
     & p134(4),p234(4),p14(4),p23(4),p34(4),
     & ampAB(2,2,2,2)
      integer q1,q2,q3,q4,eta1,eta2,h1,h2,h3,h4
      p1(:)=dcmplx(q(q1,:))
      p2(:)=dcmplx(q(q2,:))
      p3(:)=dcmplx(q(q3,:))
      p4(:)=dcmplx(q(q4,:))
      et1(:)=dcmplx(q(eta1,:))
      et2(:)=dcmplx(q(eta2,:))
      p14(:)=dcmplx(q(q1,:)+q(q4,:))
      p23(:)=dcmplx(q(q2,:)+q(q3,:))
      p34(:)=dcmplx(q(q3,:)+q(q4,:))
      p134(:)=dcmplx(q(q1,:)+q(q3,:)+q(q4,:))
      p234(:)=dcmplx(q(q2,:)+q(q3,:)+q(q4,:))
      s14=dble(p14(4)**2-p14(1)**2-p14(2)**2-p14(3)**2)
      s23=dble(p23(4)**2-p23(1)**2-p23(2)**2-p23(3)**2)
      s34=dble(p34(4)**2-p34(1)**2-p34(2)**2-p34(3)**2)
      s134=dble(p134(4)**2-p134(1)**2-p134(2)**2-p134(3)**2)
      s234=dble(p234(4)**2-p234(1)**2-p234(2)**2-p234(3)**2)
      s134mmtsq=s134-mt**2
      s234mmtsq=s234-mt**2
      s23mmtsq=s23-mt**2
      s14mmtsq=s14-mt**2
      h1=1
      call VKlSt(p1,mt,et1,h1,e1p)
c      call UbKlSt(p2,mt,et2,h1,e2p)
      call pol_real(p3,h1,e3p)
      call pol_real(p4,h1,e4p)
      h1=-1
      call VklSt(p1,mt,et1,h1,e1m)
      call UbKlst(p2,mt,et2,h1,e2m)
      call pol_real(p3,h1,e3m)
      call pol_real(p4,h1,e4m)
 
      do h1=1,2
      if (h1 .eq. 1) e1(:)=e1m(:)
      if (h1 .eq. 2) e1(:)=e1p(:)
      do h2=1,1
      if (h2 .eq. 1) e2(:)=e2m(:)
c      if (h2 .eq. 2) e2(:)=e2p(:)
      do h3=1,2
      if (h3 .eq. 1) e3(:)=e3m(:)
      if (h3 .eq. 2) e3(:)=e3p(:)
      do h4=1,2
      if (h4 .eq. 1) e4(:)=e4m(:)
      if (h4 .eq. 2) e4(:)=e4p(:)
      ampAB(h1,h2,h3,h4)= + s34**(-1)*s234mmtsq**(-1) * ( 2.D0*cdot(p3,
     &    e4)*spstrng2(e2,e3,p234,e1) - 2.D0*cdot(p4,e3)*spstrng2(e2,e4
     &    ,p234,e1) + 2.D0*cdot(e3,e4)*spstrng2(e2,p4,p234,e1) - cdot(
     &    e3,e4)*spstrng2(e2,p34,p234,e1) )
      ampAB(h1,h2,h3,h4) = ampAB(h1,h2,h3,h4) + s34**(-1)*
     & s234mmtsq**(-1)*mt * ( 2.D0*cdot(p3,e4)*spstrng1(e2,e3,e1) - 2.D0
     &    *cdot(p4,e3)*spstrng1(e2,e4,e1) + 2.D0*cdot(e3,e4)*spstrng1(
     &    e2,p4,e1) - cdot(e3,e4)*spstrng1(e2,p34,e1) )
      ampAB(h1,h2,h3,h4) = ampAB(h1,h2,h3,h4) + s34**(-1)*
     & s134mmtsq**(-1) * (  - 2.D0*cdot(p3,e4)*spstrng2(e2,p134,e3,e1)
     &     + 2.D0*cdot(p4,e3)*spstrng2(e2,p134,e4,e1) - 2.D0*cdot(e3,e4
     &    )*spstrng2(e2,p134,p4,e1) + cdot(e3,e4)*spstrng2(e2,p134,p34,
     &    e1) )
      ampAB(h1,h2,h3,h4) = ampAB(h1,h2,h3,h4) + s34**(-1)*
     & s134mmtsq**(-1)*mt * ( 2.D0*cdot(p3,e4)*spstrng1(e2,e3,e1) - 2.D0
     &    *cdot(p4,e3)*spstrng1(e2,e4,e1) + 2.D0*cdot(e3,e4)*spstrng1(
     &    e2,p4,e1) - cdot(e3,e4)*spstrng1(e2,p34,e1) )
      ampAB(h1,h2,h3,h4) = ampAB(h1,h2,h3,h4) + s234mmtsq**(-1) * ( 
     &    spstrng4(e2,e3,p23,e4,p234,e1)*s23mmtsq**(-1) )
      ampAB(h1,h2,h3,h4) = ampAB(h1,h2,h3,h4) + s234mmtsq**(-1)*mt * ( 
     &    spstrng3(e2,e3,e4,p234,e1)*s23mmtsq**(-1) + spstrng3(e2,e3,
     &    p23,e4,e1)*s23mmtsq**(-1) )
      ampAB(h1,h2,h3,h4) = ampAB(h1,h2,h3,h4) + s234mmtsq**(-1)*mt**2
     &  * ( spstrng2(e2,e3,e4,e1)*s23mmtsq**(-1) )
      ampAB(h1,h2,h3,h4) = ampAB(h1,h2,h3,h4) + s134mmtsq**(-1) * ( 
     &    spstrng4(e2,p134,e3,p14,e4,e1)*s14mmtsq**(-1) )
      ampAB(h1,h2,h3,h4) = ampAB(h1,h2,h3,h4) + s134mmtsq**(-1)*mt * ( 
     &     - spstrng3(e2,e3,p14,e4,e1)*s14mmtsq**(-1) - spstrng3(e2,
     &    p134,e3,e4,e1)*s14mmtsq**(-1) )
      ampAB(h1,h2,h3,h4) = ampAB(h1,h2,h3,h4) + s134mmtsq**(-1)*mt**2
     &  * ( spstrng2(e2,e3,e4,e1)*s14mmtsq**(-1) )
      ampAB(h1,h2,h3,h4) = ampAB(h1,h2,h3,h4) + mt * (  - spstrng3(e2,
     &    e3,p14,e4,e1)*s14mmtsq**(-1)*s23mmtsq**(-1) + spstrng3(e2,e3,
     &    p23,e4,e1)*s14mmtsq**(-1)*s23mmtsq**(-1) )
      ampAB(h1,h2,h3,h4) = ampAB(h1,h2,h3,h4) + mt**2 * ( spstrng2(e2,
     &    e3,e4,e1)*s14mmtsq**(-1)*s23mmtsq**(-1) )
      ampAB(h1,h2,h3,h4) = ampAB(h1,h2,h3,h4) - spstrng4(e2,e3,p23,p14,
     & e4,e1)*s14mmtsq**(-1)*s23mmtsq**(-1)

      enddo
      enddo
      enddo
      enddo

      do h1=1,2
      do h3=1,2
      do h4=1,2
!      h2=2
      ampAB(h1,2,h3,h4)=-(-1)**h1*dconjg(ampAB(3-h1,1,h3,h4))
      enddo
      enddo
      enddo
      return
      end
