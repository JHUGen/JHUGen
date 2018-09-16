      subroutine ttgggHamp(q,q1,q2,q3,q4,q5,eta1,eta2,ampABC)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      real(dp):: q(mxpart,4),s234mmtsq,
     & s234,s34,s45,s345,s15,s23,s15mmtsq,s23mmtsq,s1345mmtsq,
     & s2345mmtsq,s145mmtsq,s145,s1345,s2345
      complex(dp):: spstrng1,spstrng2,
     & spstrng3,spstrng4,spstrng5,spstrng6,cdot,
     & e1(4),e2(4),e3(4),e4(4),e5(4),
     & p1(4),p2(4),p3(4),p4(4),p5(4),et1(4),et2(4),
     & e1m(4),e2m(4),e3m(4),e4m(4),e5m(4),
     & e1p(4),e2p(4),e3p(4),e4p(4),e5p(4),
     & p234(4),p23(4),p34(4),p45(4),p345(4),p1345(4),
     & p2345(4),p15(4),p145(4),ampABC(2,2,2,2,2)
      integer:: q1,q2,q3,q4,q5,eta1,eta2,h1,h2,h3,h4,h5
      p1(:)=cplx2(q(q1,:))
      p2(:)=cplx2(q(q2,:))
      p3(:)=cplx2(q(q3,:))
      p4(:)=cplx2(q(q4,:))
      p5(:)=cplx2(q(q5,:))
      et1(:)=cplx2(q(eta1,:))
      et2(:)=cplx2(q(eta2,:))
      p15(:)=cplx2(q(q1,:)+q(q5,:))
      p23(:)=cplx2(q(q2,:)+q(q3,:))
      p34(:)=cplx2(q(q3,:)+q(q4,:))
      p45(:)=cplx2(q(q4,:)+q(q5,:))
      p145(:)=cplx2(q(q1,:)+q(q4,:)+q(q5,:))
      p345(:)=cplx2(q(q3,:)+q(q4,:)+q(q5,:))
      p1345(:)=cplx2(q(q1,:)+q(q3,:)+q(q4,:)+q(q5,:))
      p2345(:)=cplx2(q(q2,:)+q(q3,:)+q(q4,:)+q(q5,:))
      p234(:)=cplx2(q(q2,:)+q(q3,:)+q(q4,:))
      s45=real(p45(4)**2-p45(1)**2-p45(2)**2-p45(3)**2)
      s15=real(p15(4)**2-p15(1)**2-p15(2)**2-p15(3)**2)
      s23=real(p23(4)**2-p23(1)**2-p23(2)**2-p23(3)**2)
      s34=real(p34(4)**2-p34(1)**2-p34(2)**2-p34(3)**2)
      s145=real(p145(4)**2-p145(1)**2-p145(2)**2-p145(3)**2)
      s1345=real(p1345(4)**2-p1345(1)**2-p1345(2)**2-p1345(3)**2)
      s2345=real(p2345(4)**2-p2345(1)**2-p2345(2)**2-p2345(3)**2)
      s345=real(p345(4)**2-p345(1)**2-p345(2)**2-p345(3)**2)
      s234=real(p234(4)**2-p234(1)**2-p234(2)**2-p234(3)**2)
      s1345mmtsq=s1345-mt**2
      s2345mmtsq=s2345-mt**2
      s145mmtsq=s145-mt**2
      s234mmtsq=s234-mt**2
      s23mmtsq=s23-mt**2
      s15mmtsq=s15-mt**2
      h1=1
      call VKlSt(p1,mt,et1,h1,e1p)
c      call UbKlSt(p2,mt,et2,h1,e2p)
      call pol_real(p3,h1,e3p)
      call pol_real(p4,h1,e4p)
      call pol_real(p5,h1,e5p)
      h1=-1
      call VklSt(p1,mt,et1,h1,e1m)
      call UbKlst(p2,mt,et2,h1,e2m)
      call pol_real(p3,h1,e3m)
      call pol_real(p4,h1,e4m)
      call pol_real(p5,h1,e5m)
 
      do h1=1,2
      if (h1 == 1) e1(:)=e1m(:)
      if (h1 == 2) e1(:)=e1p(:)
      do h2=1,1
      if (h2 == 1) e2(:)=e2m(:)
c      if (h2 == 2) e2(:)=e2p(:)
      do h3=1,2
      if (h3 == 1) e3(:)=e3m(:)
      if (h3 == 2) e3(:)=e3p(:)
      do h4=1,2
      if (h4 == 1) e4(:)=e4m(:)
      if (h4 == 2) e4(:)=e4p(:)
      do h5=1,2
      if (h5 == 1) e5(:)=e5m(:)
      if (h5 == 2) e5(:)=e5p(:)
      ampABC(h1,h2,h3,h4,h5)= + s1345mmtsq**(-1) * ( 2.D0*cdot(p3,p4)*
     &    cdot(e4,e5)*spstrng2(e2,p1345,e3,e1)*s345**(-1)*s45**(-1) - 2.
     &    D0*cdot(p3,p5)*cdot(e3,e4)*spstrng2(e2,p1345,e5,e1)*s34**(-1)
     &    *s345**(-1) - 2.D0*cdot(p3,p5)*cdot(e4,e5)*spstrng2(e2,p1345,
     &    e3,e1)*s345**(-1)*s45**(-1) - 4.D0*cdot(p3,e4)*cdot(p3,e5)*
     &    spstrng2(e2,p1345,e3,e1)*s34**(-1)*s345**(-1) - 4.D0*cdot(p3,
     &    e4)*cdot(p4,e5)*spstrng2(e2,p1345,e3,e1)*s34**(-1)*s345**(-1)
     &     - 4.D0*cdot(p3,e4)*cdot(p4,e5)*spstrng2(e2,p1345,e3,e1)*
     &    s345**(-1)*s45**(-1) + 4.D0*cdot(p3,e4)*cdot(p5,e3)*spstrng2(
     &    e2,p1345,e5,e1)*s34**(-1)*s345**(-1) + 2.D0*cdot(p3,e4)*cdot(
     &    e3,e5)*spstrng2(e2,p1345,p3,e1)*s34**(-1)*s345**(-1) + 2.D0*
     &    cdot(p3,e4)*cdot(e3,e5)*spstrng2(e2,p1345,p4,e1)*s34**(-1)*
     &    s345**(-1) - 2.D0*cdot(p3,e4)*cdot(e3,e5)*spstrng2(e2,p1345,
     &    p5,e1)*s34**(-1)*s345**(-1) + 2.D0*cdot(p3,e4)*spstrng4(e2,
     &    p1345,e3,p15,e5,e1)*s34**(-1)*s15mmtsq**(-1) + 4.D0*cdot(p3,
     &    e5)*cdot(p4,e3)*spstrng2(e2,p1345,e4,e1)*s34**(-1)*s345**(-1)
     &     )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + 
     & s1345mmtsq**(-1) * ( 4.D0*cdot(p3,e5)*cdot(p5,e4)*spstrng2(e2,
     &    p1345,e3,e1)*s345**(-1)*s45**(-1) + cdot(p3,e5)*cdot(e3,e4)*
     &    spstrng2(e2,p1345,p3,e1)*s34**(-1)*s345**(-1) - 3.D0*cdot(p3,
     &    e5)*cdot(e3,e4)*spstrng2(e2,p1345,p4,e1)*s34**(-1)*s345**(-1)
     &     + cdot(p3,e5)*cdot(e3,e4)*spstrng2(e2,p1345,p5,e1)*s34**(-1)
     &    *s345**(-1) + 2.D0*cdot(p4,p5)*cdot(e3,e4)*spstrng2(e2,p1345,
     &    e5,e1)*s34**(-1)*s345**(-1) + 4.D0*cdot(p4,e3)*cdot(p4,e5)*
     &    spstrng2(e2,p1345,e4,e1)*s34**(-1)*s345**(-1) + 4.D0*cdot(p4,
     &    e3)*cdot(p4,e5)*spstrng2(e2,p1345,e4,e1)*s345**(-1)*s45**(-1)
     &     - 4.D0*cdot(p4,e3)*cdot(p5,e4)*spstrng2(e2,p1345,e5,e1)*
     &    s34**(-1)*s345**(-1) - 4.D0*cdot(p4,e3)*cdot(p5,e4)*spstrng2(
     &    e2,p1345,e5,e1)*s345**(-1)*s45**(-1) - 2.D0*cdot(p4,e3)*cdot(
     &    e4,e5)*spstrng2(e2,p1345,p3,e1)*s34**(-1)*s345**(-1) - cdot(
     &    p4,e3)*cdot(e4,e5)*spstrng2(e2,p1345,p3,e1)*s345**(-1)*
     &    s45**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + 
     & s1345mmtsq**(-1) * (  - 2.D0*cdot(p4,e3)*cdot(e4,e5)*spstrng2(e2
     &    ,p1345,p4,e1)*s34**(-1)*s345**(-1) - cdot(p4,e3)*cdot(e4,e5)*
     &    spstrng2(e2,p1345,p4,e1)*s345**(-1)*s45**(-1) + 2.D0*cdot(p4,
     &    e3)*cdot(e4,e5)*spstrng2(e2,p1345,p5,e1)*s34**(-1)*s345**(-1)
     &     + 3.D0*cdot(p4,e3)*cdot(e4,e5)*spstrng2(e2,p1345,p5,e1)*
     &    s345**(-1)*s45**(-1) - 2.D0*cdot(p4,e3)*spstrng4(e2,p1345,e4,
     &    p15,e5,e1)*s34**(-1)*s15mmtsq**(-1) + 4.D0*cdot(p4,e5)*cdot(
     &    p5,e3)*spstrng2(e2,p1345,e4,e1)*s345**(-1)*s45**(-1) + 3.D0*
     &    cdot(p4,e5)*cdot(e3,e4)*spstrng2(e2,p1345,p3,e1)*s34**(-1)*
     &    s345**(-1) + 2.D0*cdot(p4,e5)*cdot(e3,e4)*spstrng2(e2,p1345,
     &    p3,e1)*s345**(-1)*s45**(-1) - cdot(p4,e5)*cdot(e3,e4)*
     &    spstrng2(e2,p1345,p4,e1)*s34**(-1)*s345**(-1) - 2.D0*cdot(p4,
     &    e5)*cdot(e3,e4)*spstrng2(e2,p1345,p4,e1)*s345**(-1)*s45**(-1)
     &     - cdot(p4,e5)*cdot(e3,e4)*spstrng2(e2,p1345,p5,e1)*s34**(-1)
     &    *s345**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + 
     & s1345mmtsq**(-1) * (  - 2.D0*cdot(p4,e5)*cdot(e3,e4)*spstrng2(e2
     &    ,p1345,p5,e1)*s345**(-1)*s45**(-1) + 2.D0*cdot(p4,e5)*
     &    spstrng4(e2,p1345,e3,p145,e4,e1)*s145mmtsq**(-1)*s45**(-1) - 
     &    4.D0*cdot(p5,e3)*cdot(p5,e4)*spstrng2(e2,p1345,e5,e1)*
     &    s345**(-1)*s45**(-1) + cdot(p5,e3)*cdot(e4,e5)*spstrng2(e2,
     &    p1345,p3,e1)*s345**(-1)*s45**(-1) - 3.D0*cdot(p5,e3)*cdot(e4,
     &    e5)*spstrng2(e2,p1345,p4,e1)*s345**(-1)*s45**(-1) + cdot(p5,
     &    e3)*cdot(e4,e5)*spstrng2(e2,p1345,p5,e1)*s345**(-1)*s45**(-1)
     &     - 2.D0*cdot(p5,e4)*cdot(e3,e5)*spstrng2(e2,p1345,p3,e1)*
     &    s345**(-1)*s45**(-1) + 2.D0*cdot(p5,e4)*cdot(e3,e5)*spstrng2(
     &    e2,p1345,p4,e1)*s345**(-1)*s45**(-1) + 2.D0*cdot(p5,e4)*cdot(
     &    e3,e5)*spstrng2(e2,p1345,p5,e1)*s345**(-1)*s45**(-1) - 2.D0*
     &    cdot(p5,e4)*spstrng4(e2,p1345,e3,p145,e5,e1)*s145mmtsq**(-1)*
     &    s45**(-1) + cdot(e3,e4)*spstrng2(e2,p1345,e5,e1)*s345**(-1)
     &     - cdot(e3,e4)*spstrng4(e2,p1345,p3,p15,e5,e1)*s34**(-1)*
     &    s15mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + 
     & s1345mmtsq**(-1) * ( cdot(e3,e4)*spstrng4(e2,p1345,p4,p15,e5,e1)
     &    *s34**(-1)*s15mmtsq**(-1) - 2.D0*cdot(e3,e5)*spstrng2(e2,
     &    p1345,e4,e1)*s345**(-1) + cdot(e4,e5)*spstrng2(e2,p1345,e3,e1
     &    )*s345**(-1) - cdot(e4,e5)*spstrng4(e2,p1345,e3,p145,p4,e1)*
     &    s145mmtsq**(-1)*s45**(-1) + cdot(e4,e5)*spstrng4(e2,p1345,e3,
     &    p145,p5,e1)*s145mmtsq**(-1)*s45**(-1) - spstrng6(e2,p1345,e3,
     &    p145,e4,p15,e5,e1)*s145mmtsq**(-1)*s15mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + 
     & s2345mmtsq**(-1) * (  - 2.D0*cdot(p3,p4)*cdot(e4,e5)*spstrng2(e2
     &    ,e3,p2345,e1)*s345**(-1)*s45**(-1) + 2.D0*cdot(p3,p5)*cdot(e3
     &    ,e4)*spstrng2(e2,e5,p2345,e1)*s34**(-1)*s345**(-1) + 2.D0*
     &    cdot(p3,p5)*cdot(e4,e5)*spstrng2(e2,e3,p2345,e1)*s345**(-1)*
     &    s45**(-1) + 4.D0*cdot(p3,e4)*cdot(p3,e5)*spstrng2(e2,e3,p2345
     &    ,e1)*s34**(-1)*s345**(-1) + 4.D0*cdot(p3,e4)*cdot(p4,e5)*
     &    spstrng2(e2,e3,p2345,e1)*s34**(-1)*s345**(-1) + 4.D0*cdot(p3,
     &    e4)*cdot(p4,e5)*spstrng2(e2,e3,p2345,e1)*s345**(-1)*s45**(-1)
     &     - 4.D0*cdot(p3,e4)*cdot(p5,e3)*spstrng2(e2,e5,p2345,e1)*
     &    s34**(-1)*s345**(-1) - 2.D0*cdot(p3,e4)*cdot(e3,e5)*spstrng2(
     &    e2,p3,p2345,e1)*s34**(-1)*s345**(-1) - 2.D0*cdot(p3,e4)*cdot(
     &    e3,e5)*spstrng2(e2,p4,p2345,e1)*s34**(-1)*s345**(-1) + 2.D0*
     &    cdot(p3,e4)*cdot(e3,e5)*spstrng2(e2,p5,p2345,e1)*s34**(-1)*
     &    s345**(-1) + 2.D0*cdot(p3,e4)*spstrng4(e2,e3,p234,e5,p2345,e1
     &    )*s34**(-1)*s234mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + 
     & s2345mmtsq**(-1) * (  - 4.D0*cdot(p3,e5)*cdot(p4,e3)*spstrng2(e2
     &    ,e4,p2345,e1)*s34**(-1)*s345**(-1) - 4.D0*cdot(p3,e5)*cdot(p5
     &    ,e4)*spstrng2(e2,e3,p2345,e1)*s345**(-1)*s45**(-1) - cdot(p3,
     &    e5)*cdot(e3,e4)*spstrng2(e2,p3,p2345,e1)*s34**(-1)*s345**(-1)
     &     + 3.D0*cdot(p3,e5)*cdot(e3,e4)*spstrng2(e2,p4,p2345,e1)*
     &    s34**(-1)*s345**(-1) - cdot(p3,e5)*cdot(e3,e4)*spstrng2(e2,p5
     &    ,p2345,e1)*s34**(-1)*s345**(-1) - 2.D0*cdot(p4,p5)*cdot(e3,e4
     &    )*spstrng2(e2,e5,p2345,e1)*s34**(-1)*s345**(-1) - 4.D0*cdot(
     &    p4,e3)*cdot(p4,e5)*spstrng2(e2,e4,p2345,e1)*s34**(-1)*
     &    s345**(-1) - 4.D0*cdot(p4,e3)*cdot(p4,e5)*spstrng2(e2,e4,
     &    p2345,e1)*s345**(-1)*s45**(-1) + 4.D0*cdot(p4,e3)*cdot(p5,e4)
     &    *spstrng2(e2,e5,p2345,e1)*s34**(-1)*s345**(-1) + 4.D0*cdot(p4
     &    ,e3)*cdot(p5,e4)*spstrng2(e2,e5,p2345,e1)*s345**(-1)*
     &    s45**(-1) + 2.D0*cdot(p4,e3)*cdot(e4,e5)*spstrng2(e2,p3,p2345
     &    ,e1)*s34**(-1)*s345**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + 
     & s2345mmtsq**(-1) * ( cdot(p4,e3)*cdot(e4,e5)*spstrng2(e2,p3,
     &    p2345,e1)*s345**(-1)*s45**(-1) + 2.D0*cdot(p4,e3)*cdot(e4,e5)
     &    *spstrng2(e2,p4,p2345,e1)*s34**(-1)*s345**(-1) + cdot(p4,e3)*
     &    cdot(e4,e5)*spstrng2(e2,p4,p2345,e1)*s345**(-1)*s45**(-1) - 2.
     &    D0*cdot(p4,e3)*cdot(e4,e5)*spstrng2(e2,p5,p2345,e1)*s34**(-1)
     &    *s345**(-1) - 3.D0*cdot(p4,e3)*cdot(e4,e5)*spstrng2(e2,p5,
     &    p2345,e1)*s345**(-1)*s45**(-1) - 2.D0*cdot(p4,e3)*spstrng4(e2
     &    ,e4,p234,e5,p2345,e1)*s34**(-1)*s234mmtsq**(-1) - 4.D0*cdot(
     &    p4,e5)*cdot(p5,e3)*spstrng2(e2,e4,p2345,e1)*s345**(-1)*
     &    s45**(-1) - 3.D0*cdot(p4,e5)*cdot(e3,e4)*spstrng2(e2,p3,p2345
     &    ,e1)*s34**(-1)*s345**(-1) - 2.D0*cdot(p4,e5)*cdot(e3,e4)*
     &    spstrng2(e2,p3,p2345,e1)*s345**(-1)*s45**(-1) + cdot(p4,e5)*
     &    cdot(e3,e4)*spstrng2(e2,p4,p2345,e1)*s34**(-1)*s345**(-1) + 2.
     &    D0*cdot(p4,e5)*cdot(e3,e4)*spstrng2(e2,p4,p2345,e1)*
     &    s345**(-1)*s45**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + 
     & s2345mmtsq**(-1) * ( cdot(p4,e5)*cdot(e3,e4)*spstrng2(e2,p5,
     &    p2345,e1)*s34**(-1)*s345**(-1) + 2.D0*cdot(p4,e5)*cdot(e3,e4)
     &    *spstrng2(e2,p5,p2345,e1)*s345**(-1)*s45**(-1) + 2.D0*cdot(p4
     &    ,e5)*spstrng4(e2,e3,p23,e4,p2345,e1)*s45**(-1)*s23mmtsq**(-1)
     &     + 4.D0*cdot(p5,e3)*cdot(p5,e4)*spstrng2(e2,e5,p2345,e1)*
     &    s345**(-1)*s45**(-1) - cdot(p5,e3)*cdot(e4,e5)*spstrng2(e2,p3
     &    ,p2345,e1)*s345**(-1)*s45**(-1) + 3.D0*cdot(p5,e3)*cdot(e4,e5
     &    )*spstrng2(e2,p4,p2345,e1)*s345**(-1)*s45**(-1) - cdot(p5,e3)
     &    *cdot(e4,e5)*spstrng2(e2,p5,p2345,e1)*s345**(-1)*s45**(-1) + 
     &    2.D0*cdot(p5,e4)*cdot(e3,e5)*spstrng2(e2,p3,p2345,e1)*
     &    s345**(-1)*s45**(-1) - 2.D0*cdot(p5,e4)*cdot(e3,e5)*spstrng2(
     &    e2,p4,p2345,e1)*s345**(-1)*s45**(-1) - 2.D0*cdot(p5,e4)*cdot(
     &    e3,e5)*spstrng2(e2,p5,p2345,e1)*s345**(-1)*s45**(-1) - 2.D0*
     &    cdot(p5,e4)*spstrng4(e2,e3,p23,e5,p2345,e1)*s45**(-1)*
     &    s23mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + 
     & s2345mmtsq**(-1) * (  - cdot(e3,e4)*spstrng2(e2,e5,p2345,e1)*
     &    s345**(-1) - cdot(e3,e4)*spstrng4(e2,p3,p234,e5,p2345,e1)*
     &    s34**(-1)*s234mmtsq**(-1) + cdot(e3,e4)*spstrng4(e2,p4,p234,
     &    e5,p2345,e1)*s34**(-1)*s234mmtsq**(-1) + 2.D0*cdot(e3,e5)*
     &    spstrng2(e2,e4,p2345,e1)*s345**(-1) - cdot(e4,e5)*spstrng2(e2
     &    ,e3,p2345,e1)*s345**(-1) - cdot(e4,e5)*spstrng4(e2,e3,p23,p4,
     &    p2345,e1)*s45**(-1)*s23mmtsq**(-1) + cdot(e4,e5)*spstrng4(e2,
     &    e3,p23,p5,p2345,e1)*s45**(-1)*s23mmtsq**(-1) + spstrng6(e2,e3
     &    ,p23,e4,p234,e5,p2345,e1)*s234mmtsq**(-1)*s23mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt*
     & s1345mmtsq**(-1) * (  - 2.D0*cdot(p3,p4)*cdot(e4,e5)*spstrng1(e2
     &    ,e3,e1)*s345**(-1)*s45**(-1) + 2.D0*cdot(p3,p5)*cdot(e3,e4)*
     &    spstrng1(e2,e5,e1)*s34**(-1)*s345**(-1) + 2.D0*cdot(p3,p5)*
     &    cdot(e4,e5)*spstrng1(e2,e3,e1)*s345**(-1)*s45**(-1) + 4.D0*
     &    cdot(p3,e4)*cdot(p3,e5)*spstrng1(e2,e3,e1)*s34**(-1)*
     &    s345**(-1) + 4.D0*cdot(p3,e4)*cdot(p4,e5)*spstrng1(e2,e3,e1)*
     &    s34**(-1)*s345**(-1) + 4.D0*cdot(p3,e4)*cdot(p4,e5)*spstrng1(
     &    e2,e3,e1)*s345**(-1)*s45**(-1) - 4.D0*cdot(p3,e4)*cdot(p5,e3)
     &    *spstrng1(e2,e5,e1)*s34**(-1)*s345**(-1) - 2.D0*cdot(p3,e4)*
     &    cdot(e3,e5)*spstrng1(e2,p3,e1)*s34**(-1)*s345**(-1) - 2.D0*
     &    cdot(p3,e4)*cdot(e3,e5)*spstrng1(e2,p4,e1)*s34**(-1)*
     &    s345**(-1) + 2.D0*cdot(p3,e4)*cdot(e3,e5)*spstrng1(e2,p5,e1)*
     &    s34**(-1)*s345**(-1) - 2.D0*cdot(p3,e4)*spstrng3(e2,e3,p15,e5
     &    ,e1)*s34**(-1)*s15mmtsq**(-1) - 2.D0*cdot(p3,e4)*spstrng3(e2,
     &    p1345,e3,e5,e1)*s34**(-1)*s15mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt*
     & s1345mmtsq**(-1) * (  - 4.D0*cdot(p3,e5)*cdot(p4,e3)*spstrng1(e2
     &    ,e4,e1)*s34**(-1)*s345**(-1) - 4.D0*cdot(p3,e5)*cdot(p5,e4)*
     &    spstrng1(e2,e3,e1)*s345**(-1)*s45**(-1) - cdot(p3,e5)*cdot(e3
     &    ,e4)*spstrng1(e2,p3,e1)*s34**(-1)*s345**(-1) + 3.D0*cdot(p3,
     &    e5)*cdot(e3,e4)*spstrng1(e2,p4,e1)*s34**(-1)*s345**(-1) - 
     &    cdot(p3,e5)*cdot(e3,e4)*spstrng1(e2,p5,e1)*s34**(-1)*
     &    s345**(-1) - 2.D0*cdot(p4,p5)*cdot(e3,e4)*spstrng1(e2,e5,e1)*
     &    s34**(-1)*s345**(-1) - 4.D0*cdot(p4,e3)*cdot(p4,e5)*spstrng1(
     &    e2,e4,e1)*s34**(-1)*s345**(-1) - 4.D0*cdot(p4,e3)*cdot(p4,e5)
     &    *spstrng1(e2,e4,e1)*s345**(-1)*s45**(-1) + 4.D0*cdot(p4,e3)*
     &    cdot(p5,e4)*spstrng1(e2,e5,e1)*s34**(-1)*s345**(-1) + 4.D0*
     &    cdot(p4,e3)*cdot(p5,e4)*spstrng1(e2,e5,e1)*s345**(-1)*
     &    s45**(-1) + 2.D0*cdot(p4,e3)*cdot(e4,e5)*spstrng1(e2,p3,e1)*
     &    s34**(-1)*s345**(-1) + cdot(p4,e3)*cdot(e4,e5)*spstrng1(e2,p3
     &    ,e1)*s345**(-1)*s45**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt*
     & s1345mmtsq**(-1) * ( 2.D0*cdot(p4,e3)*cdot(e4,e5)*spstrng1(e2,p4
     &    ,e1)*s34**(-1)*s345**(-1) + cdot(p4,e3)*cdot(e4,e5)*spstrng1(
     &    e2,p4,e1)*s345**(-1)*s45**(-1) - 2.D0*cdot(p4,e3)*cdot(e4,e5)
     &    *spstrng1(e2,p5,e1)*s34**(-1)*s345**(-1) - 3.D0*cdot(p4,e3)*
     &    cdot(e4,e5)*spstrng1(e2,p5,e1)*s345**(-1)*s45**(-1) + 2.D0*
     &    cdot(p4,e3)*spstrng3(e2,e4,p15,e5,e1)*s34**(-1)*
     &    s15mmtsq**(-1) + 2.D0*cdot(p4,e3)*spstrng3(e2,p1345,e4,e5,e1)
     &    *s34**(-1)*s15mmtsq**(-1) - 4.D0*cdot(p4,e5)*cdot(p5,e3)*
     &    spstrng1(e2,e4,e1)*s345**(-1)*s45**(-1) - 3.D0*cdot(p4,e5)*
     &    cdot(e3,e4)*spstrng1(e2,p3,e1)*s34**(-1)*s345**(-1) - 2.D0*
     &    cdot(p4,e5)*cdot(e3,e4)*spstrng1(e2,p3,e1)*s345**(-1)*
     &    s45**(-1) + cdot(p4,e5)*cdot(e3,e4)*spstrng1(e2,p4,e1)*
     &    s34**(-1)*s345**(-1) + 2.D0*cdot(p4,e5)*cdot(e3,e4)*spstrng1(
     &    e2,p4,e1)*s345**(-1)*s45**(-1) + cdot(p4,e5)*cdot(e3,e4)*
     &    spstrng1(e2,p5,e1)*s34**(-1)*s345**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt*
     & s1345mmtsq**(-1) * ( 2.D0*cdot(p4,e5)*cdot(e3,e4)*spstrng1(e2,p5
     &    ,e1)*s345**(-1)*s45**(-1) - 2.D0*cdot(p4,e5)*spstrng3(e2,e3,
     &    p145,e4,e1)*s145mmtsq**(-1)*s45**(-1) - 2.D0*cdot(p4,e5)*
     &    spstrng3(e2,p1345,e3,e4,e1)*s145mmtsq**(-1)*s45**(-1) + 4.D0*
     &    cdot(p5,e3)*cdot(p5,e4)*spstrng1(e2,e5,e1)*s345**(-1)*
     &    s45**(-1) - cdot(p5,e3)*cdot(e4,e5)*spstrng1(e2,p3,e1)*
     &    s345**(-1)*s45**(-1) + 3.D0*cdot(p5,e3)*cdot(e4,e5)*spstrng1(
     &    e2,p4,e1)*s345**(-1)*s45**(-1) - cdot(p5,e3)*cdot(e4,e5)*
     &    spstrng1(e2,p5,e1)*s345**(-1)*s45**(-1) + 2.D0*cdot(p5,e4)*
     &    cdot(e3,e5)*spstrng1(e2,p3,e1)*s345**(-1)*s45**(-1) - 2.D0*
     &    cdot(p5,e4)*cdot(e3,e5)*spstrng1(e2,p4,e1)*s345**(-1)*
     &    s45**(-1) - 2.D0*cdot(p5,e4)*cdot(e3,e5)*spstrng1(e2,p5,e1)*
     &    s345**(-1)*s45**(-1) + 2.D0*cdot(p5,e4)*spstrng3(e2,e3,p145,
     &    e5,e1)*s145mmtsq**(-1)*s45**(-1) + 2.D0*cdot(p5,e4)*spstrng3(
     &    e2,p1345,e3,e5,e1)*s145mmtsq**(-1)*s45**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt*
     & s1345mmtsq**(-1) * (  - cdot(e3,e4)*spstrng1(e2,e5,e1)*
     &    s345**(-1) + cdot(e3,e4)*spstrng3(e2,p3,p15,e5,e1)*s34**(-1)*
     &    s15mmtsq**(-1) - cdot(e3,e4)*spstrng3(e2,p4,p15,e5,e1)*
     &    s34**(-1)*s15mmtsq**(-1) + cdot(e3,e4)*spstrng3(e2,p1345,p3,
     &    e5,e1)*s34**(-1)*s15mmtsq**(-1) - cdot(e3,e4)*spstrng3(e2,
     &    p1345,p4,e5,e1)*s34**(-1)*s15mmtsq**(-1) + 2.D0*cdot(e3,e5)*
     &    spstrng1(e2,e4,e1)*s345**(-1) - cdot(e4,e5)*spstrng1(e2,e3,e1
     &    )*s345**(-1) + cdot(e4,e5)*spstrng3(e2,e3,p145,p4,e1)*
     &    s145mmtsq**(-1)*s45**(-1) - cdot(e4,e5)*spstrng3(e2,e3,p145,
     &    p5,e1)*s145mmtsq**(-1)*s45**(-1) + cdot(e4,e5)*spstrng3(e2,
     &    p1345,e3,p4,e1)*s145mmtsq**(-1)*s45**(-1) - cdot(e4,e5)*
     &    spstrng3(e2,p1345,e3,p5,e1)*s145mmtsq**(-1)*s45**(-1) + 
     &    spstrng5(e2,e3,p145,e4,p15,e5,e1)*s145mmtsq**(-1)*
     &    s15mmtsq**(-1) + spstrng5(e2,p1345,e3,e4,p15,e5,e1)*
     &    s145mmtsq**(-1)*s15mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt*
     & s1345mmtsq**(-1) * ( spstrng5(e2,p1345,e3,p145,e4,e5,e1)*
     &    s145mmtsq**(-1)*s15mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt*
     & s2345mmtsq**(-1) * (  - 2.D0*cdot(p3,p4)*cdot(e4,e5)*spstrng1(e2
     &    ,e3,e1)*s345**(-1)*s45**(-1) + 2.D0*cdot(p3,p5)*cdot(e3,e4)*
     &    spstrng1(e2,e5,e1)*s34**(-1)*s345**(-1) + 2.D0*cdot(p3,p5)*
     &    cdot(e4,e5)*spstrng1(e2,e3,e1)*s345**(-1)*s45**(-1) + 4.D0*
     &    cdot(p3,e4)*cdot(p3,e5)*spstrng1(e2,e3,e1)*s34**(-1)*
     &    s345**(-1) + 4.D0*cdot(p3,e4)*cdot(p4,e5)*spstrng1(e2,e3,e1)*
     &    s34**(-1)*s345**(-1) + 4.D0*cdot(p3,e4)*cdot(p4,e5)*spstrng1(
     &    e2,e3,e1)*s345**(-1)*s45**(-1) - 4.D0*cdot(p3,e4)*cdot(p5,e3)
     &    *spstrng1(e2,e5,e1)*s34**(-1)*s345**(-1) - 2.D0*cdot(p3,e4)*
     &    cdot(e3,e5)*spstrng1(e2,p3,e1)*s34**(-1)*s345**(-1) - 2.D0*
     &    cdot(p3,e4)*cdot(e3,e5)*spstrng1(e2,p4,e1)*s34**(-1)*
     &    s345**(-1) + 2.D0*cdot(p3,e4)*cdot(e3,e5)*spstrng1(e2,p5,e1)*
     &    s34**(-1)*s345**(-1) + 2.D0*cdot(p3,e4)*spstrng3(e2,e3,e5,
     &    p2345,e1)*s34**(-1)*s234mmtsq**(-1) + 2.D0*cdot(p3,e4)*
     &    spstrng3(e2,e3,p234,e5,e1)*s34**(-1)*s234mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt*
     & s2345mmtsq**(-1) * (  - 4.D0*cdot(p3,e5)*cdot(p4,e3)*spstrng1(e2
     &    ,e4,e1)*s34**(-1)*s345**(-1) - 4.D0*cdot(p3,e5)*cdot(p5,e4)*
     &    spstrng1(e2,e3,e1)*s345**(-1)*s45**(-1) - cdot(p3,e5)*cdot(e3
     &    ,e4)*spstrng1(e2,p3,e1)*s34**(-1)*s345**(-1) + 3.D0*cdot(p3,
     &    e5)*cdot(e3,e4)*spstrng1(e2,p4,e1)*s34**(-1)*s345**(-1) - 
     &    cdot(p3,e5)*cdot(e3,e4)*spstrng1(e2,p5,e1)*s34**(-1)*
     &    s345**(-1) - 2.D0*cdot(p4,p5)*cdot(e3,e4)*spstrng1(e2,e5,e1)*
     &    s34**(-1)*s345**(-1) - 4.D0*cdot(p4,e3)*cdot(p4,e5)*spstrng1(
     &    e2,e4,e1)*s34**(-1)*s345**(-1) - 4.D0*cdot(p4,e3)*cdot(p4,e5)
     &    *spstrng1(e2,e4,e1)*s345**(-1)*s45**(-1) + 4.D0*cdot(p4,e3)*
     &    cdot(p5,e4)*spstrng1(e2,e5,e1)*s34**(-1)*s345**(-1) + 4.D0*
     &    cdot(p4,e3)*cdot(p5,e4)*spstrng1(e2,e5,e1)*s345**(-1)*
     &    s45**(-1) + 2.D0*cdot(p4,e3)*cdot(e4,e5)*spstrng1(e2,p3,e1)*
     &    s34**(-1)*s345**(-1) + cdot(p4,e3)*cdot(e4,e5)*spstrng1(e2,p3
     &    ,e1)*s345**(-1)*s45**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt*
     & s2345mmtsq**(-1) * ( 2.D0*cdot(p4,e3)*cdot(e4,e5)*spstrng1(e2,p4
     &    ,e1)*s34**(-1)*s345**(-1) + cdot(p4,e3)*cdot(e4,e5)*spstrng1(
     &    e2,p4,e1)*s345**(-1)*s45**(-1) - 2.D0*cdot(p4,e3)*cdot(e4,e5)
     &    *spstrng1(e2,p5,e1)*s34**(-1)*s345**(-1) - 3.D0*cdot(p4,e3)*
     &    cdot(e4,e5)*spstrng1(e2,p5,e1)*s345**(-1)*s45**(-1) - 2.D0*
     &    cdot(p4,e3)*spstrng3(e2,e4,e5,p2345,e1)*s34**(-1)*
     &    s234mmtsq**(-1) - 2.D0*cdot(p4,e3)*spstrng3(e2,e4,p234,e5,e1)
     &    *s34**(-1)*s234mmtsq**(-1) - 4.D0*cdot(p4,e5)*cdot(p5,e3)*
     &    spstrng1(e2,e4,e1)*s345**(-1)*s45**(-1) - 3.D0*cdot(p4,e5)*
     &    cdot(e3,e4)*spstrng1(e2,p3,e1)*s34**(-1)*s345**(-1) - 2.D0*
     &    cdot(p4,e5)*cdot(e3,e4)*spstrng1(e2,p3,e1)*s345**(-1)*
     &    s45**(-1) + cdot(p4,e5)*cdot(e3,e4)*spstrng1(e2,p4,e1)*
     &    s34**(-1)*s345**(-1) + 2.D0*cdot(p4,e5)*cdot(e3,e4)*spstrng1(
     &    e2,p4,e1)*s345**(-1)*s45**(-1) + cdot(p4,e5)*cdot(e3,e4)*
     &    spstrng1(e2,p5,e1)*s34**(-1)*s345**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt*
     & s2345mmtsq**(-1) * ( 2.D0*cdot(p4,e5)*cdot(e3,e4)*spstrng1(e2,p5
     &    ,e1)*s345**(-1)*s45**(-1) + 2.D0*cdot(p4,e5)*spstrng3(e2,e3,
     &    e4,p2345,e1)*s45**(-1)*s23mmtsq**(-1) + 2.D0*cdot(p4,e5)*
     &    spstrng3(e2,e3,p23,e4,e1)*s45**(-1)*s23mmtsq**(-1) + 4.D0*
     &    cdot(p5,e3)*cdot(p5,e4)*spstrng1(e2,e5,e1)*s345**(-1)*
     &    s45**(-1) - cdot(p5,e3)*cdot(e4,e5)*spstrng1(e2,p3,e1)*
     &    s345**(-1)*s45**(-1) + 3.D0*cdot(p5,e3)*cdot(e4,e5)*spstrng1(
     &    e2,p4,e1)*s345**(-1)*s45**(-1) - cdot(p5,e3)*cdot(e4,e5)*
     &    spstrng1(e2,p5,e1)*s345**(-1)*s45**(-1) + 2.D0*cdot(p5,e4)*
     &    cdot(e3,e5)*spstrng1(e2,p3,e1)*s345**(-1)*s45**(-1) - 2.D0*
     &    cdot(p5,e4)*cdot(e3,e5)*spstrng1(e2,p4,e1)*s345**(-1)*
     &    s45**(-1) - 2.D0*cdot(p5,e4)*cdot(e3,e5)*spstrng1(e2,p5,e1)*
     &    s345**(-1)*s45**(-1) - 2.D0*cdot(p5,e4)*spstrng3(e2,e3,e5,
     &    p2345,e1)*s45**(-1)*s23mmtsq**(-1) - 2.D0*cdot(p5,e4)*
     &    spstrng3(e2,e3,p23,e5,e1)*s45**(-1)*s23mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt*
     & s2345mmtsq**(-1) * (  - cdot(e3,e4)*spstrng1(e2,e5,e1)*
     &    s345**(-1) - cdot(e3,e4)*spstrng3(e2,p3,e5,p2345,e1)*
     &    s34**(-1)*s234mmtsq**(-1) - cdot(e3,e4)*spstrng3(e2,p3,p234,
     &    e5,e1)*s34**(-1)*s234mmtsq**(-1) + cdot(e3,e4)*spstrng3(e2,p4
     &    ,e5,p2345,e1)*s34**(-1)*s234mmtsq**(-1) + cdot(e3,e4)*
     &    spstrng3(e2,p4,p234,e5,e1)*s34**(-1)*s234mmtsq**(-1) + 2.D0*
     &    cdot(e3,e5)*spstrng1(e2,e4,e1)*s345**(-1) - cdot(e4,e5)*
     &    spstrng1(e2,e3,e1)*s345**(-1) - cdot(e4,e5)*spstrng3(e2,e3,p4
     &    ,p2345,e1)*s45**(-1)*s23mmtsq**(-1) + cdot(e4,e5)*spstrng3(e2
     &    ,e3,p5,p2345,e1)*s45**(-1)*s23mmtsq**(-1) - cdot(e4,e5)*
     &    spstrng3(e2,e3,p23,p4,e1)*s45**(-1)*s23mmtsq**(-1) + cdot(e4,
     &    e5)*spstrng3(e2,e3,p23,p5,e1)*s45**(-1)*s23mmtsq**(-1) + 
     &    spstrng5(e2,e3,e4,p234,e5,p2345,e1)*s234mmtsq**(-1)*
     &    s23mmtsq**(-1) + spstrng5(e2,e3,p23,e4,e5,p2345,e1)*
     &    s234mmtsq**(-1)*s23mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt*
     & s2345mmtsq**(-1) * ( spstrng5(e2,e3,p23,e4,p234,e5,e1)*
     &    s234mmtsq**(-1)*s23mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt * ( 2.D0*
     &    cdot(p3,e4)*spstrng3(e2,e3,p234,e5,e1)*s34**(-1)*
     &    s234mmtsq**(-1)*s15mmtsq**(-1) - 2.D0*cdot(p3,e4)*spstrng3(e2
     &    ,e3,p15,e5,e1)*s34**(-1)*s234mmtsq**(-1)*s15mmtsq**(-1) - 2.D0
     &    *cdot(p4,e3)*spstrng3(e2,e4,p234,e5,e1)*s34**(-1)*
     &    s234mmtsq**(-1)*s15mmtsq**(-1) + 2.D0*cdot(p4,e3)*spstrng3(e2
     &    ,e4,p15,e5,e1)*s34**(-1)*s234mmtsq**(-1)*s15mmtsq**(-1) + 2.D0
     &    *cdot(p4,e5)*spstrng3(e2,e3,p23,e4,e1)*s145mmtsq**(-1)*
     &    s45**(-1)*s23mmtsq**(-1) - 2.D0*cdot(p4,e5)*spstrng3(e2,e3,
     &    p145,e4,e1)*s145mmtsq**(-1)*s45**(-1)*s23mmtsq**(-1) - 2.D0*
     &    cdot(p5,e4)*spstrng3(e2,e3,p23,e5,e1)*s145mmtsq**(-1)*
     &    s45**(-1)*s23mmtsq**(-1) + 2.D0*cdot(p5,e4)*spstrng3(e2,e3,
     &    p145,e5,e1)*s145mmtsq**(-1)*s45**(-1)*s23mmtsq**(-1) - cdot(
     &    e3,e4)*spstrng3(e2,p3,p234,e5,e1)*s34**(-1)*s234mmtsq**(-1)*
     &    s15mmtsq**(-1) + cdot(e3,e4)*spstrng3(e2,p3,p15,e5,e1)*
     &    s34**(-1)*s234mmtsq**(-1)*s15mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt * ( cdot(e3,
     &    e4)*spstrng3(e2,p4,p234,e5,e1)*s34**(-1)*s234mmtsq**(-1)*
     &    s15mmtsq**(-1) - cdot(e3,e4)*spstrng3(e2,p4,p15,e5,e1)*
     &    s34**(-1)*s234mmtsq**(-1)*s15mmtsq**(-1) - cdot(e4,e5)*
     &    spstrng3(e2,e3,p23,p4,e1)*s145mmtsq**(-1)*s45**(-1)*
     &    s23mmtsq**(-1) + cdot(e4,e5)*spstrng3(e2,e3,p23,p5,e1)*
     &    s145mmtsq**(-1)*s45**(-1)*s23mmtsq**(-1) + cdot(e4,e5)*
     &    spstrng3(e2,e3,p145,p4,e1)*s145mmtsq**(-1)*s45**(-1)*
     &    s23mmtsq**(-1) - cdot(e4,e5)*spstrng3(e2,e3,p145,p5,e1)*
     &    s145mmtsq**(-1)*s45**(-1)*s23mmtsq**(-1) - spstrng5(e2,e3,e4,
     &    p234,p15,e5,e1)*s234mmtsq**(-1)*s15mmtsq**(-1)*s23mmtsq**(-1)
     &     + spstrng5(e2,e3,p23,e4,p234,e5,e1)*s234mmtsq**(-1)*
     &    s15mmtsq**(-1)*s23mmtsq**(-1) - spstrng5(e2,e3,p23,e4,p15,e5,
     &    e1)*s234mmtsq**(-1)*s15mmtsq**(-1)*s23mmtsq**(-1) - spstrng5(
     &    e2,e3,p23,e4,p15,e5,e1)*s145mmtsq**(-1)*s15mmtsq**(-1)*
     &    s23mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt * (  - 
     &    spstrng5(e2,e3,p23,p145,e4,e5,e1)*s145mmtsq**(-1)*
     &    s15mmtsq**(-1)*s23mmtsq**(-1) + spstrng5(e2,e3,p145,e4,p15,e5
     &    ,e1)*s145mmtsq**(-1)*s15mmtsq**(-1)*s23mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt**2*
     & s1345mmtsq**(-1) * ( 2.D0*cdot(p3,e4)*spstrng2(e2,e3,e5,e1)*
     &    s34**(-1)*s15mmtsq**(-1) - 2.D0*cdot(p4,e3)*spstrng2(e2,e4,e5
     &    ,e1)*s34**(-1)*s15mmtsq**(-1) + 2.D0*cdot(p4,e5)*spstrng2(e2,
     &    e3,e4,e1)*s145mmtsq**(-1)*s45**(-1) - 2.D0*cdot(p5,e4)*
     &    spstrng2(e2,e3,e5,e1)*s145mmtsq**(-1)*s45**(-1) - cdot(e3,e4)
     &    *spstrng2(e2,p3,e5,e1)*s34**(-1)*s15mmtsq**(-1) + cdot(e3,e4)
     &    *spstrng2(e2,p4,e5,e1)*s34**(-1)*s15mmtsq**(-1) - cdot(e4,e5)
     &    *spstrng2(e2,e3,p4,e1)*s145mmtsq**(-1)*s45**(-1) + cdot(e4,e5
     &    )*spstrng2(e2,e3,p5,e1)*s145mmtsq**(-1)*s45**(-1) - spstrng4(
     &    e2,e3,e4,p15,e5,e1)*s145mmtsq**(-1)*s15mmtsq**(-1) - 
     &    spstrng4(e2,e3,p145,e4,e5,e1)*s145mmtsq**(-1)*s15mmtsq**(-1)
     &     - spstrng4(e2,p1345,e3,e4,e5,e1)*s145mmtsq**(-1)*
     &    s15mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt**2*
     & s2345mmtsq**(-1) * ( 2.D0*cdot(p3,e4)*spstrng2(e2,e3,e5,e1)*
     &    s34**(-1)*s234mmtsq**(-1) - 2.D0*cdot(p4,e3)*spstrng2(e2,e4,
     &    e5,e1)*s34**(-1)*s234mmtsq**(-1) + 2.D0*cdot(p4,e5)*spstrng2(
     &    e2,e3,e4,e1)*s45**(-1)*s23mmtsq**(-1) - 2.D0*cdot(p5,e4)*
     &    spstrng2(e2,e3,e5,e1)*s45**(-1)*s23mmtsq**(-1) - cdot(e3,e4)*
     &    spstrng2(e2,p3,e5,e1)*s34**(-1)*s234mmtsq**(-1) + cdot(e3,e4)
     &    *spstrng2(e2,p4,e5,e1)*s34**(-1)*s234mmtsq**(-1) - cdot(e4,e5
     &    )*spstrng2(e2,e3,p4,e1)*s45**(-1)*s23mmtsq**(-1) + cdot(e4,e5
     &    )*spstrng2(e2,e3,p5,e1)*s45**(-1)*s23mmtsq**(-1) + spstrng4(
     &    e2,e3,e4,e5,p2345,e1)*s234mmtsq**(-1)*s23mmtsq**(-1) + 
     &    spstrng4(e2,e3,e4,p234,e5,e1)*s234mmtsq**(-1)*s23mmtsq**(-1)
     &     + spstrng4(e2,e3,p23,e4,e5,e1)*s234mmtsq**(-1)*
     &    s23mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt**2 * ( 2.D0*
     &    cdot(p3,e4)*spstrng2(e2,e3,e5,e1)*s34**(-1)*s234mmtsq**(-1)*
     &    s15mmtsq**(-1) - 2.D0*cdot(p4,e3)*spstrng2(e2,e4,e5,e1)*
     &    s34**(-1)*s234mmtsq**(-1)*s15mmtsq**(-1) + 2.D0*cdot(p4,e5)*
     &    spstrng2(e2,e3,e4,e1)*s145mmtsq**(-1)*s45**(-1)*
     &    s23mmtsq**(-1) - 2.D0*cdot(p5,e4)*spstrng2(e2,e3,e5,e1)*
     &    s145mmtsq**(-1)*s45**(-1)*s23mmtsq**(-1) - cdot(e3,e4)*
     &    spstrng2(e2,p3,e5,e1)*s34**(-1)*s234mmtsq**(-1)*
     &    s15mmtsq**(-1) + cdot(e3,e4)*spstrng2(e2,p4,e5,e1)*s34**(-1)*
     &    s234mmtsq**(-1)*s15mmtsq**(-1) - cdot(e4,e5)*spstrng2(e2,e3,
     &    p4,e1)*s145mmtsq**(-1)*s45**(-1)*s23mmtsq**(-1) + cdot(e4,e5)
     &    *spstrng2(e2,e3,p5,e1)*s145mmtsq**(-1)*s45**(-1)*
     &    s23mmtsq**(-1) + spstrng4(e2,e3,e4,p234,e5,e1)*
     &    s234mmtsq**(-1)*s15mmtsq**(-1)*s23mmtsq**(-1) - spstrng4(e2,
     &    e3,e4,p15,e5,e1)*s234mmtsq**(-1)*s15mmtsq**(-1)*
     &    s23mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt**2 * (  - 
     &    spstrng4(e2,e3,e4,p15,e5,e1)*s145mmtsq**(-1)*s15mmtsq**(-1)*
     &    s23mmtsq**(-1) + spstrng4(e2,e3,p23,e4,e5,e1)*s234mmtsq**(-1)
     &    *s15mmtsq**(-1)*s23mmtsq**(-1) + spstrng4(e2,e3,p23,e4,e5,e1)
     &    *s145mmtsq**(-1)*s15mmtsq**(-1)*s23mmtsq**(-1) - spstrng4(e2,
     &    e3,p145,e4,e5,e1)*s145mmtsq**(-1)*s15mmtsq**(-1)*
     &    s23mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt**3*
     & s1345mmtsq**(-1) * ( spstrng3(e2,e3,e4,e5,e1)*s145mmtsq**(-1)*
     &    s15mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt**3*
     & s2345mmtsq**(-1) * ( spstrng3(e2,e3,e4,e5,e1)*s234mmtsq**(-1)*
     &    s23mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) + mt**3 * ( 
     &    spstrng3(e2,e3,e4,e5,e1)*s234mmtsq**(-1)*s15mmtsq**(-1)*
     &    s23mmtsq**(-1) + spstrng3(e2,e3,e4,e5,e1)*s145mmtsq**(-1)*
     &    s15mmtsq**(-1)*s23mmtsq**(-1) )
      ampABC(h1,h2,h3,h4,h5) = ampABC(h1,h2,h3,h4,h5) - 2.D0*cdot(p3,e4
     & )*spstrng4(e2,e3,p234,p15,e5,e1)*s34**(-1)*s234mmtsq**(-1)*
     & s15mmtsq**(-1) + 2.D0*cdot(p4,e3)*spstrng4(e2,e4,p234,p15,e5,e1)
     &    *s34**(-1)*s234mmtsq**(-1)*s15mmtsq**(-1) - 2.D0*cdot(p4,e5)*
     &    spstrng4(e2,e3,p23,p145,e4,e1)*s145mmtsq**(-1)*s45**(-1)*
     &    s23mmtsq**(-1) + 2.D0*cdot(p5,e4)*spstrng4(e2,e3,p23,p145,e5,
     &    e1)*s145mmtsq**(-1)*s45**(-1)*s23mmtsq**(-1) + cdot(e3,e4)*
     &    spstrng4(e2,p3,p234,p15,e5,e1)*s34**(-1)*s234mmtsq**(-1)*
     &    s15mmtsq**(-1) - cdot(e3,e4)*spstrng4(e2,p4,p234,p15,e5,e1)*
     &    s34**(-1)*s234mmtsq**(-1)*s15mmtsq**(-1) + cdot(e4,e5)*
     &    spstrng4(e2,e3,p23,p145,p4,e1)*s145mmtsq**(-1)*s45**(-1)*
     &    s23mmtsq**(-1) - cdot(e4,e5)*spstrng4(e2,e3,p23,p145,p5,e1)*
     &    s145mmtsq**(-1)*s45**(-1)*s23mmtsq**(-1) - spstrng6(e2,e3,p23
     &    ,e4,p234,p15,e5,e1)*s234mmtsq**(-1)*s15mmtsq**(-1)*
     &    s23mmtsq**(-1) + spstrng6(e2,e3,p23,p145,e4,p15,e5,e1)*
     &    s145mmtsq**(-1)*s15mmtsq**(-1)*s23mmtsq**(-1)

      enddo
      enddo
      enddo
      enddo
      enddo

      do h1=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
!      h2=2
      ampABC(h1,2,h3,h4,h5)=-(-1)**h1*conjg(ampABC(3-h1,1,h3,h4,h5))
      enddo
      enddo
      enddo
      enddo
      return
      end
