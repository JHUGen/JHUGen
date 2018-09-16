      subroutine ttggZamp(q,q1,q2,q3,q4,q5,q6,eta1,eta2,ampABL,ampABR)
      implicit none
      include 'types.f'
      
C---  routine written by ttZgg/Form/losq.frm
C---  Author R.K. Ellis, May 2013
C---  Amplitude for one color ordering of
C---  t~(p1)+t(p2)+g(p3)+g(p4)+l(p5)+z(p6)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      real(dp):: q(mxpart,4),s134mmtsq,s234mmtsq,s134,s234,s34
      real(dp):: s14,s23,s14mmtsq,s23mmtsq,s56
      complex(dp):: spstrng0,spstrng1,spstrng2,
     & spstrng3,spstrng4,cdot,
     & e1(4),e2(4),e3(4),e4(4),p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),
     & et1(4),et2(4),
     & e5m(4),e5p(4),e6m(4),e6p(4),e5Cm(4),e5Cp(4),e6Cm(4),e6Cp(4),
     & e1m(4),e2m(4),e3m(4),e4m(4),e1p(4),e2p(4),e3p(4),e4p(4),
     & p134(4),p234(4),p14(4),p23(4),p34(4),p56(4),
     & ampABL(2,2,2,2,2),ampABR(2,2,2,2,2)
      integer:: q1,q2,q3,q4,q5,q6,eta1,eta2,h1,h2,h3,h4
      p1(:)=cplx2(q(q1,:))
      p2(:)=cplx2(q(q2,:))
      p3(:)=cplx2(q(q3,:))
      p4(:)=cplx2(q(q4,:))
      p5(:)=cplx2(q(q5,:))
      p6(:)=cplx2(q(q6,:))
      et1(:)=cplx2(q(eta1,:))
      et2(:)=cplx2(q(eta2,:))
      p14(:)=cplx2(q(q1,:)+q(q4,:))
      p23(:)=cplx2(q(q2,:)+q(q3,:))
      p34(:)=cplx2(q(q3,:)+q(q4,:))
      p56(:)=cplx2(q(q5,:)+q(q6,:))
      p134(:)=cplx2(q(q1,:)+q(q3,:)+q(q4,:))
      p234(:)=cplx2(q(q2,:)+q(q3,:)+q(q4,:))
      s14=real(p14(4)**2-p14(1)**2-p14(2)**2-p14(3)**2)
      s23=real(p23(4)**2-p23(1)**2-p23(2)**2-p23(3)**2)
      s34=real(p34(4)**2-p34(1)**2-p34(2)**2-p34(3)**2)
      s56=real(p56(4)**2-p56(1)**2-p56(2)**2-p56(3)**2)
      s134=real(p134(4)**2-p134(1)**2-p134(2)**2-p134(3)**2)
      s234=real(p234(4)**2-p234(1)**2-p234(2)**2-p234(3)**2)
      s134mmtsq=s134-mt**2
      s234mmtsq=s234-mt**2
      s23mmtsq=s23-mt**2
      s14mmtsq=s14-mt**2
      h1=1
      call VKlSt(p1,mt,et1,h1,e1p)
      call UbKlSt(p2,mt,et2,h1,e2p)
      call pol_real(p3,h1,e3p)
      call pol_real(p4,h1,e4p)
      call ubarspinor0(p5,h1,e5p)
      call ubarspinor0(p6,h1,e6p)
      call uspinor0(p6,h1,e6Cp)
      call uspinor0(p5,h1,e5Cp)
      h1=-1
      call VklSt(p1,mt,et1,h1,e1m)
      call UbKlst(p2,mt,et2,h1,e2m)
      call pol_real(p3,h1,e3m)
      call pol_real(p4,h1,e4m)
      call ubarspinor0(p5,h1,e5m)
      call ubarspinor0(p6,h1,e6m)
      call uspinor0(p6,h1,e6Cm)
      call uspinor0(p5,h1,e5Cm)
 
      do h1=1,2
      if (h1 == 1) e1(:)=e1m(:)
      if (h1 == 2) e1(:)=e1p(:)
      do h2=1,2
      if (h2 == 1) e2(:)=e2m(:)
      if (h2 == 2) e2(:)=e2p(:)
      do h3=1,2
      if (h3 == 1) e3(:)=e3m(:)
      if (h3 == 2) e3(:)=e3p(:)
      do h4=1,2
      if (h4 == 1) e4(:)=e4m(:)
      if (h4 == 2) e4(:)=e4p(:)
      ampABL(h1,h2,h3,h4,1)= + s34**(-1)*s234mmtsq**(-1) * ( 4.D0*cdot(
     &    p3,e4)*spstrng0(e6p,e1)*spstrng2(e2,e3,p234,e5Cp)*s56**(-1)
     &     - 4.D0*cdot(p4,e3)*spstrng0(e6p,e1)*spstrng2(e2,e4,p234,e5Cp
     &    )*s56**(-1) + 4.D0*cdot(e3,e4)*spstrng0(e6p,e1)*spstrng2(e2,
     &    p4,p234,e5Cp)*s56**(-1) - 2.D0*cdot(e3,e4)*spstrng0(e6p,e1)*
     &    spstrng2(e2,p34,p234,e5Cp)*s56**(-1) )
      ampABL(h1,h2,h3,h4,1) = ampABL(h1,h2,h3,h4,1) + s34**(-1)*
     & s234mmtsq**(-1)*mt * ( 4.D0*cdot(p3,e4)*spstrng0(e6p,e1)*
     &    spstrng1(e2,e3,e5Cp)*s56**(-1) - 4.D0*cdot(p4,e3)*spstrng0(
     &    e6p,e1)*spstrng1(e2,e4,e5Cp)*s56**(-1) + 4.D0*cdot(e3,e4)*
     &    spstrng0(e6p,e1)*spstrng1(e2,p4,e5Cp)*s56**(-1) - 2.D0*cdot(
     &    e3,e4)*spstrng0(e6p,e1)*spstrng1(e2,p34,e5Cp)*s56**(-1) )
      ampABL(h1,h2,h3,h4,1) = ampABL(h1,h2,h3,h4,1) + s34**(-1)*
     & s134mmtsq**(-1) * (  - 4.D0*cdot(p3,e4)*spstrng0(e2,e5Cp)*
     &    spstrng2(e6p,p134,e3,e1)*s56**(-1) + 4.D0*cdot(p4,e3)*
     &    spstrng0(e2,e5Cp)*spstrng2(e6p,p134,e4,e1)*s56**(-1) - 4.D0*
     &    cdot(e3,e4)*spstrng0(e2,e5Cp)*spstrng2(e6p,p134,p4,e1)*
     &    s56**(-1) + 2.D0*cdot(e3,e4)*spstrng0(e2,e5Cp)*spstrng2(e6p,
     &    p134,p34,e1)*s56**(-1) )
      ampABL(h1,h2,h3,h4,1) = ampABL(h1,h2,h3,h4,1) + s34**(-1)*
     & s134mmtsq**(-1)*mt * ( 4.D0*cdot(p3,e4)*spstrng0(e2,e5Cp)*
     &    spstrng1(e6p,e3,e1)*s56**(-1) - 4.D0*cdot(p4,e3)*spstrng0(e2,
     &    e5Cp)*spstrng1(e6p,e4,e1)*s56**(-1) + 4.D0*cdot(e3,e4)*
     &    spstrng0(e2,e5Cp)*spstrng1(e6p,p4,e1)*s56**(-1) - 2.D0*cdot(
     &    e3,e4)*spstrng0(e2,e5Cp)*spstrng1(e6p,p34,e1)*s56**(-1) )
      ampABL(h1,h2,h3,h4,1) = ampABL(h1,h2,h3,h4,1) + s234mmtsq**(-1)
     &  * ( 2.D0*spstrng0(e6p,e1)*spstrng4(e2,e3,p23,e4,p234,e5Cp)*
     &    s23mmtsq**(-1)*s56**(-1) )
      ampABL(h1,h2,h3,h4,1) = ampABL(h1,h2,h3,h4,1) + s234mmtsq**(-1)*
     & mt * ( 2.D0*spstrng0(e6p,e1)*spstrng3(e2,e3,e4,p234,e5Cp)*
     &    s23mmtsq**(-1)*s56**(-1) + 2.D0*spstrng0(e6p,e1)*spstrng3(e2,
     &    e3,p23,e4,e5Cp)*s23mmtsq**(-1)*s56**(-1) )
      ampABL(h1,h2,h3,h4,1) = ampABL(h1,h2,h3,h4,1) + s234mmtsq**(-1)*
     & mt**2 * ( 2.D0*spstrng0(e6p,e1)*spstrng2(e2,e3,e4,e5Cp)*
     &    s23mmtsq**(-1)*s56**(-1) )
      ampABL(h1,h2,h3,h4,1) = ampABL(h1,h2,h3,h4,1) + s134mmtsq**(-1)
     &  * ( 2.D0*spstrng0(e2,e5Cp)*spstrng4(e6p,p134,e3,p14,e4,e1)*
     &    s14mmtsq**(-1)*s56**(-1) )
      ampABL(h1,h2,h3,h4,1) = ampABL(h1,h2,h3,h4,1) + s134mmtsq**(-1)*
     & mt * (  - 2.D0*spstrng0(e2,e5Cp)*spstrng3(e6p,e3,p14,e4,e1)*
     &    s14mmtsq**(-1)*s56**(-1) - 2.D0*spstrng0(e2,e5Cp)*spstrng3(
     &    e6p,p134,e3,e4,e1)*s14mmtsq**(-1)*s56**(-1) )
      ampABL(h1,h2,h3,h4,1) = ampABL(h1,h2,h3,h4,1) + s134mmtsq**(-1)*
     & mt**2 * ( 2.D0*spstrng0(e2,e5Cp)*spstrng2(e6p,e3,e4,e1)*
     &    s14mmtsq**(-1)*s56**(-1) )
      ampABL(h1,h2,h3,h4,1) = ampABL(h1,h2,h3,h4,1) + mt * (  - 2.D0*
     &    spstrng1(e2,e3,e5Cp)*spstrng2(e6p,p14,e4,e1)*s14mmtsq**(-1)*
     &    s23mmtsq**(-1)*s56**(-1) + 2.D0*spstrng1(e6p,e4,e1)*spstrng2(
     &    e2,e3,p23,e5Cp)*s14mmtsq**(-1)*s23mmtsq**(-1)*s56**(-1) )
      ampABL(h1,h2,h3,h4,1) = ampABL(h1,h2,h3,h4,1) + mt**2 * ( 2.D0*
     &    spstrng1(e2,e3,e5Cp)*spstrng1(e6p,e4,e1)*s14mmtsq**(-1)*
     &    s23mmtsq**(-1)*s56**(-1) )
      ampABL(h1,h2,h3,h4,1) = ampABL(h1,h2,h3,h4,1) - 2.D0*spstrng2(e2,
     & e3,p23,e5Cp)*spstrng2(e6p,p14,e4,e1)*s14mmtsq**(-1)*
     & s23mmtsq**(-1)*s56**(-1)

      ampABR(h1,h2,h3,h4,1)= + s34**(-1)*s234mmtsq**(-1) * ( 4.D0*cdot(
     &    p3,e4)*spstrng0(e5m,e1)*spstrng2(e2,e3,p234,e6Cm)*s56**(-1)
     &     - 4.D0*cdot(p4,e3)*spstrng0(e5m,e1)*spstrng2(e2,e4,p234,e6Cm
     &    )*s56**(-1) + 4.D0*cdot(e3,e4)*spstrng0(e5m,e1)*spstrng2(e2,
     &    p4,p234,e6Cm)*s56**(-1) - 2.D0*cdot(e3,e4)*spstrng0(e5m,e1)*
     &    spstrng2(e2,p34,p234,e6Cm)*s56**(-1) )
      ampABR(h1,h2,h3,h4,1) = ampABR(h1,h2,h3,h4,1) + s34**(-1)*
     & s234mmtsq**(-1)*mt * ( 4.D0*cdot(p3,e4)*spstrng0(e5m,e1)*
     &    spstrng1(e2,e3,e6Cm)*s56**(-1) - 4.D0*cdot(p4,e3)*spstrng0(
     &    e5m,e1)*spstrng1(e2,e4,e6Cm)*s56**(-1) + 4.D0*cdot(e3,e4)*
     &    spstrng0(e5m,e1)*spstrng1(e2,p4,e6Cm)*s56**(-1) - 2.D0*cdot(
     &    e3,e4)*spstrng0(e5m,e1)*spstrng1(e2,p34,e6Cm)*s56**(-1) )
      ampABR(h1,h2,h3,h4,1) = ampABR(h1,h2,h3,h4,1) + s34**(-1)*
     & s134mmtsq**(-1) * (  - 4.D0*cdot(p3,e4)*spstrng0(e2,e6Cm)*
     &    spstrng2(e5m,p134,e3,e1)*s56**(-1) + 4.D0*cdot(p4,e3)*
     &    spstrng0(e2,e6Cm)*spstrng2(e5m,p134,e4,e1)*s56**(-1) - 4.D0*
     &    cdot(e3,e4)*spstrng0(e2,e6Cm)*spstrng2(e5m,p134,p4,e1)*
     &    s56**(-1) + 2.D0*cdot(e3,e4)*spstrng0(e2,e6Cm)*spstrng2(e5m,
     &    p134,p34,e1)*s56**(-1) )
      ampABR(h1,h2,h3,h4,1) = ampABR(h1,h2,h3,h4,1) + s34**(-1)*
     & s134mmtsq**(-1)*mt * ( 4.D0*cdot(p3,e4)*spstrng0(e2,e6Cm)*
     &    spstrng1(e5m,e3,e1)*s56**(-1) - 4.D0*cdot(p4,e3)*spstrng0(e2,
     &    e6Cm)*spstrng1(e5m,e4,e1)*s56**(-1) + 4.D0*cdot(e3,e4)*
     &    spstrng0(e2,e6Cm)*spstrng1(e5m,p4,e1)*s56**(-1) - 2.D0*cdot(
     &    e3,e4)*spstrng0(e2,e6Cm)*spstrng1(e5m,p34,e1)*s56**(-1) )
      ampABR(h1,h2,h3,h4,1) = ampABR(h1,h2,h3,h4,1) + s234mmtsq**(-1)
     &  * ( 2.D0*spstrng0(e5m,e1)*spstrng4(e2,e3,p23,e4,p234,e6Cm)*
     &    s23mmtsq**(-1)*s56**(-1) )
      ampABR(h1,h2,h3,h4,1) = ampABR(h1,h2,h3,h4,1) + s234mmtsq**(-1)*
     & mt * ( 2.D0*spstrng0(e5m,e1)*spstrng3(e2,e3,e4,p234,e6Cm)*
     &    s23mmtsq**(-1)*s56**(-1) + 2.D0*spstrng0(e5m,e1)*spstrng3(e2,
     &    e3,p23,e4,e6Cm)*s23mmtsq**(-1)*s56**(-1) )
      ampABR(h1,h2,h3,h4,1) = ampABR(h1,h2,h3,h4,1) + s234mmtsq**(-1)*
     & mt**2 * ( 2.D0*spstrng0(e5m,e1)*spstrng2(e2,e3,e4,e6Cm)*
     &    s23mmtsq**(-1)*s56**(-1) )
      ampABR(h1,h2,h3,h4,1) = ampABR(h1,h2,h3,h4,1) + s134mmtsq**(-1)
     &  * ( 2.D0*spstrng0(e2,e6Cm)*spstrng4(e5m,p134,e3,p14,e4,e1)*
     &    s14mmtsq**(-1)*s56**(-1) )
      ampABR(h1,h2,h3,h4,1) = ampABR(h1,h2,h3,h4,1) + s134mmtsq**(-1)*
     & mt * (  - 2.D0*spstrng0(e2,e6Cm)*spstrng3(e5m,e3,p14,e4,e1)*
     &    s14mmtsq**(-1)*s56**(-1) - 2.D0*spstrng0(e2,e6Cm)*spstrng3(
     &    e5m,p134,e3,e4,e1)*s14mmtsq**(-1)*s56**(-1) )
      ampABR(h1,h2,h3,h4,1) = ampABR(h1,h2,h3,h4,1) + s134mmtsq**(-1)*
     & mt**2 * ( 2.D0*spstrng0(e2,e6Cm)*spstrng2(e5m,e3,e4,e1)*
     &    s14mmtsq**(-1)*s56**(-1) )
      ampABR(h1,h2,h3,h4,1) = ampABR(h1,h2,h3,h4,1) + mt * ( 2.D0*
     &    spstrng1(e5m,e4,e1)*spstrng2(e2,e3,p23,e6Cm)*s14mmtsq**(-1)*
     &    s23mmtsq**(-1)*s56**(-1) - 2.D0*spstrng1(e2,e3,e6Cm)*
     &    spstrng2(e5m,p14,e4,e1)*s14mmtsq**(-1)*s23mmtsq**(-1)*
     &    s56**(-1) )
      ampABR(h1,h2,h3,h4,1) = ampABR(h1,h2,h3,h4,1) + mt**2 * ( 2.D0*
     &    spstrng1(e5m,e4,e1)*spstrng1(e2,e3,e6Cm)*s14mmtsq**(-1)*
     &    s23mmtsq**(-1)*s56**(-1) )
      ampABR(h1,h2,h3,h4,1) = ampABR(h1,h2,h3,h4,1) - 2.D0*spstrng2(e5m
     & ,p14,e4,e1)*spstrng2(e2,e3,p23,e6Cm)*s14mmtsq**(-1)*
     & s23mmtsq**(-1)*s56**(-1)

      enddo
      enddo
      enddo
      enddo

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      ampABR(h1,h2,h3,h4,2)=
     & -(-1)**(h1+h2)*conjg(ampABL(3-h1,3-h2,h3,h4,1))
      ampABL(h1,h2,h3,h4,2)=
     & -(-1)**(h1+h2)*conjg(ampABR(3-h1,3-h2,h3,h4,1))
      enddo
      enddo
      enddo
      enddo
      return
      end
