      subroutine ttqqgHampsq(q,q1,q2,q3,q4,q5,eta1,eta2,ampsq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      real(dp):: q(mxpart,4),s134mmtsq,s234mmtsq,
     & s134,s234,s34,s35,s45,s345,s15,s25,s15mmtsq,
     & s25mmtsq,s1345mmtsq,
     & s2345mmtsq,s1345,s2345,ampsq
      complex(dp):: spstrng0,spstrng1,spstrng2,
     & spstrng3,cdot,m(4),
     & p1(4),p2(4),p3(4),p4(4),p5(4),et1(4),et2(4),
     & e1C(4),e1Cp(4),e1Cm(4),e2(4),e2p(4),e2m(4),e5m(4),e5p(4),e5(4),
     & e3(4),e3C(4),e4(4),e4C(4),
     & e3m(4),e4m(4),e3p(4),e4p(4),e3Cm(4),e4Cm(4),e3Cp(4),e4Cp(4),
     & p234(4),p34(4),p25(4),p35(4),p45(4),p345(4),p1345(4),
     & p2345(4),p15(4),p134(4),T34x21(2,2,2,2,2),
     & T21x34(2,2,2,2,2),T24x31(2,2,2,2,2),T31x24(2,2,2,2,2)
      integer:: q1,q2,q3,q4,q5,eta1,eta2,h1,h2,h3,h4,h5
      p1(:)=cplx2(q(q1,:))
      p2(:)=cplx2(q(q2,:))
      p3(:)=cplx2(q(q3,:))
      p4(:)=cplx2(q(q4,:))
      p5(:)=cplx2(q(q5,:))
      et1(:)=cplx2(q(eta1,:))
      et2(:)=cplx2(q(eta2,:))
      p15(:)=cplx2(q(q1,:)+q(q5,:))
      p25(:)=cplx2(q(q2,:)+q(q5,:))
      p34(:)=cplx2(q(q3,:)+q(q4,:))
      p45(:)=cplx2(q(q4,:)+q(q5,:))
      p35(:)=cplx2(q(q3,:)+q(q5,:))
      p345(:)=cplx2(q(q3,:)+q(q4,:)+q(q5,:))
      p1345(:)=cplx2(q(q1,:)+q(q3,:)+q(q4,:)+q(q5,:))
      p2345(:)=cplx2(q(q2,:)+q(q3,:)+q(q4,:)+q(q5,:))
      p234(:)=cplx2(q(q2,:)+q(q3,:)+q(q4,:))
      p134(:)=cplx2(q(q1,:)+q(q3,:)+q(q4,:))
      s45=real(p45(4)**2-p45(1)**2-p45(2)**2-p45(3)**2)
      s15=real(p15(4)**2-p15(1)**2-p15(2)**2-p15(3)**2)
      s25=real(p25(4)**2-p25(1)**2-p25(2)**2-p25(3)**2)
      s34=real(p34(4)**2-p34(1)**2-p34(2)**2-p34(3)**2)
      s35=real(p35(4)**2-p35(1)**2-p35(2)**2-p35(3)**2)
      s134=real(p134(4)**2-p134(1)**2-p134(2)**2-p134(3)**2)
      s234=real(p234(4)**2-p234(1)**2-p234(2)**2-p234(3)**2)
      s345=real(p345(4)**2-p345(1)**2-p345(2)**2-p345(3)**2)
      s1345=real(p1345(4)**2-p1345(1)**2-p1345(2)**2-p1345(3)**2)
      s2345=real(p2345(4)**2-p2345(1)**2-p2345(2)**2-p2345(3)**2)
      s15mmtsq=s15-mt**2
      s25mmtsq=s25-mt**2
      s134mmtsq=s134-mt**2
      s234mmtsq=s234-mt**2
      s1345mmtsq=s1345-mt**2
      s2345mmtsq=s2345-mt**2
      h1=-1
      call VklSt(p1,mt,p3,h1,e1Cm)
      call UbKlst(p2,mt,p3,h1,e2m)
      call ubarspinor0(p3,h1,e3m)
      call uspinor0(p3,-h1,e3Cm)
      call ubarspinor0(p4,h1,e4m)
      call uspinor0(p4,-h1,e4Cm)
      call pol_real(p5,h1,e5m)
      h1=1
      call VKlSt(p1,mt,p3,h1,e1Cp)
      call UbKlSt(p2,mt,p3,h1,e2p)
      call ubarspinor0(p3,h1,e3p)
      call uspinor0(p3,-h1,e3Cp)
      call ubarspinor0(p4,h1,e4p)
      call uspinor0(p4,-h1,e4Cp)
      call pol_real(p5,h1,e5p)
      do h3=1,2
      if (h3 == 1) then
         e3(:)=e3m(:)
         e3C(:)=e3Cm(:)
      elseif (h3 == 2) then
         e3(:)=e3p(:)
         e3C(:)=e3Cp(:)
      endif
      h4=3-h3
      if (h4 == 1) then
         e4(:)=e4m(:)
         e4C(:)=e4Cm(:)
      elseif (h4 == 2) then
         e4(:)=e4p(:)
         e4C(:)=e4Cp(:)
      endif
      do h5=1,2
      if (h5 == 1) e5(:)=e5m(:)
      if (h5 == 2) e5(:)=e5p(:)
      do h1=1,2
      if (h1 == 1) then
          e1C(:)=e1Cm(:)
      elseif (h1 == 2) then
          e1C(:)=e1Cp(:)
      endif
      do h2=1,2
      if (h2 == 1) then
          e2(:)=e2m(:)
      elseif (h2 == 2) then
          e2(:)=e2p(:)
      endif
      T34x21(h1,h2,h3,h4,h5)= + xn**(-1)*two*s345**(-1) * (  - 
     &    spstrng0(e2,e3C)*spstrng3(e4,e5,p45,p2345,e1C)*
     &    s2345mmtsq**(-1)*s45**(-1) + spstrng0(e2,e4C)*spstrng3(e3,e5,
     &    p35,p2345,e1C)*s2345mmtsq**(-1)*s35**(-1) + spstrng0(e3,e1C)*
     &    spstrng3(e2,p1345,p45,e5,e4C)*s1345mmtsq**(-1)*s45**(-1) - 
     &    spstrng0(e4,e1C)*spstrng3(e2,p1345,p35,e5,e3C)*
     &    s1345mmtsq**(-1)*s35**(-1) + spstrng1(e2,p1345,e3C)*spstrng2(
     &    e4,e5,p45,e1C)*s1345mmtsq**(-1)*s45**(-1) - spstrng1(e2,p1345
     &    ,e4C)*spstrng2(e3,e5,p35,e1C)*s1345mmtsq**(-1)*s35**(-1) - 
     &    spstrng1(e3,p2345,e1C)*spstrng2(e2,p45,e5,e4C)*
     &    s2345mmtsq**(-1)*s45**(-1) + spstrng1(e4,p2345,e1C)*spstrng2(
     &    e2,p35,e5,e3C)*s2345mmtsq**(-1)*s35**(-1) )
      T34x21(h1,h2,h3,h4,h5) = T34x21(h1,h2,h3,h4,h5) + xn**(-1)*two*mt
     & *s345**(-1) * (  - spstrng0(e2,e3C)*spstrng2(e4,e5,p45,e1C)*
     &    s1345mmtsq**(-1)*s45**(-1) - spstrng0(e2,e3C)*spstrng2(e4,e5,
     &    p45,e1C)*s2345mmtsq**(-1)*s45**(-1) + spstrng0(e2,e4C)*
     &    spstrng2(e3,e5,p35,e1C)*s1345mmtsq**(-1)*s35**(-1) + 
     &    spstrng0(e2,e4C)*spstrng2(e3,e5,p35,e1C)*s2345mmtsq**(-1)*
     &    s35**(-1) - spstrng0(e3,e1C)*spstrng2(e2,p45,e5,e4C)*
     &    s1345mmtsq**(-1)*s45**(-1) - spstrng0(e3,e1C)*spstrng2(e2,p45
     &    ,e5,e4C)*s2345mmtsq**(-1)*s45**(-1) + spstrng0(e4,e1C)*
     &    spstrng2(e2,p35,e5,e3C)*s1345mmtsq**(-1)*s35**(-1) + 
     &    spstrng0(e4,e1C)*spstrng2(e2,p35,e5,e3C)*s2345mmtsq**(-1)*
     &    s35**(-1) )

      T21x34(h1,h2,h3,h4,h5)= + s34**(-1)*xn**(-1)*two * ( spstrng0(e2,
     &    e3C)*spstrng3(e4,p234,e5,p2345,e1C)*s234mmtsq**(-1)*
     &    s2345mmtsq**(-1) - spstrng0(e2,e3C)*spstrng3(e4,p234,p15,e5,
     &    e1C)*s234mmtsq**(-1)*s15mmtsq**(-1) + spstrng0(e2,e4C)*
     &    spstrng3(e3,p234,e5,p2345,e1C)*s234mmtsq**(-1)*
     &    s2345mmtsq**(-1) - spstrng0(e2,e4C)*spstrng3(e3,p234,p15,e5,
     &    e1C)*s234mmtsq**(-1)*s15mmtsq**(-1) - spstrng0(e3,e1C)*
     &    spstrng3(e2,e5,p25,p134,e4C)*s134mmtsq**(-1)*s25mmtsq**(-1)
     &     + spstrng0(e3,e1C)*spstrng3(e2,p1345,e5,p134,e4C)*
     &    s134mmtsq**(-1)*s1345mmtsq**(-1) - spstrng0(e4,e1C)*spstrng3(
     &    e2,e5,p25,p134,e3C)*s134mmtsq**(-1)*s25mmtsq**(-1) + 
     &    spstrng0(e4,e1C)*spstrng3(e2,p1345,e5,p134,e3C)*
     &    s134mmtsq**(-1)*s1345mmtsq**(-1) + spstrng1(e2,p1345,e3C)*
     &    spstrng2(e4,p15,e5,e1C)*s1345mmtsq**(-1)*s15mmtsq**(-1) + 
     &    spstrng1(e2,p1345,e4C)*spstrng2(e3,p15,e5,e1C)*
     &    s1345mmtsq**(-1)*s15mmtsq**(-1) )
      T21x34(h1,h2,h3,h4,h5) = T21x34(h1,h2,h3,h4,h5) + s34**(-1)*
     & xn**(-1)*two * ( spstrng1(e3,p2345,e1C)*spstrng2(e2,e5,p25,e4C)*
     &    s2345mmtsq**(-1)*s25mmtsq**(-1) + spstrng1(e4,p2345,e1C)*
     &    spstrng2(e2,e5,p25,e3C)*s2345mmtsq**(-1)*s25mmtsq**(-1) )
      T21x34(h1,h2,h3,h4,h5) = T21x34(h1,h2,h3,h4,h5) + s34**(-1)*
     & xn**(-1)*two*mt * ( spstrng0(e2,e3C)*spstrng2(e4,e5,p2345,e1C)*
     &    s234mmtsq**(-1)*s2345mmtsq**(-1) - spstrng0(e2,e3C)*spstrng2(
     &    e4,p15,e5,e1C)*s234mmtsq**(-1)*s15mmtsq**(-1) - spstrng0(e2,
     &    e3C)*spstrng2(e4,p15,e5,e1C)*s1345mmtsq**(-1)*s15mmtsq**(-1)
     &     + spstrng0(e2,e3C)*spstrng2(e4,p234,e5,e1C)*s234mmtsq**(-1)*
     &    s2345mmtsq**(-1) + spstrng0(e2,e3C)*spstrng2(e4,p234,e5,e1C)*
     &    s234mmtsq**(-1)*s15mmtsq**(-1) + spstrng0(e2,e4C)*spstrng2(e3
     &    ,e5,p2345,e1C)*s234mmtsq**(-1)*s2345mmtsq**(-1) - spstrng0(e2
     &    ,e4C)*spstrng2(e3,p15,e5,e1C)*s234mmtsq**(-1)*s15mmtsq**(-1)
     &     - spstrng0(e2,e4C)*spstrng2(e3,p15,e5,e1C)*s1345mmtsq**(-1)*
     &    s15mmtsq**(-1) + spstrng0(e2,e4C)*spstrng2(e3,p234,e5,e1C)*
     &    s234mmtsq**(-1)*s2345mmtsq**(-1) + spstrng0(e2,e4C)*spstrng2(
     &    e3,p234,e5,e1C)*s234mmtsq**(-1)*s15mmtsq**(-1) + spstrng0(e3,
     &    e1C)*spstrng2(e2,e5,p25,e4C)*s134mmtsq**(-1)*s25mmtsq**(-1)
     &     + spstrng0(e3,e1C)*spstrng2(e2,e5,p25,e4C)*s2345mmtsq**(-1)*
     &    s25mmtsq**(-1) )
      T21x34(h1,h2,h3,h4,h5) = T21x34(h1,h2,h3,h4,h5) + s34**(-1)*
     & xn**(-1)*two*mt * (  - spstrng0(e3,e1C)*spstrng2(e2,e5,p134,e4C)
     &    *s134mmtsq**(-1)*s1345mmtsq**(-1) - spstrng0(e3,e1C)*
     &    spstrng2(e2,e5,p134,e4C)*s134mmtsq**(-1)*s25mmtsq**(-1) - 
     &    spstrng0(e3,e1C)*spstrng2(e2,p1345,e5,e4C)*s134mmtsq**(-1)*
     &    s1345mmtsq**(-1) + spstrng0(e4,e1C)*spstrng2(e2,e5,p25,e3C)*
     &    s134mmtsq**(-1)*s25mmtsq**(-1) + spstrng0(e4,e1C)*spstrng2(e2
     &    ,e5,p25,e3C)*s2345mmtsq**(-1)*s25mmtsq**(-1) - spstrng0(e4,
     &    e1C)*spstrng2(e2,e5,p134,e3C)*s134mmtsq**(-1)*
     &    s1345mmtsq**(-1) - spstrng0(e4,e1C)*spstrng2(e2,e5,p134,e3C)*
     &    s134mmtsq**(-1)*s25mmtsq**(-1) - spstrng0(e4,e1C)*spstrng2(e2
     &    ,p1345,e5,e3C)*s134mmtsq**(-1)*s1345mmtsq**(-1) + spstrng1(e2
     &    ,e5,e3C)*spstrng1(e4,p2345,e1C)*s2345mmtsq**(-1)*
     &    s25mmtsq**(-1) + spstrng1(e2,e5,e4C)*spstrng1(e3,p2345,e1C)*
     &    s2345mmtsq**(-1)*s25mmtsq**(-1) - spstrng1(e2,p1345,e3C)*
     &    spstrng1(e4,e5,e1C)*s1345mmtsq**(-1)*s15mmtsq**(-1) )
      T21x34(h1,h2,h3,h4,h5) = T21x34(h1,h2,h3,h4,h5) + s34**(-1)*
     & xn**(-1)*two*mt * (  - spstrng1(e2,p1345,e4C)*spstrng1(e3,e5,e1C
     &    )*s1345mmtsq**(-1)*s15mmtsq**(-1) )
      T21x34(h1,h2,h3,h4,h5) = T21x34(h1,h2,h3,h4,h5) + s34**(-1)*
     & xn**(-1)*two*mt**2 * ( spstrng0(e2,e3C)*spstrng1(e4,e5,e1C)*
     &    s234mmtsq**(-1)*s2345mmtsq**(-1) + spstrng0(e2,e3C)*spstrng1(
     &    e4,e5,e1C)*s234mmtsq**(-1)*s15mmtsq**(-1) + spstrng0(e2,e3C)*
     &    spstrng1(e4,e5,e1C)*s1345mmtsq**(-1)*s15mmtsq**(-1) + 
     &    spstrng0(e2,e4C)*spstrng1(e3,e5,e1C)*s234mmtsq**(-1)*
     &    s2345mmtsq**(-1) + spstrng0(e2,e4C)*spstrng1(e3,e5,e1C)*
     &    s234mmtsq**(-1)*s15mmtsq**(-1) + spstrng0(e2,e4C)*spstrng1(e3
     &    ,e5,e1C)*s1345mmtsq**(-1)*s15mmtsq**(-1) + spstrng0(e3,e1C)*
     &    spstrng1(e2,e5,e4C)*s134mmtsq**(-1)*s1345mmtsq**(-1) + 
     &    spstrng0(e3,e1C)*spstrng1(e2,e5,e4C)*s134mmtsq**(-1)*
     &    s25mmtsq**(-1) + spstrng0(e3,e1C)*spstrng1(e2,e5,e4C)*
     &    s2345mmtsq**(-1)*s25mmtsq**(-1) + spstrng0(e4,e1C)*spstrng1(
     &    e2,e5,e3C)*s134mmtsq**(-1)*s1345mmtsq**(-1) + spstrng0(e4,e1C
     &    )*spstrng1(e2,e5,e3C)*s134mmtsq**(-1)*s25mmtsq**(-1) + 
     &    spstrng0(e4,e1C)*spstrng1(e2,e5,e3C)*s2345mmtsq**(-1)*
     &    s25mmtsq**(-1) )

      T24x31(h1,h2,h3,h4,h5)= + s34**(-1)*two*s345**(-1) * ( cdot(e5,
     &    p34)*spstrng0(e2,e3C)*spstrng1(e4,p2345,e1C)*s2345mmtsq**(-1)
     &     + cdot(e5,p34)*spstrng0(e2,e4C)*spstrng1(e3,p2345,e1C)*
     &    s2345mmtsq**(-1) - cdot(e5,p34)*spstrng0(e3,e1C)*spstrng1(e2,
     &    p1345,e4C)*s1345mmtsq**(-1) - cdot(e5,p34)*spstrng0(e4,e1C)*
     &    spstrng1(e2,p1345,e3C)*s1345mmtsq**(-1) + cdot(e5,p345)*
     &    spstrng0(e2,e3C)*spstrng1(e4,p2345,e1C)*s2345mmtsq**(-1) + 
     &    cdot(e5,p345)*spstrng0(e2,e4C)*spstrng1(e3,p2345,e1C)*
     &    s2345mmtsq**(-1) - cdot(e5,p345)*spstrng0(e3,e1C)*spstrng1(e2
     &    ,p1345,e4C)*s1345mmtsq**(-1) - cdot(e5,p345)*spstrng0(e4,e1C)
     &    *spstrng1(e2,p1345,e3C)*s1345mmtsq**(-1) - spstrng1(e3,e5,e4C
     &    )*spstrng2(e2,p34,p2345,e1C)*s2345mmtsq**(-1) + spstrng1(e3,
     &    e5,e4C)*spstrng2(e2,p1345,p34,e1C)*s1345mmtsq**(-1) + 1.D0/2.D
     &    0*spstrng1(e3,p34,e4C)*spstrng2(e2,e5,p2345,e1C)*
     &    s2345mmtsq**(-1) - 1.D0/2.D0*spstrng1(e3,p34,e4C)*spstrng2(e2
     &    ,p1345,e5,e1C)*s1345mmtsq**(-1) )
      T24x31(h1,h2,h3,h4,h5) = T24x31(h1,h2,h3,h4,h5) + s34**(-1)*two*
     & s345**(-1) * (  - spstrng1(e3,p345,e4C)*spstrng2(e2,e5,p2345,e1C
     &    )*s2345mmtsq**(-1) + spstrng1(e3,p345,e4C)*spstrng2(e2,p1345,
     &    e5,e1C)*s1345mmtsq**(-1) )
      T24x31(h1,h2,h3,h4,h5) = T24x31(h1,h2,h3,h4,h5) + s34**(-1)*two
     &  * ( spstrng0(e3,e1C)*spstrng3(e2,e5,p25,p134,e4C)*
     &    s134mmtsq**(-1)*s25mmtsq**(-1) - spstrng0(e3,e1C)*spstrng3(e2
     &    ,p1345,e5,p134,e4C)*s134mmtsq**(-1)*s1345mmtsq**(-1) + 
     &    spstrng0(e4,e1C)*spstrng3(e2,e5,p25,p134,e3C)*s134mmtsq**(-1)
     &    *s25mmtsq**(-1) - spstrng0(e4,e1C)*spstrng3(e2,p1345,e5,p134,
     &    e3C)*s134mmtsq**(-1)*s1345mmtsq**(-1) - spstrng1(e3,p2345,e1C
     &    )*spstrng2(e2,e5,p25,e4C)*s2345mmtsq**(-1)*s25mmtsq**(-1) - 
     &    spstrng1(e4,p2345,e1C)*spstrng2(e2,e5,p25,e3C)*
     &    s2345mmtsq**(-1)*s25mmtsq**(-1) )
      T24x31(h1,h2,h3,h4,h5) = T24x31(h1,h2,h3,h4,h5) + s34**(-1)*two*
     & mt*s345**(-1) * ( cdot(e5,p34)*spstrng0(e2,e3C)*spstrng0(e4,e1C)
     &    *s1345mmtsq**(-1) + cdot(e5,p34)*spstrng0(e2,e3C)*spstrng0(e4
     &    ,e1C)*s2345mmtsq**(-1) + cdot(e5,p34)*spstrng0(e2,e4C)*
     &    spstrng0(e3,e1C)*s1345mmtsq**(-1) + cdot(e5,p34)*spstrng0(e2,
     &    e4C)*spstrng0(e3,e1C)*s2345mmtsq**(-1) + cdot(e5,p345)*
     &    spstrng0(e2,e3C)*spstrng0(e4,e1C)*s1345mmtsq**(-1) + cdot(e5,
     &    p345)*spstrng0(e2,e3C)*spstrng0(e4,e1C)*s2345mmtsq**(-1) + 
     &    cdot(e5,p345)*spstrng0(e2,e4C)*spstrng0(e3,e1C)*
     &    s1345mmtsq**(-1) + cdot(e5,p345)*spstrng0(e2,e4C)*spstrng0(e3
     &    ,e1C)*s2345mmtsq**(-1) + 1.D0/2.D0*spstrng1(e2,e5,e1C)*
     &    spstrng1(e3,p34,e4C)*s1345mmtsq**(-1) + 1.D0/2.D0*spstrng1(e2
     &    ,e5,e1C)*spstrng1(e3,p34,e4C)*s2345mmtsq**(-1) - spstrng1(e2,
     &    e5,e1C)*spstrng1(e3,p345,e4C)*s1345mmtsq**(-1) - spstrng1(e2,
     &    e5,e1C)*spstrng1(e3,p345,e4C)*s2345mmtsq**(-1) - spstrng1(e2,
     &    p34,e1C)*spstrng1(e3,e5,e4C)*s1345mmtsq**(-1) )
      T24x31(h1,h2,h3,h4,h5) = T24x31(h1,h2,h3,h4,h5) + s34**(-1)*two*
     & mt*s345**(-1) * (  - spstrng1(e2,p34,e1C)*spstrng1(e3,e5,e4C)*
     &    s2345mmtsq**(-1) )
      T24x31(h1,h2,h3,h4,h5) = T24x31(h1,h2,h3,h4,h5) + s34**(-1)*two*
     & mt * (  - spstrng0(e3,e1C)*spstrng2(e2,e5,p25,e4C)*
     &    s134mmtsq**(-1)*s25mmtsq**(-1) - spstrng0(e3,e1C)*spstrng2(e2
     &    ,e5,p25,e4C)*s2345mmtsq**(-1)*s25mmtsq**(-1) + spstrng0(e3,
     &    e1C)*spstrng2(e2,e5,p134,e4C)*s134mmtsq**(-1)*
     &    s1345mmtsq**(-1) + spstrng0(e3,e1C)*spstrng2(e2,e5,p134,e4C)*
     &    s134mmtsq**(-1)*s25mmtsq**(-1) + spstrng0(e3,e1C)*spstrng2(e2
     &    ,p1345,e5,e4C)*s134mmtsq**(-1)*s1345mmtsq**(-1) - spstrng0(e4
     &    ,e1C)*spstrng2(e2,e5,p25,e3C)*s134mmtsq**(-1)*s25mmtsq**(-1)
     &     - spstrng0(e4,e1C)*spstrng2(e2,e5,p25,e3C)*s2345mmtsq**(-1)*
     &    s25mmtsq**(-1) + spstrng0(e4,e1C)*spstrng2(e2,e5,p134,e3C)*
     &    s134mmtsq**(-1)*s1345mmtsq**(-1) + spstrng0(e4,e1C)*spstrng2(
     &    e2,e5,p134,e3C)*s134mmtsq**(-1)*s25mmtsq**(-1) + spstrng0(e4,
     &    e1C)*spstrng2(e2,p1345,e5,e3C)*s134mmtsq**(-1)*
     &    s1345mmtsq**(-1) - spstrng1(e2,e5,e3C)*spstrng1(e4,p2345,e1C)
     &    *s2345mmtsq**(-1)*s25mmtsq**(-1) )
      T24x31(h1,h2,h3,h4,h5) = T24x31(h1,h2,h3,h4,h5) + s34**(-1)*two*
     & mt * (  - spstrng1(e2,e5,e4C)*spstrng1(e3,p2345,e1C)*
     &    s2345mmtsq**(-1)*s25mmtsq**(-1) )
      T24x31(h1,h2,h3,h4,h5) = T24x31(h1,h2,h3,h4,h5) + s34**(-1)*two*
     & mt**2 * (  - spstrng0(e3,e1C)*spstrng1(e2,e5,e4C)*
     &    s134mmtsq**(-1)*s1345mmtsq**(-1) - spstrng0(e3,e1C)*spstrng1(
     &    e2,e5,e4C)*s134mmtsq**(-1)*s25mmtsq**(-1) - spstrng0(e3,e1C)*
     &    spstrng1(e2,e5,e4C)*s2345mmtsq**(-1)*s25mmtsq**(-1) - 
     &    spstrng0(e4,e1C)*spstrng1(e2,e5,e3C)*s134mmtsq**(-1)*
     &    s1345mmtsq**(-1) - spstrng0(e4,e1C)*spstrng1(e2,e5,e3C)*
     &    s134mmtsq**(-1)*s25mmtsq**(-1) - spstrng0(e4,e1C)*spstrng1(e2
     &    ,e5,e3C)*s2345mmtsq**(-1)*s25mmtsq**(-1) )
      T24x31(h1,h2,h3,h4,h5) = T24x31(h1,h2,h3,h4,h5) + two*s345**(-1)
     &  * ( spstrng0(e2,e3C)*spstrng3(e4,e5,p45,p2345,e1C)*
     &    s2345mmtsq**(-1)*s45**(-1) - spstrng0(e3,e1C)*spstrng3(e2,
     &    p1345,p45,e5,e4C)*s1345mmtsq**(-1)*s45**(-1) - spstrng1(e2,
     &    p1345,e3C)*spstrng2(e4,e5,p45,e1C)*s1345mmtsq**(-1)*s45**(-1)
     &     + spstrng1(e3,p2345,e1C)*spstrng2(e2,p45,e5,e4C)*
     &    s2345mmtsq**(-1)*s45**(-1) )
      T24x31(h1,h2,h3,h4,h5) = T24x31(h1,h2,h3,h4,h5) + two*mt*
     & s345**(-1) * ( spstrng0(e2,e3C)*spstrng2(e4,e5,p45,e1C)*
     &    s1345mmtsq**(-1)*s45**(-1) + spstrng0(e2,e3C)*spstrng2(e4,e5,
     &    p45,e1C)*s2345mmtsq**(-1)*s45**(-1) + spstrng0(e3,e1C)*
     &    spstrng2(e2,p45,e5,e4C)*s1345mmtsq**(-1)*s45**(-1) + 
     &    spstrng0(e3,e1C)*spstrng2(e2,p45,e5,e4C)*s2345mmtsq**(-1)*
     &    s45**(-1) )

      T31x24(h1,h2,h3,h4,h5)= + s34**(-1)*two*s345**(-1) * (  - cdot(e5
     &    ,p34)*spstrng0(e2,e3C)*spstrng1(e4,p2345,e1C)*
     &    s2345mmtsq**(-1) - cdot(e5,p34)*spstrng0(e2,e4C)*spstrng1(e3,
     &    p2345,e1C)*s2345mmtsq**(-1) + cdot(e5,p34)*spstrng0(e3,e1C)*
     &    spstrng1(e2,p1345,e4C)*s1345mmtsq**(-1) + cdot(e5,p34)*
     &    spstrng0(e4,e1C)*spstrng1(e2,p1345,e3C)*s1345mmtsq**(-1) - 
     &    cdot(e5,p345)*spstrng0(e2,e3C)*spstrng1(e4,p2345,e1C)*
     &    s2345mmtsq**(-1) - cdot(e5,p345)*spstrng0(e2,e4C)*spstrng1(e3
     &    ,p2345,e1C)*s2345mmtsq**(-1) + cdot(e5,p345)*spstrng0(e3,e1C)
     &    *spstrng1(e2,p1345,e4C)*s1345mmtsq**(-1) + cdot(e5,p345)*
     &    spstrng0(e4,e1C)*spstrng1(e2,p1345,e3C)*s1345mmtsq**(-1) + 
     &    spstrng1(e3,e5,e4C)*spstrng2(e2,p34,p2345,e1C)*
     &    s2345mmtsq**(-1) - spstrng1(e3,e5,e4C)*spstrng2(e2,p1345,p34,
     &    e1C)*s1345mmtsq**(-1) - 1.D0/2.D0*spstrng1(e3,p34,e4C)*
     &    spstrng2(e2,e5,p2345,e1C)*s2345mmtsq**(-1) + 1.D0/2.D0*
     &    spstrng1(e3,p34,e4C)*spstrng2(e2,p1345,e5,e1C)*
     &    s1345mmtsq**(-1) )
      T31x24(h1,h2,h3,h4,h5) = T31x24(h1,h2,h3,h4,h5) + s34**(-1)*two*
     & s345**(-1) * ( spstrng1(e3,p345,e4C)*spstrng2(e2,e5,p2345,e1C)*
     &    s2345mmtsq**(-1) - spstrng1(e3,p345,e4C)*spstrng2(e2,p1345,e5
     &    ,e1C)*s1345mmtsq**(-1) )
      T31x24(h1,h2,h3,h4,h5) = T31x24(h1,h2,h3,h4,h5) + s34**(-1)*two
     &  * (  - spstrng0(e2,e3C)*spstrng3(e4,p234,e5,p2345,e1C)*
     &    s234mmtsq**(-1)*s2345mmtsq**(-1) + spstrng0(e2,e3C)*spstrng3(
     &    e4,p234,p15,e5,e1C)*s234mmtsq**(-1)*s15mmtsq**(-1) - 
     &    spstrng0(e2,e4C)*spstrng3(e3,p234,e5,p2345,e1C)*
     &    s234mmtsq**(-1)*s2345mmtsq**(-1) + spstrng0(e2,e4C)*spstrng3(
     &    e3,p234,p15,e5,e1C)*s234mmtsq**(-1)*s15mmtsq**(-1) - 
     &    spstrng1(e2,p1345,e3C)*spstrng2(e4,p15,e5,e1C)*
     &    s1345mmtsq**(-1)*s15mmtsq**(-1) - spstrng1(e2,p1345,e4C)*
     &    spstrng2(e3,p15,e5,e1C)*s1345mmtsq**(-1)*s15mmtsq**(-1) )
      T31x24(h1,h2,h3,h4,h5) = T31x24(h1,h2,h3,h4,h5) + s34**(-1)*two*
     & mt*s345**(-1) * (  - cdot(e5,p34)*spstrng0(e2,e3C)*spstrng0(e4,
     &    e1C)*s1345mmtsq**(-1) - cdot(e5,p34)*spstrng0(e2,e3C)*
     &    spstrng0(e4,e1C)*s2345mmtsq**(-1) - cdot(e5,p34)*spstrng0(e2,
     &    e4C)*spstrng0(e3,e1C)*s1345mmtsq**(-1) - cdot(e5,p34)*
     &    spstrng0(e2,e4C)*spstrng0(e3,e1C)*s2345mmtsq**(-1) - cdot(e5,
     &    p345)*spstrng0(e2,e3C)*spstrng0(e4,e1C)*s1345mmtsq**(-1) - 
     &    cdot(e5,p345)*spstrng0(e2,e3C)*spstrng0(e4,e1C)*
     &    s2345mmtsq**(-1) - cdot(e5,p345)*spstrng0(e2,e4C)*spstrng0(e3
     &    ,e1C)*s1345mmtsq**(-1) - cdot(e5,p345)*spstrng0(e2,e4C)*
     &    spstrng0(e3,e1C)*s2345mmtsq**(-1) - 1.D0/2.D0*spstrng1(e2,e5,
     &    e1C)*spstrng1(e3,p34,e4C)*s1345mmtsq**(-1) - 1.D0/2.D0*
     &    spstrng1(e2,e5,e1C)*spstrng1(e3,p34,e4C)*s2345mmtsq**(-1) + 
     &    spstrng1(e2,e5,e1C)*spstrng1(e3,p345,e4C)*s1345mmtsq**(-1) + 
     &    spstrng1(e2,e5,e1C)*spstrng1(e3,p345,e4C)*s2345mmtsq**(-1) + 
     &    spstrng1(e2,p34,e1C)*spstrng1(e3,e5,e4C)*s1345mmtsq**(-1) )
      T31x24(h1,h2,h3,h4,h5) = T31x24(h1,h2,h3,h4,h5) + s34**(-1)*two*
     & mt*s345**(-1) * ( spstrng1(e2,p34,e1C)*spstrng1(e3,e5,e4C)*
     &    s2345mmtsq**(-1) )
      T31x24(h1,h2,h3,h4,h5) = T31x24(h1,h2,h3,h4,h5) + s34**(-1)*two*
     & mt * (  - spstrng0(e2,e3C)*spstrng2(e4,e5,p2345,e1C)*
     &    s234mmtsq**(-1)*s2345mmtsq**(-1) + spstrng0(e2,e3C)*spstrng2(
     &    e4,p15,e5,e1C)*s234mmtsq**(-1)*s15mmtsq**(-1) + spstrng0(e2,
     &    e3C)*spstrng2(e4,p15,e5,e1C)*s1345mmtsq**(-1)*s15mmtsq**(-1)
     &     - spstrng0(e2,e3C)*spstrng2(e4,p234,e5,e1C)*s234mmtsq**(-1)*
     &    s2345mmtsq**(-1) - spstrng0(e2,e3C)*spstrng2(e4,p234,e5,e1C)*
     &    s234mmtsq**(-1)*s15mmtsq**(-1) - spstrng0(e2,e4C)*spstrng2(e3
     &    ,e5,p2345,e1C)*s234mmtsq**(-1)*s2345mmtsq**(-1) + spstrng0(e2
     &    ,e4C)*spstrng2(e3,p15,e5,e1C)*s234mmtsq**(-1)*s15mmtsq**(-1)
     &     + spstrng0(e2,e4C)*spstrng2(e3,p15,e5,e1C)*s1345mmtsq**(-1)*
     &    s15mmtsq**(-1) - spstrng0(e2,e4C)*spstrng2(e3,p234,e5,e1C)*
     &    s234mmtsq**(-1)*s2345mmtsq**(-1) - spstrng0(e2,e4C)*spstrng2(
     &    e3,p234,e5,e1C)*s234mmtsq**(-1)*s15mmtsq**(-1) + spstrng1(e2,
     &    p1345,e3C)*spstrng1(e4,e5,e1C)*s1345mmtsq**(-1)*
     &    s15mmtsq**(-1) )
      T31x24(h1,h2,h3,h4,h5) = T31x24(h1,h2,h3,h4,h5) + s34**(-1)*two*
     & mt * ( spstrng1(e2,p1345,e4C)*spstrng1(e3,e5,e1C)*
     &    s1345mmtsq**(-1)*s15mmtsq**(-1) )
      T31x24(h1,h2,h3,h4,h5) = T31x24(h1,h2,h3,h4,h5) + s34**(-1)*two*
     & mt**2 * (  - spstrng0(e2,e3C)*spstrng1(e4,e5,e1C)*
     &    s234mmtsq**(-1)*s2345mmtsq**(-1) - spstrng0(e2,e3C)*spstrng1(
     &    e4,e5,e1C)*s234mmtsq**(-1)*s15mmtsq**(-1) - spstrng0(e2,e3C)*
     &    spstrng1(e4,e5,e1C)*s1345mmtsq**(-1)*s15mmtsq**(-1) - 
     &    spstrng0(e2,e4C)*spstrng1(e3,e5,e1C)*s234mmtsq**(-1)*
     &    s2345mmtsq**(-1) - spstrng0(e2,e4C)*spstrng1(e3,e5,e1C)*
     &    s234mmtsq**(-1)*s15mmtsq**(-1) - spstrng0(e2,e4C)*spstrng1(e3
     &    ,e5,e1C)*s1345mmtsq**(-1)*s15mmtsq**(-1) )
      T31x24(h1,h2,h3,h4,h5) = T31x24(h1,h2,h3,h4,h5) + two*s345**(-1)
     &  * (  - spstrng0(e2,e4C)*spstrng3(e3,e5,p35,p2345,e1C)*
     &    s2345mmtsq**(-1)*s35**(-1) + spstrng0(e4,e1C)*spstrng3(e2,
     &    p1345,p35,e5,e3C)*s1345mmtsq**(-1)*s35**(-1) + spstrng1(e2,
     &    p1345,e4C)*spstrng2(e3,e5,p35,e1C)*s1345mmtsq**(-1)*s35**(-1)
     &     - spstrng1(e4,p2345,e1C)*spstrng2(e2,p35,e5,e3C)*
     &    s2345mmtsq**(-1)*s35**(-1) )
      T31x24(h1,h2,h3,h4,h5) = T31x24(h1,h2,h3,h4,h5) + two*mt*
     & s345**(-1) * (  - spstrng0(e2,e4C)*spstrng2(e3,e5,p35,e1C)*
     &    s1345mmtsq**(-1)*s35**(-1) - spstrng0(e2,e4C)*spstrng2(e3,e5,
     &    p35,e1C)*s2345mmtsq**(-1)*s35**(-1) - spstrng0(e4,e1C)*
     &    spstrng2(e2,p35,e5,e3C)*s1345mmtsq**(-1)*s35**(-1) - 
     &    spstrng0(e4,e1C)*spstrng2(e2,p35,e5,e3C)*s2345mmtsq**(-1)*
     &    s35**(-1) )

      enddo
      enddo
      enddo
      enddo
      ampsq=0d0
      do h1=1,2
      do h2=1,2
      do h3=1,2
      h4=3-h3
      do h5=1,2
      m(1)=T34x21(h1,h2,h3,h4,h5)
      m(2)=T21x34(h1,h2,h3,h4,h5)
      m(3)=T24x31(h1,h2,h3,h4,h5)
      m(4)=T31x24(h1,h2,h3,h4,h5)
C     Factor of V removed from color sum
      ampsq=ampsq+xn*real(
     & +m(1)*conjg(m(1))+m(2)*conjg(m(2))
     & +m(3)*conjg(m(3))+m(4)*conjg(m(4)))
     & +real((m(1)+m(2))*conjg(m(3)+m(4))
     &      +(m(3)+m(4))*conjg(m(1)+m(2)))
      enddo
      enddo
      enddo
      enddo
      return
      end
