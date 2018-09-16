      subroutine ttqqZampsq(p,j34,i1,i2,i3,i4,i5,i6,ampsq)
      implicit none
      include 'types.f'
      
c      t(i2)+t~(i1)+q(i3)+q~(i4)+e-(i5)+e+(i6)
c      j34 is flavour of 34 line
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      real(dp):: p(mxpart,4),s134mmtsq,s234mmtsq,s124,s123,
     & s134,s234,s12,s34,s56,ampsq
      complex(dp):: spstrng0,spstrng1,
     & e1C(4),e2(4),p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),
     & e1Cp(4),e1Cm(4),e2p(4),e2m(4),
     & e3m(4),e3Cm(4),e3p(4),e3Cp(4),e4m(4),e4Cm(4),e4p(4),e4Cp(4),
     & e5m(4),e5Cp(4),e6Cm(4),e6p(4),
     & p123(4),p124(4),p134(4),p234(4),
     & p12(4),p34(4),p56(4),
     & amp(2,2,2,2),propL(2,2),propR(2,2),prop
      integer:: i1,i2,i3,i4,q5,q6,i5,i6,h1,h2,h34,h56,jswap,j,j34
      p1(:)=cplx2(p(i1,:))
      p2(:)=cplx2(p(i2,:))
      p3(:)=cplx2(p(i3,:))
      p4(:)=cplx2(p(i4,:))
      p123(:)=cplx2(p(i1,:)+p(i2,:)+p(i3,:))
      p124(:)=cplx2(p(i1,:)+p(i2,:)+p(i4,:))
      p134(:)=cplx2(p(i1,:)+p(i3,:)+p(i4,:))
      p234(:)=cplx2(p(i2,:)+p(i3,:)+p(i4,:))
      p12(:)=cplx2(p(i1,:)+p(i2,:))
      p34(:)=cplx2(p(i3,:)+p(i4,:))
      p56(:)=cplx2(p(i5,:)+p(i6,:))
      s12=real(p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2)
      s34=real(p34(4)**2-p34(1)**2-p34(2)**2-p34(3)**2)
      s56=real(p56(4)**2-p56(1)**2-p56(2)**2-p56(3)**2)
      prop=s56/cplx2((s56-zmass**2),zmass*zwidth)

      do j=1,2
C---first index is flavour, second is helicity of 56 line
      propL(j,1)=(Q(j)*q1+L(j)*l1*prop)
      propL(j,2)=(Q(j)*q1+L(j)*r1*prop)
      propR(j,1)=(Q(j)*q1+R(j)*l1*prop)
      propR(j,2)=(Q(j)*q1+R(j)*r1*prop)
      enddo
      s134=real(p134(4)**2-p134(1)**2-p134(2)**2-p134(3)**2)
      s234=real(p234(4)**2-p234(1)**2-p234(2)**2-p234(3)**2)
      s123=real(p123(4)**2-p123(1)**2-p123(2)**2-p123(3)**2)
      s124=real(p124(4)**2-p124(1)**2-p124(2)**2-p124(3)**2)
      s134mmtsq=s134-mt**2
      s234mmtsq=s234-mt**2
      h1=1
      call UbKlSt(p2,mt,p3,h1,e2p)
      call VKlSt(p1,mt,p3,h1,e1Cp)
      call ubarspinor0(p3,h1,e3p)
      call ubarspinor0(p4,h1,e4p)
      call uspinor0(p3,h1,e3Cp)
      call uspinor0(p4,h1,e4Cp)
      h1=-1
      call UbKlst(p2,mt,p3,h1,e2m)
      call VklSt(p1,mt,p3,h1,e1Cm)
      call ubarspinor0(p3,h1,e3m)
      call ubarspinor0(p4,h1,e4m)
      call uspinor0(p4,h1,e4Cm)
      call uspinor0(p3,h1,e3Cm)
      
      do jswap=1,2
      if (jswap == 1) then
      q5=i5
      q6=i6
      elseif (jswap == 2) then
      q5=i6
      q6=i5
      endif
      p5(:)=cplx2(p(q5,:))
      p6(:)=cplx2(p(q6,:))
      
      h1=1
      call ubarspinor0(p6,h1,e6p)
      call uspinor0(p5,h1,e5Cp)
      h1=-1
      call ubarspinor0(p5,h1,e5m)
      call uspinor0(p6,h1,e6Cm)
      do h1=1,2
      if (h1 == 1) e1C(:)=e1Cm(:)
      if (h1 == 2) e1C(:)=e1Cp(:)
      do h2=1,2
      if (h2 == 1) e2(:)=e2m(:)
      if (h2 == 2) e2(:)=e2p(:)
      amp(h1,h2,1,jswap)= + s34**(-1)*s56**(-1) * (  - 4.D0*propL(2,
     &    jswap)*spstrng0(e4p,e1C)*spstrng0(e2,e5Cp)*spstrng1(e6p,p134,
     &    e3Cp)*s134mmtsq**(-1) + 4.D0*propL(2,jswap)*spstrng0(e2,e3Cp)
     &    *spstrng0(e6p,e1C)*spstrng1(e4p,p234,e5Cp)*s234mmtsq**(-1) - 
     &    4.D0*propR(2,jswap)*spstrng0(e3m,e1C)*spstrng0(e2,e6Cm)*
     &    spstrng1(e5m,p134,e4Cm)*s134mmtsq**(-1) + 4.D0*propR(2,jswap)
     &    *spstrng0(e5m,e1C)*spstrng0(e2,e4Cm)*spstrng1(e3m,p234,e6Cm)*
     &    s234mmtsq**(-1) )
      amp(h1,h2,1,jswap) = amp(h1,h2,1,jswap) + s34**(-1)*mt*s56**(-1)
     &  * ( 4.D0*propL(2,jswap)*spstrng0(e3m,e5Cp)*spstrng0(e2,e4Cm)*
     &    spstrng0(e6p,e1C)*s234mmtsq**(-1) + 4.D0*propL(2,jswap)*
     &    spstrng0(e3m,e1C)*spstrng0(e2,e5Cp)*spstrng0(e6p,e4Cm)*
     &    s134mmtsq**(-1) + 4.D0*propR(2,jswap)*spstrng0(e4p,e6Cm)*
     &    spstrng0(e5m,e1C)*spstrng0(e2,e3Cp)*s234mmtsq**(-1) + 4.D0*
     &    propR(2,jswap)*spstrng0(e4p,e1C)*spstrng0(e5m,e3Cp)*spstrng0(
     &    e2,e6Cm)*s134mmtsq**(-1) )
      amp(h1,h2,1,jswap) = amp(h1,h2,1,jswap) + s56**(-1)*s12**(-1)
     &  * ( 4.D0*propL(j34,jswap)*spstrng0(e3m,e5Cp)*spstrng0(e4p,e1C)*
     &    spstrng1(e2,p124,e6Cm)*s124**(-1) - 4.D0*propL(j34,jswap)*
     &    spstrng0(e3m,e5Cp)*spstrng0(e2,e4Cm)*spstrng1(e6p,p124,e1C)*
     &    s124**(-1) + 4.D0*propL(j34,jswap)*spstrng0(e3m,e1C)*
     &    spstrng0(e6p,e4Cm)*spstrng1(e2,p123,e5Cp)*s123**(-1) - 4.D0*
     &    propL(j34,jswap)*spstrng0(e2,e3Cp)*spstrng0(e6p,e4Cm)*
     &    spstrng1(e5m,p123,e1C)*s123**(-1) )

      amp(h1,h2,2,jswap)= + s34**(-1)*s56**(-1) * (  - 4.D0*propL(2,
     &    jswap)*spstrng0(e3p,e1C)*spstrng0(e2,e5Cp)*spstrng1(e6p,p134,
     &    e4Cp)*s134mmtsq**(-1) + 4.D0*propL(2,jswap)*spstrng0(e2,e4Cp)
     &    *spstrng0(e6p,e1C)*spstrng1(e3p,p234,e5Cp)*s234mmtsq**(-1) - 
     &    4.D0*propR(2,jswap)*spstrng0(e4m,e1C)*spstrng0(e2,e6Cm)*
     &    spstrng1(e5m,p134,e3Cm)*s134mmtsq**(-1) + 4.D0*propR(2,jswap)
     &    *spstrng0(e5m,e1C)*spstrng0(e2,e3Cm)*spstrng1(e4m,p234,e6Cm)*
     &    s234mmtsq**(-1) )
      amp(h1,h2,2,jswap) = amp(h1,h2,2,jswap) + s34**(-1)*mt*s56**(-1)
     &  * ( 4.D0*propL(2,jswap)*spstrng0(e4m,e5Cp)*spstrng0(e2,e3Cm)*
     &    spstrng0(e6p,e1C)*s234mmtsq**(-1) + 4.D0*propL(2,jswap)*
     &    spstrng0(e4m,e1C)*spstrng0(e2,e5Cp)*spstrng0(e6p,e3Cm)*
     &    s134mmtsq**(-1) + 4.D0*propR(2,jswap)*spstrng0(e3p,e6Cm)*
     &    spstrng0(e5m,e1C)*spstrng0(e2,e4Cp)*s234mmtsq**(-1) + 4.D0*
     &    propR(2,jswap)*spstrng0(e3p,e1C)*spstrng0(e5m,e4Cp)*spstrng0(
     &    e2,e6Cm)*s134mmtsq**(-1) )
      amp(h1,h2,2,jswap) = amp(h1,h2,2,jswap) + s56**(-1)*s12**(-1)
     &  * ( 4.D0*propR(j34,jswap)*spstrng0(e3p,e6Cm)*spstrng0(e4m,e1C)*
     &    spstrng1(e2,p124,e5Cp)*s124**(-1) - 4.D0*propR(j34,jswap)*
     &    spstrng0(e3p,e6Cm)*spstrng0(e2,e4Cp)*spstrng1(e5m,p124,e1C)*
     &    s124**(-1) + 4.D0*propR(j34,jswap)*spstrng0(e3p,e1C)*
     &    spstrng0(e5m,e4Cp)*spstrng1(e2,p123,e6Cm)*s123**(-1) - 4.D0*
     &    propR(j34,jswap)*spstrng0(e5m,e4Cp)*spstrng0(e2,e3Cm)*
     &    spstrng1(e6p,p123,e1C)*s123**(-1) )

      enddo
      enddo
      enddo
      ampsq=0d0
      do h1=1,2
      do h2=1,2
      do h34=1,2
      do h56=1,2
      ampsq=ampsq
     & +real(amp(h1,h2,h34,h56)*conjg(amp(h1,h2,h34,h56)))
      enddo
      enddo
      enddo
      enddo
      return
      end
