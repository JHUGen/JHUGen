      subroutine ttqqHampsq(q,q1,q2,q3,q4,ampsqf)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      real(dp):: q(mxpart,4),s134mmtsq,s234mmtsq,s134,s234,s34
      real(dp):: p3Dp4,ampsqf
      complex(dp):: spstrng0,spstrng1,p1(4),p2(4),p3(4),p4(4),
     & e1C(4),e1Cp(4),e1Cm(4),e2(4),e2p(4),e2m(4),
     & e3(4),e3C(4),e4(4),e4C(4),
     & e3m(4),e4m(4),e3p(4),e4p(4),e3Cm(4),e4Cm(4),e3Cp(4),e4Cp(4),
     & p134(4),p234(4),amp(2,2,2,2)
      integer:: q1,q2,q3,q4,h1,h2,h3,h4
      p3Dp4=
     & q(q3,4)*q(q4,4)-q(q3,1)*q(q4,1)-q(q3,2)*q(q4,2)-q(q3,3)*q(q4,3)
      p1(:)=cmplx(q(q1,:),kind=dp)
      p2(:)=cmplx(q(q2,:),kind=dp)
      p3(:)=cmplx(q(q3,:),kind=dp)
      p4(:)=cmplx(q(q4,:),kind=dp)
      p134(:)=cmplx(q(q1,:)+q(q3,:)+q(q4,:),kind=dp)
      p234(:)=cmplx(q(q2,:)+q(q3,:)+q(q4,:),kind=dp)
      s134=real(p134(4)**2-p134(1)**2-p134(2)**2-p134(3)**2)
      s234=real(p234(4)**2-p234(1)**2-p234(2)**2-p234(3)**2)
      s34=2d0*p3Dp4
      s134mmtsq=s134-mt**2
      s234mmtsq=s234-mt**2
      h1=-1
      call VklSt(p1,mt,p3,h1,e1Cm)
      call UbKlst(p2,mt,p3,h1,e2m)
      call ubarspinor0(p3,h1,e3m)
      call uspinor0(p3,-h1,e3Cm)
      call ubarspinor0(p4,h1,e4m)
      call uspinor0(p4,-h1,e4Cm)
      h1=1
      call VKlSt(p1,mt,p3,h1,e1Cp)
      call UbKlSt(p2,mt,p3,h1,e2p)
      call ubarspinor0(p3,h1,e3p)
      call uspinor0(p3,-h1,e3Cp)
      call ubarspinor0(p4,h1,e4p)
      call uspinor0(p4,-h1,e4Cp)
      ampsqf=0d0
      amp(:,:,:,:)=czip
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
      amp(h1,h2,h3,h4)= + two*s34**(-1)*s134mmtsq**(-1) * (  - 
     &    spstrng0(e3,e1C)*spstrng1(e2,p134,e4C) - spstrng0(e4,e1C)*
     &    spstrng1(e2,p134,e3C) )
      amp(h1,h2,h3,h4) = amp(h1,h2,h3,h4) + two*s34**(-1)*
     & s234mmtsq**(-1) * ( spstrng0(e2,e3C)*spstrng1(e4,p234,e1C) + 
     &    spstrng0(e2,e4C)*spstrng1(e3,p234,e1C) )
      amp(h1,h2,h3,h4) = amp(h1,h2,h3,h4) + mt*two*s34**(-1)*
     & s134mmtsq**(-1) * ( spstrng0(e2,e3C)*spstrng0(e4,e1C) + 
     &    spstrng0(e2,e4C)*spstrng0(e3,e1C) )
      amp(h1,h2,h3,h4) = amp(h1,h2,h3,h4) + mt*two*s34**(-1)*
     & s234mmtsq**(-1) * ( spstrng0(e2,e3C)*spstrng0(e4,e1C) + 
     &    spstrng0(e2,e4C)*spstrng0(e3,e1C) )

      ampsqf=ampsqf
     & +real(amp(h1,h2,h3,h4)*conjg(amp(h1,h2,h3,h4)))
      enddo
      enddo
      enddo
      return
      end
