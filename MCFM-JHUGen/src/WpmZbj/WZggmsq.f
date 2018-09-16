      function WZggmsq(p1,p2,p3,p4,p5,p6,p7,p8)
      implicit none
      include 'types.f'
      real(dp):: WZggmsq

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'nwz.f'
      complex(dp)::prop34,prop56,prop3456
      complex(dp)::AB3456(2,2),BA3456(2,2),AB5634(2,2),BA5634(2,2),
     & NAB(2,2),NBA(2,2),ab2w(2,2),ba2w(2,2),absra(2,2),basra(2,2),
     & absrl(2,2),basrl(2,2),
     & AB1(2,2,2),AB2(2,2,2),NAB1(2,2,2),NAB2w(2,2,2),
     & NABsra(2,2,2),NABsrl(2,2,2),
     & BA1(2,2,2),BA2(2,2,2),NBA1(2,2,2),NBA2w(2,2,2),
     & NBAsra(2,2,2),NBAsrl(2,2,2),
     & AmpAB(2,2,2),AmpBA(2,2,2),AmpQED(2,2,2)
      real(dp):: s3456,s4,fac,q3,q4,l3,l4
      integer:: jd,ju,p1,p2,p3,p4,p5,p6,p7,p8,h56,h7,h8
C     statement functions
      s4(p1,p2,p3,p4)=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &               +s(p2,p3)+s(p2,p4)+s(p3,p4)
C     end statement functions

      fac=4._dp*V*xn*gwsq**2*esq**2*gsq**2


      prop56=s(p5,p6)/cplx2(s(p5,p6)-zmass**2,zmass*zwidth)
      prop34=s(p3,p4)/cplx2(s(p3,p4)-wmass**2,wmass*wwidth)
      s3456=s4(p3,p4,p5,p6)
      prop3456=s3456/cplx2(s3456-wmass**2,wmass*wwidth)
C------Setup Amplitudes for LH 56 line;
      call TWZggAB(p1,p2,p3,p4,p5,p6,p7,p8,AB3456)
      call TWZggAB(p1,p2,p3,p4,p5,p6,p8,p7,BA3456)
      call TWZggAB(p1,p2,p5,p6,p3,p4,p7,p8,AB5634)
      call TWZggAB(p1,p2,p5,p6,p3,p4,p8,p7,BA5634)

      call TWZggNAB(p1,p2,p3,p4,p5,p6,p7,p8,NAB)
      call TWZggNAB(p1,p2,p3,p4,p5,p6,p8,p7,NBA)

c      call TWZgg2w(p1,p2,p3,p4,p5,p6,p7,p8,AB2w)
c      call TWZgg2w(p1,p2,p3,p4,p5,p6,p8,p7,BA2w)
      if (nwz == 1) then
      call TWZggSRL(p1,p2,p5,p6,p3,p4,p7,p8,AB2w)
      call TWZggSRL(p1,p2,p5,p6,p3,p4,p8,p7,BA2W)
      elseif (nwz == -1) then
      call TWZggSRA(p1,p2,p5,p6,p3,p4,p7,p8,AB2w)
      call TWZggSRA(p1,p2,p5,p6,p3,p4,p8,p7,BA2W)
      endif
      call TWZggSRA(p1,p2,p3,p4,p5,p6,p7,p8,ABsra)
      call TWZggSRA(p1,p2,p3,p4,p5,p6,p8,p7,BAsra)
      call TWZggSRL(p1,p2,p3,p4,p5,p6,p7,p8,ABsrl)
      call TWZggSRL(p1,p2,p3,p4,p5,p6,p8,p7,BAsrl)


C     First index is helicity of Z (56) decay line
      do h7=1,2
      do h8=1,2
      AB1(1,h7,h8)=AB3456(h7,h8)
      AB2(1,h7,h8)=AB5634(h7,h8)
      NAB1(1,h7,h8)=NAB(h7,h8)
      NAB2w(1,h7,h8)=AB2w(h7,h8)
      NABsra(1,h7,h8)=ABsra(h7,h8)
      NABsrl(1,h7,h8)=ABsrl(h7,h8)

      BA1(1,h7,h8)=BA3456(h8,h7)
      BA2(1,h7,h8)=BA5634(h8,h7)
      NBA1(1,h7,h8)=NBA(h8,h7)
      NBA2w(1,h7,h8)=BA2w(h8,h7)
      NBAsra(1,h7,h8)=BAsra(h8,h7)
      NBAsrl(1,h7,h8)=BAsrl(h8,h7)
      enddo
      enddo

C------Setup Amplitudes for RH 56 line;
      call TWZggAB(p1,p2,p3,p4,p6,p5,p7,p8,AB3456)
      call TWZggAB(p1,p2,p6,p5,p3,p4,p7,p8,AB5634)
      call TWZggNAB(p1,p2,p3,p4,p6,p5,p7,p8,NAB)
      call TWZggAB(p1,p2,p3,p4,p6,p5,p8,p7,BA3456)
      call TWZggAB(p1,p2,p6,p5,p3,p4,p8,p7,BA5634)
      call TWZggNAB(p1,p2,p3,p4,p6,p5,p8,p7,NBA)

      call TWZggSRA(p1,p2,p3,p4,p6,p5,p7,p8,ABsra)
      call TWZggSRA(p1,p2,p3,p4,p6,p5,p8,p7,BAsra)
      call TWZggSRL(p1,p2,p3,p4,p6,p5,p7,p8,ABsrl)
      call TWZggSRL(p1,p2,p3,p4,p6,p5,p8,p7,BAsrl)


      do h7=1,2
      do h8=1,2
      AB1(2,h7,h8)=AB3456(h7,h8)
      AB2(2,h7,h8)=AB5634(h7,h8)
      NAB1(2,h7,h8)=NAB(h7,h8)
      NAB2w(2,h7,h8)=czip
      NBA2w(2,h7,h8)=czip
      NABsra(2,h7,h8)=ABsra(h7,h8)
      NABsrl(2,h7,h8)=ABsrl(h7,h8)

      BA1(2,h7,h8)=BA3456(h8,h7)
      BA2(2,h7,h8)=BA5634(h8,h7)
      NBA1(2,h7,h8)=NBA(h8,h7)
      NBAsra(2,h7,h8)=BAsra(h8,h7)
      NBAsrl(2,h7,h8)=BAsrl(h8,h7)
      enddo
      enddo

c      AB1=czip
c      AB2=czip
c      NAB1=czip
c      NAB2w=czip
c      NABsra=czip
c      NABsrl=czip

c      BA1=czip
c      BA2=czip
c      NBA1=czip
c      NBA2w=czip
c      NBAsra=czip
c      NBAsrl=czip


C     For u->dW^+ process Z is emitted after W for AB3456 and before W for AB5634
      if (nwz == 1) then
      ju=2
      jd=1
      q3=0._dp
      l3=ln
      q4=qe
      l4=le
      elseif (nwz == -1) then
      jd=2
      ju=1
      q3=qe
      l3=le
      q4=0._dp
      l4=ln
      endif

      do h7=1,2
      do h8=1,2
      ampAB(1,h7,h8)=
     &  ((Q(jd)*qe+L(jd)*le*prop56)*AB1(1,h7,h8)
     &  +(Q(ju)*qe+L(ju)*le*prop56)*AB2(1,h7,h8)
     &  +((Q(ju)-Q(jd))*qe+(L(ju)-L(jd))*le*prop56)
     &  *NAB1(1,h7,h8)*prop3456
     &  +0.5_dp/xw*NAB2w(1,h7,h8)*prop3456)*prop34
     &  +(q3*qe+l3*le*prop56)*NABsrl(1,h7,h8)*prop3456
     &  +(q4*qe+l4*le*prop56)*NABsra(1,h7,h8)*prop3456

      ampAB(2,h7,h8)=
     &  ((Q(jd)*qe+L(jd)*re*prop56)*AB1(2,h7,h8)
     &  +(Q(ju)*qe+L(ju)*re*prop56)*AB2(2,h7,h8)
     &  +((Q(ju)-Q(jd))*qe+(L(ju)-L(jd))*re*prop56)
     &  *NAB1(2,h7,h8)*prop3456
     &  +0.5_dp/xw*NAB2w(2,h7,h8)*prop3456)*prop34
     &  +(q3*qe+l3*re*prop56)*NABsrl(2,h7,h8)*prop3456
     &  +(q4*qe+l4*re*prop56)*NABsra(2,h7,h8)*prop3456

      ampBA(1,h7,h8)=
     &  ((Q(jd)*qe+L(jd)*le*prop56)*BA1(1,h7,h8)
     &  +(Q(ju)*qe+L(ju)*le*prop56)*BA2(1,h7,h8)
     &  +((Q(ju)-Q(jd))*qe+(L(ju)-L(jd))*le*prop56)
     &  *NBA1(1,h7,h8)*prop3456
     &  +0.5_dp/xw*NBA2w(1,h7,h8)*prop3456)*prop34
     &  +(q3*qe+l3*le*prop56)*NBAsrl(1,h7,h8)*prop3456
     &  +(q4*qe+l4*le*prop56)*NBAsra(1,h7,h8)*prop3456

      ampBA(2,h7,h8)=
     &  ((Q(jd)*qe+L(jd)*re*prop56)*BA1(2,h7,h8)
     &  +(Q(ju)*qe+L(ju)*re*prop56)*BA2(2,h7,h8)
     &  +((Q(ju)-Q(jd))*qe+(L(ju)-L(jd))*re*prop56)
     &  *NBA1(2,h7,h8)*prop3456
     &  +0.5_dp/xw*NBA2w(2,h7,h8)*prop3456)*prop34
     &  +(q3*qe+l3*re*prop56)*NBAsrl(2,h7,h8)*prop3456
     &  +(q4*qe+l4*re*prop56)*NBAsra(2,h7,h8)*prop3456
      enddo
      enddo

      WZggmsq=0._dp
      do h56=1,2
      do h7=1,2
      do h8=1,2
      AmpQED(h56,h7,h8)=AmpAB(h56,h7,h8)+AmpBA(h56,h7,h8)
      WZggmsq=WZggmsq
     & +fac*(+abs(AmpAB(h56,h7,h8))**2+abs(AmpBA(h56,h7,h8))**2
     & -1._dp/xn**2*abs(AmpQED(h56,h7,h8))**2)
      enddo
      enddo
      enddo
      return
      end
