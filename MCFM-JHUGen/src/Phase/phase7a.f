      subroutine phase7a(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,wt,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'zerowidth.f'
      include 'kprocess.f'
      include 'breit.f'
c******* generate phase space for 2-->5 process
c******* r(mxdim),p1(4),p2(4) are inputs reversed in sign from physical values 
c---- phase space for -p1-p2 --> p5+p6+p3+p4+p7+p8+p9
c---- with all 2 pi's (ie 1/(2*pi)^17)
      logical:: oldzerowidth
      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p5(4),p6(4),p3(4),p4(4),p7(4),
     & p8(4),p9(4)
      real(dp):: p12(4),p789(4),p34(4),p78(4),p56(4),p3456(4)
      real(dp):: wt,wt0,wt12,wt789,wt34,wt78,wt3456,wt56
      integer:: j
      parameter(wt0=1._dp/twopi**5)

c--- written for real contribution to qq->H(->WW)+qq only
      if  ((kcase==kqq_HWW) .or. (kcase==kqq_HZZ)
     ..or. (kcase==kHWW2jt) .or. (kcase==kHZZ2jt)
     ..or. (kcase==kHWW3jt) .or. (kcase==kHZZ3jt)
     ..or. (kcase==kWpWp3j)) then
        continue
      else 
        write(6,*) 'Phase space routine not correct - needs updating.'
        write(6,*) 'kcase',kcase
        stop
      endif

      oldzerowidth=zerowidth

      wt=0._dp
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo

c--- In the case of HWW/HZZ+2jets, we should generate s3456 according
c--- to a Breit-Wigner at mH

      if (kcase==kWpWp3j) then
        n2=0
        mass2=0
        width2=0
        n3=0
      else
        n2=1
        mass2=hmass
        width2=hwidth
        n3=0
      endif

      call phi1_2(r(1),r(2),r(3),r(4),p12,p3456,p789,wt12,*99)

      if     ((kcase==kHWW2jt) .or. (kcase==kqq_HWW)
     &   .or. (kcase==kHWW3jt)) then
        n2=1
        mass2=wmass
        width2=wwidth
        n3=1
        mass3=wmass
        width3=wwidth
      elseif ((kcase==kHZZ2jt) .or. (kcase==kqq_HZZ)
     &   .or. (kcase==kHZZ3jt)) then
        n2=1
        mass2=zmass
        width2=zwidth
        n3=1
        mass3=zmass
        width3=zwidth
      elseif (kcase==kWpWp3j) then
        n2=1
        mass2=wmass
        width2=wwidth
        n3=1
        mass3=wmass
        width3=wwidth
      endif

      call phi1_2(r(5),r(6),r(7),r(8),p3456,p56,p34,wt3456,*99)
      call phi3m0(r(13),r(14),p34,p3,p4,wt34,*99)
      call phi3m0(r(11),r(12),p56,p5,p6,wt56,*99)
      call phi1_2m_nobw(zip,r(15),r(16),r(17),zip,p789,p9,p78,wt789,*99)
      call phi3m0(r(18),r(19),p78,p7,p8,wt78,*99)
      wt=wt0*wt12*wt3456*wt34*wt56*wt789*wt78
      return
      
 99   wt=0._dp
      zerowidth=oldzerowidth
      return 1
      end

