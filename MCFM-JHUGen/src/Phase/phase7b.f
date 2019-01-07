      subroutine phase7b(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,wt,*)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'zerowidth.f'
      include 'process.f'
      include 'breit.f'
c******* generate phase space for 2-->5 process
c******* r(mxdim),p1(4),p2(4) are inputs reversed in sign from physical values 
c---- phase space for -p1-p2 --> p5+p6+p3+p4+p7+p8+p9
c---- with all 2 pi's (ie 1/(2*pi)^17)
      logical oldzerowidth
      double precision r(mxdim)
      double precision p1(4),p2(4),p5(4),p6(4),p3(4),p4(4),p7(4),
     . p8(4),p9(4)
      double precision p12(4),p349(4),p34(4),
     . p78(4),p56(4),p5678(4)
      double precision wt,wt0,wt12,wt349,wt34,wt78,wt5678,wt56
      integer j
      parameter(wt0=1d0/twopi**5)

c--- written for real contribution to qq->WH(->WW)+g (and ZH) only
      if    ((case .ne. 'WH__WW')
     & .and. (case .ne. 'ZH__WW')
     & .and. (case .ne. 'WH__ZZ')
     & .and. (case .ne. 'ZH__ZZ')) then
        write(6,*) 'Phase space routine not correct - needs updating.'
        stop
      endif

      oldzerowidth=zerowidth

      wt=0d0
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo

c--- In these cases, we should generate s5678 according to
c--- a Breit-Wigner at mH
      n2=0
      mass3=hmass
      width3=hwidth
      n3=1
      call phi1_2(r(1),r(2),r(3),r(4),p12,p349,p5678,wt12,*99)
      if  ((case .eq. 'WH__WW') .or. (case .eq. 'ZH__WW')) then
      n2=1
      mass2=wmass
      width2=wwidth
      n3=1
      mass3=wmass
      width3=wwidth
      elseif ((case .eq. 'WH__ZZ') .or. (case .eq. 'ZH__ZZ')) then
      n2=1
      mass2=zmass
      width2=zwidth
      n3=1
      mass3=zmass
      width3=zwidth
      endif
      call phi1_2(r(5),r(6),r(7),r(8),p5678,p56,p78,wt5678,*99)
      call phi3m0(r(11),r(12),p56,p5,p6,wt56,*99)
      call phi3m0(r(13),r(14),p78,p7,p8,wt78,*99)
      n3=1
      if ((case .eq. 'WH__WW') .or. (case .eq. 'WH__ZZ')) then
        mass3=wmass
        width3=wwidth
      else
        mass3=zmass
        width3=zwidth
      endif
      call phi1_2m(zip,r(15),r(16),r(17),zip,p349,p9,p34,wt349,*99)
      call phi3m0(r(18),r(19),p34,p3,p4,wt34,*99)
      wt=wt0*wt12*wt5678*wt56*wt78*wt349*wt34
      return
      
 99   wt=0d0
      zerowidth=oldzerowidth
      return 1
      end

