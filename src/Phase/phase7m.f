      subroutine phase7m(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,wt,*)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'breit.f'
c******* generate phase space for 2-->4 process
c******* r(mxdim),p1(4),p2(4) are inputs reversed in sign from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6+p7+p8+p9+p10
c---- with all 2 pi's (ie 1/(2*pi)^20)
      integer j
      double precision r(mxdim)
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4),
     . p9(4),p12(4),
     . p345(4),p678(4),p45(4),p78(4),p345678(4),
     . smin,wt,wt0,wt129,wt345678,wt345,wt678,wt45,wt78
      parameter(wt0=1d0/twopi**5)
      wt=0d0
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      smin=100d0

      n2=1
      n3=1

      mass3=hmass
      width3=hwidth
      call phi1_2m(0d0,r(1),r(2),r(3),smin,p12,p9,p345678,wt129,*99)
      
      n2=0
      n3=0

      mass2=mtau
      width2=tauwidth
      mass3=mtau
      width3=tauwidth

      call phi1_2(r(4),r(5),r(6),r(7),p345678,p345,p678,wt345678,*99)
      
      n2=1
      n3=1

      mass3=wmass
      width3=wwidth
      call phi1_2m(0d0,r(8),r(9),r(10),smin,p345,p3,p45,wt345,*99)
      call phi1_2m(0d0,r(11),r(12),r(13),smin,p678,p6,p78,wt678,*99)
      
      if ((p3(4).le.0d0).or.(p6(4).le.0d0)) goto 99
      call phi3m0(r(14),r(15),p45,p4,p5,wt45,*99)
      if ((p4(4).le.0d0).or.(p5(4).le.0d0)) goto 99
      call phi3m0(r(16),r(17),p78,p7,p8,wt78,*99)
      if ((p7(4).le.0d0).or.(p8(4).le.0d0)) goto 99

      wt=wt0*wt129*wt345678*wt345*wt678*wt45*wt78

      return
 99   wt=0d0
c      write(*,*) 'wt129',wt129
c      write(*,*) 'wt345678',wt345678
c      write(*,*) 'wt345',wt345
c      write(*,*) 'wt678',wt678
c      write(*,*) 'wt45',wt45
c      write(*,*) 'wt78',wt78
      return 1
      end

