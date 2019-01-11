      subroutine phase7(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,wt,*)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'process.f'
      include 'breit.f'
c******* generate phase space for 2-->4 process
c******* r(mxdim),p1(4),p2(4) are inputs reversed in sign from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6+p7+p8+p9+p10
c---- with all 2 pi's (ie 1/(2*pi)^20)
      integer nu,iflip,j
      double precision r(mxdim)
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4),
     . p9(4),p12(4),pa(4),pb(4),
     . p345(4),p678(4),p34(4),p78(4),
     . smin,wt,wt0,wt12,wt345,wt678,wt34,wt78,wt9
      parameter(wt0=1d0/twopi**5)
      data iflip/0/
      save iflip
      wt=0d0
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      smin=100d0

      n2=0
      n3=0
      if  ((case .eq. 'qq_ttg') 
     . .or.(case .eq. 'tt_bbl') 
     . .or.(case .eq. 'tt_bbh') 
     . .or.(case .eq. 'tt_bbu')) then 
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth
      elseif (case .eq. 'hlljet') then
        mass2=mtau
        width2=tauwidth
        mass3=mtau
        width3=tauwidth
      else
        write(*,*) 'Bad case in phase7.f'
        stop
      endif
      
c     massive particle p12 decaying into pa pb
c     with invariant mass 
c     of particle two s2 and particle three s3 integrated over.
c     vectors returned p2 and p3 are in the same frame as p1 is supplied
      call phi1_2(r(1),r(2),r(3),r(4),p12,pa,pb,wt12,*99)

      if (case .eq. 'qq_ttg') then
c--- branchings are all of the form A->B+C with B^2=0 and C^2 a B.-W.
      n2=0
      n3=1
      endif
      
      if (iflip .eq. 0) then
C----emission of parton p9
        iflip=1
        call phi1_2m(0d0,r(5),r(6),r(7),smin,pa,p9,p345,wt9,*99)
        do nu=1,4
        p678(nu)=pb(nu)
        enddo
      elseif (iflip .eq. 1) then 
C----emission of parton p9
        iflip=0
        call phi1_2m(0d0,r(5),r(6),r(7),smin,pa,p9,p678,wt9,*99)
        do nu=1,4
        p345(nu)=pb(nu)
        enddo
      endif 

      mass3=wmass
      width3=wwidth
c--decay of p345 into p5 and p34
      call phi1_2m(mb,r(8),r(9),r(10),smin,p345,p5,p34,wt345,*99)
c--decay of p678 into p6 and p78
      call phi1_2m(mb,r(11),r(12),r(13),smin,p678,p6,p78,wt678,*99)

      if ((p5(4).le.0d0).or.(p6(4).le.0d0)) goto 99
      call phi3m0(r(14),r(15),p34,p3,p4,wt34,*99)
      if ((p3(4).le.0d0).or.(p4(4).le.0d0)) goto 99
      call phi3m0(r(16),r(17),p78,p7,p8,wt78,*99)
      if ((p7(4).le.0d0).or.(p8(4).le.0d0)) goto 99

      wt=wt0*wt12*wt9*wt345*wt678*wt34*wt78
      
      return
 99   wt=0d0
      return 1
      end

