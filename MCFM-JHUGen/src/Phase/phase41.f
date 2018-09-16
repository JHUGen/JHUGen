      subroutine phase41(r,p1,p2,p3,p4,p5,p6,wt,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'limits.f'
      include 'mxdim.f'
      include 'zerowidth.f'
      include 'breit.f'
c********* generate phase space for 2-->4 process
c********* r(mxdim),p1(4),p2(4) are inputs reversed in sign from physical values 
c---- phase space for -p1-p2 --> p4+p5+p6+p7
c---- with all 2 pi's (ie 1/(2*pi)^8)
      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      real(dp):: p12(4),p345(4),p34(4),s345min
      real(dp):: wt,wt3456,wt345,wt34,wt0
      logical:: oldzerowidth
      real(dp):: Mbbsq
      
      integer:: j
      parameter(wt0=1._dp/twopi**2)
      wt=0._dp
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      s345min=mb**2
c---- calculate momenta of top and bbbar

      n2=0
      n3=1
      mass3=mt
      width3=twidth
      
      oldzerowidth=zerowidth
      zerowidth=.true.

      call phi1_2m(mb,r(1),r(2),r(3),s345min,p12,p6,p345,wt3456,*99)

      zerowidth=oldzerowidth

      n3=1
      mass3=wmass
      width3=wwidth
      call phi1_2m(mb,r(4),r(5),r(6),zip,p345,p5,p34,wt345,*99)

      Mbbsq=2*(mb**2+p5(4)*p6(4)-p5(1)*p6(1)-p5(2)*p6(2)-p5(3)*p6(3))
      if ((Mbbsq > bbsqmax) .or. (Mbbsq < bbsqmin)) return 1
      
      call phi3m0(r(7),r(8),p34,p3,p4,wt34,*99)

      wt=wt0*wt3456*wt345*wt34

      return
 99   wt=0._dp
      zerowidth=oldzerowidth
      return 1
      end

