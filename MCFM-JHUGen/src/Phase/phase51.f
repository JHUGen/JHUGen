      subroutine phase51(r,p1,p2,p3,p4,p5,p6,p7,wt)
      implicit none
      include 'types.f'
c----phase space for qg--> t(\nu(3) b(5) e^+(4))  bbar(6) q'(7) 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'kprocess.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'breit.f'
c********* generate phase space for 2-->5 process
c********* r(mxdim),p1(4),p2(4) are inputs 
c--------- incoming p1 and p2 reversed in sign from physical values 
c---- i.e. phase space for -p1-p2 --> p3+p4+p5+p6+p7
c---- with all 2 pi's (ie 1/(2*pi)^11)

      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      real(dp):: p127(4),p12(4),p345(4),p34(4),smin
      real(dp):: wt,wt127,wt3456,wt345,wt34,wt0,bmass
      integer:: j
      parameter(wt0=1._dp/twopi**3)

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo

      if (kcase==kt_bbar) then
      bmass=0._dp
      smin=0._dp
      else
      bmass=mb
      smin=4._dp*bmass**2
      endif
      call 
     & phi1_2m_nobw(0._dp,r(13),r(12),r(11),smin,p12,p7,p127,wt127,*99)

      n2=0
      n3=1
      mass3=mt
      width3=twidth
      call phi1_2m(bmass,r(1),r(2),r(3),smin,p127,p6,p345,wt3456,*99)
      n3=1
      mass3=wmass
      width3=wwidth

      call phi1_2m(bmass,r(4),r(5),r(6),smin,p345,p5,p34,wt345,*99)

      call phi3m0(r(7),r(8),p34,p3,p4,wt34,*99)
      wt=wt0*wt127*wt3456*wt345*wt34

      if (debug) write(6,*) 'wt123',wt127
      if (debug) write(6,*) 'wt3456',wt3456
      if (debug) write(6,*) 'wt345',wt345
      if (debug) write(6,*) 'wt34',wt34

      return
 99   continue
      wt=0._dp
      return
      end

