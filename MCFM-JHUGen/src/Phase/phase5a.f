      subroutine phase5a(r,p1,p2,p3,p4,p5,p6,p7,wt)
      implicit none
      include 'types.f'
c----phase space for signal
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'mxdim.f'
c********* generate phase space for 2-->5 process
c********* r(mxdim),p1(4),p2(4) are inputs 
c--------- incoming p1 and p2 reversed in sign from physical values 
c---- i.e. phase space for -p1-p2 --> p3+p4+p5+p6+p7
c---- with all 2 pi's (ie 1/(2*pi)^11)

      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      real(dp):: p12(4),p56(4),p34(4),p567(4),smin
      real(dp):: wt,wt34567,wt567,wt34,wt56,wt0
      integer:: j

      parameter(wt0=1._dp/twopi**3)

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      smin=0._dp

c--- Note that r(9) and r(10) are used in gen5 to generate x1 and x2
      call phi1_2(r(1),r(2),r(3),r(4),p12,p567,p34,wt34567,*99)
      call phi3m0(r(5),r(6),p34,p3,p4,wt34,*99)
      call phi1_2m_bw(zip,r(7),r(8),r(11),smin,p567,p7,p56,
     &   wmass,wwidth,wt567,*99)
      call phi3m0(r(12),r(13),p56,p5,p6,wt56,*99)
      
      wt=wt0*wt34567*wt34*wt567*wt56

      return
 99   continue
      wt=0._dp
      return
      end

