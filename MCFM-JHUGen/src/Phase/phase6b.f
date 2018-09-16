      subroutine phase6b(r,p1,p2,p3,p4,p5,p6,p7,p8,wt,*)
      implicit none
      include 'types.f'
c----phase space for signal
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'mxdim.f'
      include 'breit.f'
c********* generate phase space for 2-->6 process
c********* r(mxdim),p1(4),p2(4) are inputs 
c--------- incoming p1 and p2 reversed in sign from physical values 
c---- i.e. phase space for -p1-p2 --> p3+p4+p5+p6+p7+p8
c---- with all 2 pi's (ie 1/(2*pi)^11)

      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4)
      real(dp):: p12(4),p56(4),p34(4),p5678(4),p78(4)
      real(dp):: wt,wt12,wt5678,wt567,wt34,wt56,wt0
      integer:: j

      parameter(wt0=1._dp/twopi**4)

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      
c--- Hack to generate points with p8 soft
c      r(7)=r(7)*1.e-8_dp
c--- Hack to generate points with p1.p8 small
c      r(2)=r(2)*1.e-8_dp
c--- Hack to generate points with p2.p8 small
c      r(2)=1._dp-r(2)*1.e-8_dp
      
c--- Note that r(9) and r(10) are used in gen6 to generate x1 and x2
c      call phi1_2(r(1),r(2),r(3),r(4),p12,p5678,p34,wt12,*99)  
c      call phi3m0(r(5),r(6),p34,p3,p4,wt34,*99)
c      call phi1_2m_nobw(zip,r(7),r(8),r(11),smin,p5678,p8,p567,
c     &   wt5678,*99)
c      call phi1_2m_bw(zip,r(12),r(13),r(14),smin,p567,p7,p56,
c     &   wmass,wwidth,wt567,*99)
c      call phi3m0(r(15),r(16),p56,p5,p6,wt56,*99)
c     
c      wt=wt0*wt12*wt34*wt5678*wt567*wt56

c--- Alternative generation
      call phi1_2(r(1),r(2),r(3),r(4),p12,p5678,p34,wt12,*99)
      call phi3m0(r(5),r(6),p34,p3,p4,wt34,*99)
      n2=0
      call phi1_2(r(7),r(8),r(11),r(12),p5678,p78,p56,wt5678,*99)  
      n2=1
      call phi3m0(r(13),r(14),p56,p5,p6,wt567,*99)
      call phi3m0(r(15),r(16),p78,p7,p8,wt56,*99)
      
      wt=wt0*wt12*wt34*wt5678*wt567*wt56

      return

 99   continue
      n2=1
      wt=0._dp
      return 1
      end

