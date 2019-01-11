      subroutine phase6a(r,p1,p2,p3,p4,p5,p6,p7,p8,wt,*)
c----phase space for signal
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'process.f'
      include 'zerowidth.f'
c********* generate phase space for 2-->6 process
c********* r(mxdim),p1(4),p2(4) are inputs 
c--------- incoming p1 and p2 reversed in sign from physical values 
c---- i.e. phase space for -p1-p2 --> p3+p4+p5+p6+p7+p8
c---- with all 2 pi's (ie 1/(2*pi)^11)

      double precision r(mxdim)
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4)
      double precision p12(4),p56(4),p34(4),p567(4),p34567(4),smin,s3min
      double precision wt,wt12,wt34567,wt567,wt34,wt56,wt0
      integer j

      parameter(wt0=1d0/twopi**4)

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      smin=mb
      s3min=mt**2
      if (zerowidth) then
        s3min=(mt+wmass)**2
      endif
c--- the parameters above are appropriate for W+t processes only
      if (case .eq. 'W_cwdk') then
        smin=zip
      s3min=mc**2
      endif


c--- Hack to generate points with p8 soft
c      r(1)=1d0-r(1)*1d-8
c--- Hack to generate points with p1.p8 small
c       r(2)=r(2)*1d-8
c--- Hack to generate points with p2.p8 small
c      r(2)=1d0-r(2)*1d-8
      
c--- Note that r(9) and r(10) are used in gen6 to generate x1 and x2
      call phi1_2m_nobw(smin,r(1),r(2),r(3),s3min,
     . p12,p8,p34567,wt12,*99)  
      call phi1_2(r(4),r(5),r(6),r(7),p34567,p567,p34,wt34567,*99)
      call phi3m0(r(8),r(11),p34,p3,p4,wt34,*99)
      call phi1_2m_bw(zip,r(12),r(13),r(14),zip,p567,p7,p56,
     .   wmass,wwidth,wt567,*99)
      call phi3m0(r(15),r(16),p56,p5,p6,wt56,*99)
      
      wt=wt0*wt12*wt34567*wt34*wt567*wt56

      return

 99   continue
      wt=0d0
      return 1
      end

