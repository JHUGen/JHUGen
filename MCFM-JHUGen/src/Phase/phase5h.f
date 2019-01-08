      subroutine phase5h(r,p1,p2,p3,p4,p5,p6,p7,wt)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'breit.f'
c********* generate phase space for 2-->5 process
c********* r(mxdim),p1(4),p2(4) are inputs 
c--------- incoming p1 and p2 reversed in sign from physical values 
c---- i.e. phase space for -p1-p2 --> p3+p4+p5+p6+p7
c---- with all 2 pi's (ie 1/(2*pi)^11)

      double precision r(mxdim),wt,wt34567,wt34,wt567,wt67,wt0,tmp
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      double precision p567(4),p12(4),p67(4),p34(4)
      integer j
      parameter(wt0=1d0/twopi**3)

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo

      n3=1
      call phi1_2(r(1),r(2),r(3),r(4),p12,p567,p34,wt34567,*99)
      n3=0
      call phi1_2m(0d0,r(5),r(6),r(7),0d0,p567,p5,p67,wt567,*99)
      call phi3m0(r(8),r(11),p34,p3,p4,wt34,*99)
      call phi3m0(r(12),r(13),p67,p6,p7,wt67,*99)

      wt=wt0*wt34567*wt567*wt67*wt34
  
      return
 99   continue
      wt=0d0
      return
      end

