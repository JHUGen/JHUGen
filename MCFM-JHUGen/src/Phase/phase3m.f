      subroutine phase3m(r,p1,p2,p3,p4,p5,p6,p7,m3,m4,m5,wt)
      implicit none
      include 'types.f'
c----generate phase space for 2-->3 process with masses m3,m4,m5
c----r(mxdim),p1(4),p2(4) are inputs 
c----incoming p1 and p2 reversed in sign from physical values 
c----i.e. phase space for -p1-p2 --> p3+p4+p5
c----with all 2 pi's (ie 1/(2*pi)^5)
c----(p4,p5) are dummies
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'mxdim.f'

      integer:: j
      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      real(dp):: p12(4),p34(4),smin
      real(dp):: wt,wt125,wt34,wt0,m3,m4,m5
      parameter(wt0=1._dp/twopi)

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      p6(j)=0._dp
      p7(j)=0._dp
      enddo
      smin=(m3+m4)**2

c---generate p5 and p34, 
c---smin is the minimum inv mass of 34 system
c---m5 is the mass of p5
      call phi1_2m(m5,r(1),r(2),r(3),smin,p12,p5,p34,wt125,*99)

c---decay 34-system
      call phi3m(r(4),r(5),p34,p3,p4,m3,m4,wt34,*99)

      wt=wt0*wt125*wt34
      return
 99   continue
      wt=0._dp
      return
      end

