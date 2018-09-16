      subroutine phase7m_alt(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,m3,m4,m5,wt)
      implicit none
      include 'types.f'
c----generate phase space for 2-->7 process with masses m3,m4,m5
c----r(mxdim),p1(4),p2(4) are inputs 
c----incoming p1 and p2 reversed in sign from physical values 
c----i.e. phase space for -p1-p2 --> (p3+p4+p5)+(p6+p7+p8)+p9
c----with all 2 pi's (ie 1/(2*pi)^5)
c----(p4,p5) are dummies
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'zerowidth.f'
      include 'breit.f'
      integer:: j
      logical:: oldzerowidth
      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4),
     & p9(4),p12(4),p345678(4),p345(4),p678(4),p34(4),p78(4),smin
      real(dp):: wt0,wt12,wt345678,wt345,wt678,wt34,wt78,wt,
     & m3,m4,m5
      parameter(wt0=1._dp/twopi**5)

      oldzerowidth=zerowidth

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      p6(j)=0._dp
      p7(j)=0._dp
      enddo
      smin=(m3+m4)**2

      n2=0
      n3=0
      
c---generate p9 and p345678, 
c---smin is the minimum inv mass of 345678 system
c---m5 is the mass of p9

c--- DEBUG: make p9 soft
c      r(1)=1._dp-r(1)/1d3
c--- END DEBUG      
      call phi1_2m(m5,r(1),r(2),r(3),smin,p12,p9,p345678,wt12,*99)

      n2=1
      n3=1
      mass2=mt
      width2=twidth
      mass3=mt
      width3=twidth

c---decay 345678-system
      zerowidth=.true.
      call phi1_2(r(4),r(5),r(6),r(7),p345678,p345,p678,wt345678,*99)
      zerowidth=oldzerowidth

      n2=0
      n3=1
      mass3=wmass
      width3=wwidth

      smin=0._dp
c--decay of p345 into p5 and p34
      call phi1_2m(mb,r(8),r(9),r(10),smin,p345,p5,p34,wt345,*99)
c--decay of p678 into p6 and p78
      call phi1_2m(mb,r(11),r(12),r(13),smin,p678,p6,p78,wt678,*99)

      if ((p5(4)<=0._dp).or.(p6(4)<=0._dp)) goto 99
      call phi3m0(r(14),r(15),p34,p3,p4,wt34,*99)
      if ((p3(4)<=0._dp).or.(p4(4)<=0._dp)) goto 99
      call phi3m0(r(16),r(17),p78,p7,p8,wt78,*99)
      if ((p7(4)<=0._dp).or.(p8(4)<=0._dp)) goto 99

      wt=wt0*wt12*wt345678*wt345*wt678*wt34*wt78
      
      return
      
 99   continue
      zerowidth=oldzerowidth
      wt=0._dp

      return
      end

