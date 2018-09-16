      subroutine phase7dk(r,p1,p2,p3,p4,p5,p6,p7,p8,p9,wt,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'zerowidth.f'
      include 'plabel.f'
      include 'kprocess.f'
      include 'decay1q2a.f'
      include 'breit.f'
c******* generate phase space for 2-->7 process
c******* r(mxdim),p1(4),p2(4) are inputs reversed in sign from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6+p7+p8+p9
c---- with all 2 pi's (ie 1/(2*pi)^14)
      logical:: oldzerowidth
      integer:: j,nu
      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p5(4),p6(4),p3(4),p4(4),p7(4),p8(4),
     & p9(4)
      real(dp):: p12(4),p3459(4),p678(4),p78(4),p34(4),smin,
     & p59(4),p349(4)
      real(dp):: wt,wt0,wt12,wt678,wt3459,wt34,wt78,
     & wt59,wt349,tmp
      parameter(wt0=1._dp/twopi**5)

c--- alternate radiation between decay of top (=1) and antitop (=2) quarks;
c--- use additional (always uniform) variable to determine choice
      if (r(20) < 0.5_dp) then
        decay1q2a=1
      else
        decay1q2a=2
      endif

c--- DEBUG
c      decay1q2a=2 ! always radiation in anti-top decay
c      decay1q2a=1 ! DEBUG: always radiation in top decay

      oldzerowidth=zerowidth

      wt=0._dp
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      smin=zip

c---- calculate momenta of top and bbbar
      n2=1
      n3=1
      
c--- specific to ttbar including radiation in decay
      mass2=mt
      width2=twidth
      mass3=mt
      width3=twidth
      zerowidth=.true.
      call phi1_2(r(1),r(2),r(3),r(4),p12,p3459,p678,wt12,*99)
      zerowidth=oldzerowidth

      n2=0
      n3=1
      mass3=wmass
      width3=wwidth

      if (kcase==ktthWdk) then
c--- radiation from hadronic decay of W
        if (plabel(3) == 'pp') then
         decay1q2a=1
      else
         decay1q2a=2
      endif
c--- top decay including radiation
        n3=1
        call phi1_2m(mb,r(5),r(6),r(7),smin,p3459,p5,p349,wt3459,*99)
        n3=0
        call phi1_2m(zip,r(8),r(18),r(19),smin,p349,p9,p34,wt349,*99)
        call phi3m0(r(11),r(12),p34,p3,p4,wt34,*99)
c--- anti-top decay
        n3=1
        call phi1_2m(mb,r(13),r(14),r(15),smin,p678,p6,p78,wt678,*99)
        call phi3m0(r(16),r(17),p78,p7,p8,wt78,*99)      
        wt=wt0*wt12*wt3459*wt349*wt34*wt678*wt78
      else
c--- radiation from top (not W)
        call phi1_2(r(5),r(6),r(7),r(8),p3459,p59,p34,wt3459,*99)
        call phi3m(r(18),r(19),p59,p5,p9,mb,zip,wt59,*99)
        call phi3m0(r(11),r(12),p34,p3,p4,wt34,*99)
c--- anti-top decay
        call phi1_2m(mb,r(13),r(14),r(15),smin,p678,p6,p78,wt678,*99)
        call phi3m0(r(16),r(17),p78,p7,p8,wt78,*99)      
        wt=wt0*wt12*wt3459*wt34*wt59*wt678*wt78      
c--- multiply by a factor of two since we are including radiation
c--- from top and anti-top quarks at the same time
        wt=wt*2._dp
      endif
             
      if (decay1q2a == 2) then
c--- perform swap (3,4,5) -> (8,7,6)
      do nu=1,4
      tmp=p3(nu)
      p3(nu)=p8(nu)
      p8(nu)=tmp
      tmp=p4(nu)
      p4(nu)=p7(nu)
      p7(nu)=tmp
      tmp=p5(nu)
      p5(nu)=p6(nu)
      p6(nu)=tmp
      enddo
      endif
 
      return
 99   wt=0._dp
      zerowidth=oldzerowidth
      
      return 1
      end

