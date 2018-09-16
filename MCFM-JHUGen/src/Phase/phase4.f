      subroutine phase4(r,p1,p2,p3,p4,p5,p6,wt,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'heavyflav.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'kprocess.f'
      include 'breit.f'
      include 'zerowidth.f'
      include 'limits.f'
c---- generate phase space for 2-->4 process
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign 
c---- from physical values 
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^8)
      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      real(dp):: p12(4),p34(4),p56(4),p345(4)
      real(dp):: wt,wt3456,wt34,wt56,wt0,mtbsq
      real(dp):: p35(4),p46(4),wt35,wt46
      integer:: j,n2save,n3save
      parameter(wt0=1._dp/twopi**2)

      n2save=n2
      n3save=n3

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
c p56 is the b-bbar system
c--- p3 and p4 are normally massless, except for single top
c--- production with an explicit b - then use mt (p3) and mb (p4)
      if ((kcase==kqg_tbq) .or. (kcase==kqq_tbg)
     ..or.(kcase==kqqtbgg) .or. (kcase==kqgtbqq)) then
c--- hack to generate small s56
c        r(1)=1._dp-r(1)/1d3
c        r(2)=r(2)/1d3
c--- end of hack

c--- alternative generation - to generate small s36 and s46
c        r(1)=r(1)/1d5
c        call phi1_2(r(1),r(2),r(3),r(4),p12,p34,p56,wt3456,*99)
c--- Next two lines for small s36
c        call phi3m(r(7),r(8),p34,p3,p6,mt,zip,wt34,*99)
c        call phi3m(r(5),r(6),p56,p4,p5,mb,zip,wt56,*99)
c--- Next two lines for small s46
c        call phi3m(r(7),r(8),p34,p4,p6,mb,zip,wt34,*99)
c        call phi3m(r(5),r(6),p56,p3,p5,mt,zip,wt56,*99)
c        wt=wt0*wt3456*wt34*wt56
c        return
c--- end hack      
c--- New-style PS generation
        mtbsq=(mt+mb)**2
        call phi1_2m(0._dp,r(1),r(2),r(3),mtbsq,p12,p6,p345,wt3456,*99)
        call phi1_2m(0._dp,r(4),r(5),r(6),mtbsq,p345,p5,p34,wt56,*99)
        call phi3m(r(7),r(8),p34,p3,p4,mt,mb,wt34,*99)
        wt=wt0*wt3456*wt34*wt56
        return

      elseif ((kcase==kH_tjet) .or. (kcase==kH_tdkj)) then
c--- New-style PS generation
        mtbsq=(mt+hmass)**2
        call phi1_2m(0._dp,r(1),r(2),r(3),mtbsq,p12,p6,p345,wt3456,*99)
        n2=0
        n3=1
        mtbsq=(hmass-20._dp*hwidth)**2
        call phi1_2m(mt,r(4),r(5),r(6),mtbsq,p345,p5,p34,wt56,*99)
        n2=n2save
        n3=n3save
        call phi3m(r(7),r(8),p34,p3,p4,mb,mb,wt34,*99)
        wt=wt0*wt3456*wt34*wt56
        return

      elseif ((kcase==kZ_tjet) .or. (kcase==kZ_tdkj)) then
c--- New-style PS generation
        if (zerowidth) then
          mtbsq=(mt+zmass)**2
        else
          mtbsq=(mt+sqrt(wsqmin))**2
        endif
        call phi1_2m(0._dp,r(1),r(2),r(3),mtbsq,p12,p6,p345,wt3456,*99)
        n2=0
        n3=1
        mtbsq=wsqmin
        call phi1_2m(mt,r(4),r(5),r(6),mtbsq,p345,p5,p34,wt56,*99)
        n2=n2save
        n3=n3save
        call phi3m(r(7),r(8),p34,p3,p4,mb,mb,wt34,*99)
        wt=wt0*wt3456*wt34*wt56
        return

c--- This branch is good for testing cancellation of IR singularities
c      elseif ((kcase=='qq_tbgIR') .or. (kcase=='qqtbggIR')) then
c        r(1)=r(1)*1.e-5_dp ! to check small s35
c        call phi1_2(r(1),r(2),r(3),r(4),p12,p35,p46,wt3456,*99)
c        call phi3m(r(5),r(6),p35,p3,p5,mt,0._dp,wt35,*99)
c        call phi3m(r(7),r(8),p46,p4,p6,mb,0._dp,wt46,*99)
c        wt=wt0*wt3456*wt35*wt46
c        return

c--- Default case
      else
        call phi1_2(r(1),r(2),r(3),r(4),p12,p56,p34,wt3456,*99)
        call phi3m0(r(7),r(8),p34,p3,p4,wt34,*99)
      endif
      
      if ( ((kcase==kWbbmas) .and. (flav == 5))
     ..or. ((kcase==kZbbmas) .and. (flav == 5))
     ..or. (kcase==kWHbbar)
     ..or. (kcase==kZHbbar)
     ..or. (kcase==kZccmas) .or. (kcase==kvlchkm)) then
        call phi3m(r(5),r(6),p56,p5,p6,mb,mb,wt56,*99)
      elseif ((kcase==kWbbmas) .and. (flav == 4)) then
        call phi3m(r(5),r(6),p56,p5,p6,mc,mc,wt56,*99)
      elseif (((kcase==kZbbmas) .and. (flav == 6))
     &    .or. (kcase==kqq_ttz) .or. (kcase==kqqtthz)) then
        call phi3m(r(5),r(6),p56,p5,p6,mt,mt,wt56,*99)
      elseif ((kcase==kqq_ttw) .or. (kcase==kttwldk) 
     &   .or. (kcase==kWttmas)) then
        call phi3m(r(5),r(6),p56,p5,p6,mt,mt,wt56,*99)
      elseif (kcase==kW_cjet) then
        call phi3m(r(5),r(6),p56,p5,p6,mc,zip,wt56,*99)
      elseif (kcase==kH_tjet) then
        call phi3m(r(5),r(6),p56,p5,p6,mt,zip,wt56,*99)
       
      elseif (kcase==kWbfrmc) then
        call phi3m(r(5),r(6),p56,p5,p6,mb,zip,wt56,*99)
      elseif (kcase==kW_tndk) then
        call phi3m(r(5),r(6),p56,p5,p6,mt,zip,wt56,*99)
      elseif (kcase==kWtbndk) then
        call phi3m(r(5),r(6),p56,p5,p6,mt,mb,wt56,*99)
      else
        call phi3m0(r(5),r(6),p56,p5,p6,wt56,*99)
      endif
      wt=wt0*wt3456*wt34*wt56
      if (debug) write(6,*) 'wt in phase4',wt
      return
      n2=n2save
      n3=n3save
 99   wt=0._dp
      return 1
      end

