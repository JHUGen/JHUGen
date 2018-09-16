      subroutine phase5(r,p1,p2,p3,p4,p5,p6,p7,wt)
      implicit none
      include 'types.f'
c----phase space for signal
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'heavyflav.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'kprocess.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'breit.f'
      include 'zerowidth.f'
      include 'limits.f'
c********* generate phase space for 2-->5 process
c********* r(mxdim),p1(4),p2(4) are inputs 
c--------- incoming p1 and p2 reversed in sign from physical values 
c---- i.e. phase space for -p1-p2 --> p3+p4+p5+p6+p7
c---- with all 2 pi's (ie 1/(2*pi)^11)

      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),pswidth
      real(dp):: p127(4),p12(4),p56(4),p34(4),p3456(4),p345(4)
      real(dp):: wt,wt127,wt345,wt3456,wt34,wt56,wt0,smin,mtbsq
      integer:: j,n2save,n3save

      parameter(wt0=1._dp/twopi**3)

      n2save=n2
      n3save=n3

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo


      smin=mb**2

      if ( (kcase==kZ_tjet) .or. (kcase==kZt2jet) 
     &.or. (kcase==kZ_tdkj) .or. (kcase==kZtdk2j)) then
c--- New-style PS generation
        if (zerowidth) then
        mtbsq=(mt+zmass)**2
      else
        mtbsq=(mt+sqrt(wsqmin))**2
      endif
        call phi1_2m(0._dp,r(1),r(2),r(3),mtbsq,p12,p7,p3456,wt3456,*99)
        call phi1_2m(0._dp,r(4),r(5),r(6),mtbsq,p3456,p6,p345,wt345,*99)
        n2=0
        n3=1
        if (zerowidth) then
        mtbsq=zmass**2
      else
        mtbsq=wsqmin
      endif
        call phi1_2m(mt,r(7),r(8),r(11),mtbsq,p345,p5,p34,wt56,*99)
        n2=n2save
        n3=n3save
        call phi3m(r(12),r(13),p34,p3,p4,0._dp,0._dp,wt34,*99)
        wt=wt0*wt3456*wt345*wt56*wt34
        return
      endif


      if ( (kcase==kH_tjet) .or. (kcase==kHt2jet) 
     &.or. (kcase==kH_tdkj)) then
c--- New-style PS generation
        if (zerowidth) then
        mtbsq=(mt+hmass)**2
      else
        mtbsq=(mt+sqrt(wsqmin))**2
      endif
        call phi1_2m(0._dp,r(1),r(2),r(3),mtbsq,p12,p7,p3456,wt3456,*99)
        call phi1_2m(0._dp,r(4),r(5),r(6),mtbsq,p3456,p6,p345,wt345,*99)
        n2=0
        n3=1
        if (zerowidth) then
        mtbsq=hmass**2
      else
        mtbsq=wsqmin
      endif
        call phi1_2m(mt,r(7),r(8),r(11),mtbsq,p345,p5,p34,wt56,*99)
        n2=n2save
        n3=n3save
        call phi3m(r(12),r(13),p34,p3,p4,0._dp,0._dp,wt34,*99)
        wt=wt0*wt3456*wt345*wt56*wt34
        return
      endif

c--- In the case of HVV_4l, we should generate s127 according to
c--- a Breit-Wigner at mH, otherwise just linearly      
      if (  (kcase==kHWW_4l) 
     & .or. (kcase==kHWW2lq)
     & .or. (kcase==kHWW_tb)
     & .or. (kcase==kHWWint)
     & .or. (kcase==kHZZ_4l)
     & .or. (kcase==kHZZ_tb)
     & .or. (kcase==kHZZint)
     & .or. (kcase==kHmZZ4l)
     & .or. (kcase==kHWWjet)
     & .or. (kcase==kHZZjet)
     & ) then
        call phi1_2m_bw(zip,r(13),r(12),r(11),smin,p12,p7,p127,
     &   hmass,hwidth,wt127,*99)
      elseif (kcase==kHZZqgI) then
c--- width to use in generation of PS: if too far from threshold, just
c--- use a width of 10 GeV in the B.W. to sample PS adequately
        if (hmass < mass2+mass3-hwidth*5._dp) then
          pswidth=10._dp
        else
          pswidth=hwidth
        endif
        call phi1_2m_bw(zip,r(13),r(12),r(11),smin,p12,p7,p127,
     &   hmass,pswidth,wt127,*99)
      else
        call phi1_2m_nobw(zip,r(13),r(12),r(11),
     &   smin,p12,p7,p127,wt127,*99)
      endif
      
      call phi1_2(r(1),r(2),r(3),r(4),p127,p56,p34,wt3456,*99)
      if (   (kcase==kWbbjet) .or. (kcase==kWbbjem)
     &  .or. ((kcase==kWbbmas) .and. (flav == 5))  
     &  .or. (kcase==kWHbbar)  
     &  .or. (kcase==kZHbbar)  
     &  .or. ((kcase==kW_bjet) .and. (mb > 0._dp)) ) then
        call phi3m(r(5),r(6),p56,p5,p6,mb,mb,wt56,*99)
      elseif ((kcase==kWbbmas) .and. (flav == 4)) then
        call phi3m(r(5),r(6),p56,p5,p6,mc,mc,wt56,*99)
      elseif ((kcase==kWttmas) .or. (kcase==kqq_ttw)) then
        call phi3m(r(5),r(6),p56,p5,p6,mt,mt,wt56,*99)
      elseif (((kcase==kWH1jet) .or. (kcase==kZH1jet)) 
     &  .and. (hdecaymode == 'bqba')) then
        call phi3m(r(5),r(6),p56,p5,p6,mb,mb,wt56,*99)
      elseif (((kcase==kWH1jet) .or. (kcase==kZH1jet)) 
     &  .and. (hdecaymode == 'tlta')) then
        call phi3m(r(5),r(6),p56,p5,p6,mtau,mtau,wt56,*99)
        elseif (((kcase==kWH1jet) .or. (kcase==kZH1jet)) 
     &  .and. (hdecaymode == 'gaga')) then
        call phi3m0(r(5),r(6),p56,p5,p6,wt56,*99)
      else
        call phi3m0(r(5),r(6),p56,p5,p6,wt56,*99)
      endif
      call phi3m0(r(7),r(8),p34,p3,p4,wt34,*99)
      wt=wt0*wt127*wt3456*wt56*wt34

  
c      if (debug) write(6,*) 'wt127',wt127
c      if (debug) write(6,*) 'wt3456',wt3456
c      if (debug) write(6,*) 'wt34',wt34
c      if (debug) write(6,*) 'wt56',wt56

      return
 99   continue
      wt=0._dp
      n2=n2save
      n3=n3save
      return
      end

