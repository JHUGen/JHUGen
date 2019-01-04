      subroutine phase5(r,p1,p2,p3,p4,p5,p6,p7,wt)
c----phase space for signal
      implicit none
      include 'constants.f'
      include 'heavyflav.f'
      include 'masses.f'
      include 'process.f'
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

      double precision r(mxdim)
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),pswidth
      double precision p127(4),p12(4),p56(4),p34(4),p3456(4),p345(4)
      double precision wt,wt127,wt345,wt3456,wt34,wt56,wt0,smin,mtbsq
      integer j,n2save,n3save

      parameter(wt0=1d0/twopi**3)

      n2save=n2
      n3save=n3

      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo


      smin=mb**2

      if ( (case .eq. 'Z_tjet') .or. (case .eq. 'Zt2jet') 
     &.or. (case .eq. 'Z_tdkj') .or. (case .eq. 'Ztdk2j')) then
c--- New-style PS generation
        if (zerowidth) then
        mtbsq=(mt+zmass)**2
      else
        mtbsq=(mt+dsqrt(wsqmin))**2
      endif
        call phi1_2m(0d0,r(1),r(2),r(3),mtbsq,p12,p7,p3456,wt3456,*99)
        call phi1_2m(0d0,r(4),r(5),r(6),mtbsq,p3456,p6,p345,wt345,*99)
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
        call phi3m(r(12),r(13),p34,p3,p4,0d0,0d0,wt34,*99)
        wt=wt0*wt3456*wt345*wt56*wt34
        return
      endif


      if ( (case .eq. 'H_tjet') .or. (case .eq. 'Ht2jet') 
     &.or. (case .eq. 'H_tdkj')) then
c--- New-style PS generation
        if (zerowidth) then
        mtbsq=(mt+hmass)**2
      else
        mtbsq=(mt+dsqrt(wsqmin))**2
      endif
        call phi1_2m(0d0,r(1),r(2),r(3),mtbsq,p12,p7,p3456,wt3456,*99)
        call phi1_2m(0d0,r(4),r(5),r(6),mtbsq,p3456,p6,p345,wt345,*99)
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
        call phi3m(r(12),r(13),p34,p3,p4,0d0,0d0,wt34,*99)
        wt=wt0*wt3456*wt345*wt56*wt34
        return
      endif

c--- In the case of HVV_4l, we should generate s127 according to
c--- a Breit-Wigner at mH, otherwise just linearly      
      if (  (case .eq. 'HWW_4l') 
     . .or. (case .eq. 'HWW2lq')
     . .or. (case .eq. 'HWW_tb')
     . .or. (case .eq. 'HWWint')
     . .or. (case .eq. 'HZZ_4l')
     . .or. (case .eq. 'HZZ_tb')
     . .or. (case .eq. 'HZZint')
     . .or. (case .eq. 'HmZZ4l')
     . .or. (case .eq. 'HWWjet')
     . .or. (case .eq. 'HZZjet')
     . ) then
        call phi1_2m_bw(zip,r(13),r(12),r(11),smin,p12,p7,p127,
     .   hmass,hwidth,wt127,*99)
      elseif (case .eq. 'HZZqgI') then
c--- width to use in generation of PS: if too far from threshold, just
c--- use a width of 10 GeV in the B.W. to sample PS adequately
        if (hmass .lt. mass2+mass3-hwidth*5d0) then
          pswidth=10d0
        else
          pswidth=hwidth
        endif
        call phi1_2m_bw(zip,r(13),r(12),r(11),smin,p12,p7,p127,
     .   hmass,pswidth,wt127,*99)
      else
        call phi1_2m_nobw(zip,r(13),r(12),r(11),
     .   smin,p12,p7,p127,wt127,*99)
      endif
      
      call phi1_2(r(1),r(2),r(3),r(4),p127,p56,p34,wt3456,*99)

      if (   (case .eq. 'Wbbjet') .or. (case .eq. 'Wbbjem')
     .  .or. ((case .eq. 'Wbbmas') .and. (flav .eq. 5))  
     .  .or. (case .eq. 'WHbbar')  
     .  .or. (case .eq. 'ZHbbar')  
     .  .or. ((case .eq. 'W_bjet') .and. (mb .gt. 0d0)) ) then
        call phi3m(r(5),r(6),p56,p5,p6,mb,mb,wt56,*99)
      elseif ((case .eq. 'Wbbmas') .and. (flav .eq. 4)) then
        call phi3m(r(5),r(6),p56,p5,p6,mc,mc,wt56,*99)
      elseif ((case .eq. 'Wttmas') .or. (case .eq. 'qq_ttw')) then
        call phi3m(r(5),r(6),p56,p5,p6,mt,mt,wt56,*99)
      else
        call phi3m0(r(5),r(6),p56,p5,p6,wt56,*99)
      endif
      call phi3m0(r(7),r(8),p34,p3,p4,wt34,*99)
      wt=wt0*wt127*wt3456*wt56*wt34

  
      if (debug) write(6,*) 'wt127',wt127
      if (debug) write(6,*) 'wt3456',wt3456
      if (debug) write(6,*) 'wt34',wt34
      if (debug) write(6,*) 'wt56',wt56

      return
 99   continue
      wt=0d0
      n2=n2save
      n3=n3save
      return
      end

