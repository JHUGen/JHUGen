      subroutine phase6(r,p1,p2,p3,p4,p5,p6,p7,p8,wt,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'zerowidth.f'
      include 'hdecaymode.f'
      include 'kprocess.f'
      include 'breit.f'
      include 'limits.f'
      include 'ipsgen.f'
c******* generate phase space for 2-->4 process
c******* r(mxdim),p1(4),p2(4) are inputs reversed in sign from physical values 
c---- phase space for -p1-p2 --> p5+p6+p3+p4+p7+p8
c---- with all 2 pi's (ie 1/(2*pi)^14)
      logical:: oldzerowidth,bw34_56
      real(dp):: r(mxdim)
      real(dp):: p1(4),p2(4),p5(4),p6(4),p3(4),p4(4),p7(4),p8(4)
      real(dp):: p12(4),p345(4),p678(4),p78(4),p34(4),smin,
     & p348(4),p567(4),p56(4),p3456(4),p5678(4),tmp(4)
      real(dp):: wt,wt0,wt12,wt678,wt345,wt34,wt78,wt3456,wt56,
     & wt5678,s34,s56,s36,s45,wt_zz,wt_ww
      integer:: j
      parameter(wt0=1._dp/twopi**4)
      oldzerowidth=zerowidth

      wt=0._dp
      do j=1,4
      p12(j)=-p1(j)-p2(j)
      enddo
      smin=mb**2

c---- calculate momenta of top and bbbar      
      n2=1
      n3=1
      if (
     &        (kcase==ktt_bbh)
     &   .or. (kcase==ktt_bbl)
     &   .or. (kcase==ktt_ldk)
     &   .or. (kcase==ktt_hdk)
     &   .or. (kcase==ktt_udk)
     &   .or. (kcase==ktthWdk)
     &   .or. (kcase==kttZbbl)
     &   .or. (kcase==ktt_bbu)
     &   .or. (kcase==kttbdkl)
     &   .or. (kcase==kttbdkh)
     &   .or. (kcase==kvlchk6)) then
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth
      zerowidth=.true.
        call phi1_2(r(1),r(2),r(3),r(4),p12,p345,p678,wt12,*99)
      zerowidth=oldzerowidth
      elseif ((kcase==kWtbwdk) .or. (kcase==kW_twdk)) then
        mass2=mt
        width2=twidth
        mass3=mt
        width3=twidth
        call phi1_2(r(1),r(2),r(3),r(4),p12,p348,p567,wt12,*99)
      elseif (kcase==ktautau) then
        mass2=mtau
        width2=tauwidth
        mass3=mtau
        width3=tauwidth
        zerowidth=.true.
        call phi1_2(r(1),r(2),r(3),r(4),p12,p345,p678,wt12,*99)
        zerowidth=oldzerowidth
      elseif ((kcase==kqqZZqq) .or. (kcase==kqqWWqq)) then
c--- Three different branches depending on cuts
        if (m3456max < 200._dp) then
c------ all events should be close to Higgs peak: use Higgs BW with 0.1 GeV width
          call phi1_2bw(r(1),r(2),r(3),r(4),
     &                  p12,p3456,p78,hmass,0.1_dp,wt12,*99)
        else
          if (m3456min > 135._dp) then
c------ range does not include peak so do not use BW form at all 
          call phi1_2nobw(r(1),r(2),r(3),r(4),p12,p3456,p78,wt12,*99)
          else
c------ entire range of m3456, including peak:  use BW at Higgs mass with 10 GeV width 
          call phi1_2bw(r(1),r(2),r(3),r(4),
     &                  p12,p3456,p78,hmass,10._dp,wt12,*99)
          endif
        endif
        call phi1_2(r(5),r(6),r(7),r(8),p3456,p56,p34,wt3456,*99)
        call phi3m0(r(13),r(14),p34,p3,p4,wt34,*99)
        call phi3m0(r(11),r(12),p56,p5,p6,wt56,*99)
        call phi3m0(r(15),r(16),p78,p7,p8,wt78,*99)
        wt=wt0*wt12*wt3456*wt34*wt56*wt78
c------ for qqWWqq, interchange 4<->6
        if (kcase==kqqWWqq) then
          tmp(:)=p4(:)
          p4(:)=p6(:)
          p6(:)=tmp(:)
        endif
        return
      elseif (kcase==kqqVVqq) then
c--- Three different branches depending on cuts
        if (m3456max < 200._dp) then
c------ all events should be close to Higgs peak: use Higgs BW with 0.1 GeV width
          call phi1_2bw(r(1),r(2),r(3),r(4),
     &                  p12,p3456,p78,hmass,0.1_dp,wt12,*99)
        else
          if (m3456min > 135._dp) then
c------ range does not include peak so do not use BW form at all 
          call phi1_2nobw(r(1),r(2),r(3),r(4),p12,p3456,p78,wt12,*99)
          else
c------ entire range of m3456, including peak:  use BW at Higgs mass with 10 GeV width 
          call phi1_2bw(r(1),r(2),r(3),r(4),
     &                  p12,p3456,p78,hmass,10._dp,wt12,*99)
          endif
        endif
c--- Further branching of p3456 depends on ipsgen
        if (ipsgen == 1) then
c--- generating Z(34) Z(56)
          bw34_56=.true.
          mass2=zmass
          mass3=zmass
          width2=zwidth
          width3=zwidth
        else
c--- generating W(36) W(54) --> exchange of p4 and p6 done later
          bw34_56=.false.
          mass2=wmass
          mass3=wmass
          width2=wwidth
          width3=wwidth
        endif
        call phi1_2(r(5),r(6),r(7),r(8),p3456,p56,p34,wt3456,*99)
        call phi3m0(r(13),r(14),p34,p3,p4,wt34,*99)
        call phi3m0(r(11),r(12),p56,p5,p6,wt56,*99)
        call phi3m0(r(15),r(16),p78,p7,p8,wt78,*99)
        wt=wt0*wt12*wt3456*wt34*wt56*wt78
c-----------------------------> now implement exchange of p4 and p6 if required
        if (bw34_56 .eqv. .false.) then
          tmp(:)=p4(:)
          p4(:)=p6(:)
          p6(:)=tmp(:)
        endif
c--- compute Jacobian factor to suppress poorly-sampled regions
        s34=2._dp*(p3(4)*p4(4)-p3(1)*p4(1)-p3(2)*p4(2)-p3(3)*p4(3))
        s56=2._dp*(p5(4)*p6(4)-p5(1)*p6(1)-p5(2)*p6(2)-p5(3)*p6(3))
        s36=2._dp*(p3(4)*p6(4)-p3(1)*p6(1)-p3(2)*p6(2)-p3(3)*p6(3))
        s45=2._dp*(p5(4)*p4(4)-p5(1)*p4(1)-p5(2)*p4(2)-p5(3)*p4(3))
c--- wt_zz must also suppress possible photon pole for Z(34), in addition to BW
        wt_zz=(((s34-zmass**2)**2+(zmass*zwidth)**2)
     &        *((s56-zmass**2)**2+(zmass*zwidth)**2))
     &        *(s34/zmass**2)**2
        wt_ww=(((s36-wmass**2)**2+(wmass*wwidth)**2)
     &        *((s45-wmass**2)**2+(wmass*wwidth)**2))
        if (bw34_56) then
          wt=wt*(wt_ww/(wt_zz+wt_ww))
        else
          wt=wt*(wt_zz/(wt_zz+wt_ww))
        endif
        return
      elseif ((kcase==kWpWp2j) .or. (kcase==kqqWWss)) then
        n2=0
        n3=0
        call phi1_2(r(1),r(2),r(3),r(4),p12,p3456,p78,wt12,*99)
        n2=1
        mass2=wmass
        width2=wwidth
        n3=1
        mass3=wmass
        width3=wwidth
        call phi1_2(r(5),r(6),r(7),r(8),p3456,p56,p34,wt3456,*99)
        call phi3m0(r(13),r(14),p34,p3,p4,wt34,*99)
        call phi3m0(r(11),r(12),p56,p5,p6,wt56,*99)
        call phi3m0(r(15),r(16),p78,p7,p8,wt78,*99)
        wt=wt0*wt12*wt3456*wt34*wt56*wt78
        return
      elseif ((kcase==kWpmZbj) .or. (kcase==kWpmZbb)
     &   .or. (kcase==kWpmZjj) .or. (kcase==kqqWZqq)) then
        n2=0
        n3=0
        call phi1_2(r(1),r(2),r(3),r(4),p12,p3456,p78,wt12,*99)
        n2=1
        mass2=zmass
        width2=zwidth
        n3=1
        mass3=wmass
        width3=wwidth
        call phi1_2(r(5),r(6),r(7),r(8),p3456,p56,p34,wt3456,*99)
        call phi3m0(r(13),r(14),p34,p3,p4,wt34,*99)
        call phi3m0(r(11),r(12),p56,p5,p6,wt56,*99)
        call phi3m0(r(15),r(16),p78,p7,p8,wt78,*99)
        wt=wt0*wt12*wt3456*wt34*wt56*wt78
      return
      elseif (kcase==kWH1jet) then
        n2=0
        n3=0
        call phi1_2(r(1),r(2),r(3),r(4),p12,p3456,p78,wt12,*99)
        n2=1
        mass2=hmass
        width2=hwidth
        n3=1
        mass3=wmass
        width3=wwidth
        call phi1_2(r(5),r(6),r(7),r(8),p3456,p56,p34,wt3456,*99)
        call phi3m0(r(13),r(14),p34,p3,p4,wt34,*99)


        if (hdecaymode == 'tlta') then
        call phi3m(r(11),r(12),p56,p5,p6,mtau,mtau,wt56,*99)
        elseif (hdecaymode == 'bqba') then
        call phi3m(r(11),r(12),p56,p5,p6,mb,mb,wt56,*99)
        elseif (hdecaymode == 'gaga') then
        call phi3m0(r(11),r(12),p56,p5,p6,wt56,*99)
        else
        write(6,*) 'Undefined decay option in phase6.f'
        write(6,*) 'hdecaymode',hdecaymode
        stop
        endif

        call phi3m0(r(15),r(16),p78,p7,p8,wt78,*99)
        wt=wt0*wt12*wt3456*wt34*wt56*wt78
      return

      elseif (kcase==kZH1jet) then
        n2=0
        n3=0
        call phi1_2(r(1),r(2),r(3),r(4),p12,p3456,p78,wt12,*99)
        n2=1
        mass2=hmass
        width2=hwidth
        n3=1
        mass3=zmass
        width3=zwidth
        call phi1_2(r(5),r(6),r(7),r(8),p3456,p56,p34,wt3456,*99)
        call phi3m0(r(13),r(14),p34,p3,p4,wt34,*99)

        if (hdecaymode == 'tlta') then
        call phi3m(r(11),r(12),p56,p5,p6,mtau,mtau,wt56,*99)
        elseif (hdecaymode == 'bqba') then
        call phi3m(r(11),r(12),p56,p5,p6,mb,mb,wt56,*99)
        elseif (hdecaymode == 'gaga') then
        call phi3m0(r(11),r(12),p56,p5,p6,wt56,*99)
        else
        write(6,*) 'Undefined decay option in phase6.f'
        write(6,*) 'hdecaymode',hdecaymode
        stop
        endif

        call phi3m0(r(15),r(16),p78,p7,p8,wt78,*99)
        wt=wt0*wt12*wt3456*wt34*wt56*wt78
      return


      elseif ((kcase==kHWWjet) 
     & .or. (kcase==kHWW2jt)
     & .or. (kcase==kqq_HWW)) then
c--- In the case of HWWjet, we should generate s3456 according to
c--- a Breit-Wigner at mH
        n2=1
        mass2=hmass
        width2=hwidth
        n3=0
        call phi1_2(r(1),r(2),r(3),r(4),p12,p3456,p78,wt12,*99)
        n2=1
        mass2=wmass
        width2=wwidth
        n3=1
        mass3=wmass
        width3=wwidth
        call phi1_2(r(5),r(6),r(7),r(8),p3456,p56,p34,wt3456,*99)
        call phi3m0(r(13),r(14),p34,p3,p4,wt34,*99)
        call phi3m0(r(11),r(12),p56,p5,p6,wt56,*99)
        call phi3m0(r(15),r(16),p78,p7,p8,wt78,*99)
        wt=wt0*wt12*wt3456*wt34*wt56*wt78
        return
      elseif ((kcase==kHZZjet)
     & .or. (kcase==kHZZ2jt)
     & .or. (kcase==kqq_HZZ)) then
c--- In the case of HZZjet, we should generate s3456 according to
c--- a Breit-Wigner at mH
        n2=1
        mass2=hmass
        width2=hwidth
        n3=0
        call phi1_2(r(1),r(2),r(3),r(4),p12,p3456,p78,wt12,*99)
        n2=1
        mass2=zmass
        width2=zwidth
        n3=1
        mass3=zmass
        width3=zwidth
        call phi1_2(r(5),r(6),r(7),r(8),p3456,p56,p34,wt3456,*99)
        call phi3m0(r(13),r(14),p34,p3,p4,wt34,*99)
        call phi3m0(r(11),r(12),p56,p5,p6,wt56,*99)
        call phi3m0(r(15),r(16),p78,p7,p8,wt78,*99)
        wt=wt0*wt12*wt3456*wt34*wt56*wt78
        return
      elseif ((kcase==kWH__WW) .or. (kcase==kZH__WW)) then
c--- In the case of W+H(->WW), we should generate s5678 according to
c--- a Breit-Wigner at mH
        n2=1
      if (kcase==kWH__WW) then
        mass2=wmass
        width2=wwidth
      else
        mass2=zmass
        width2=zwidth
      endif
        n3=1
        mass3=hmass
        width3=hwidth
        call phi1_2(r(1),r(2),r(3),r(4),p12,p34,p5678,wt12,*99)
        n2=1
        mass2=wmass
        width2=wwidth
        n3=1
        mass3=wmass
        width3=wwidth
        call phi1_2(r(5),r(6),r(7),r(8),p5678,p56,p78,wt5678,*99)
        call phi3m0(r(13),r(14),p34,p3,p4,wt34,*99)
        call phi3m0(r(11),r(12),p56,p5,p6,wt56,*99)
        call phi3m0(r(15),r(16),p78,p7,p8,wt78,*99)
        wt=wt0*wt12*wt5678*wt34*wt56*wt78
        return
      elseif ((kcase==kWH__ZZ) .or. (kcase==kZH__ZZ)) then
c--- In the case of W+H(->WW), we should generate s5678 according to
c--- a Breit-Wigner at mH
        n2=1
      if (kcase==kWH__ZZ) then
        mass2=wmass
        width2=wwidth
      else
        mass2=zmass
        width2=zwidth
      endif
        n3=1
        mass3=hmass
        width3=hwidth
        call phi1_2(r(1),r(2),r(3),r(4),p12,p34,p5678,wt12,*99)
        n2=1
        mass2=zmass
        width2=zwidth
        n3=1
        mass3=zmass
        width3=zwidth
        call phi1_2(r(5),r(6),r(7),r(8),p5678,p56,p78,wt5678,*99)
        call phi3m0(r(13),r(14),p34,p3,p4,wt34,*99)
        call phi3m0(r(11),r(12),p56,p5,p6,wt56,*99)
        call phi3m0(r(15),r(16),p78,p7,p8,wt78,*99)
        wt=wt0*wt12*wt5678*wt34*wt56*wt78
        return
      else
        write(*,*) 'Case not supported in phase6.f'
        stop
      endif

c      write(6,*) '345',(p345(4)**2-p345(3)**2-p345(2)**2-p345(1)**2)
c      write(6,*) '678',(p678(4)**2-p678(3)**2-p678(2)**2-p678(1)**2)

      mass3=wmass
      width3=wwidth
      if (   (kcase==ktt_bbh)  
     &  .or. (kcase==ktt_bbl)
     &  .or. (kcase==ktt_ldk)
     &  .or. (kcase==ktt_hdk)
     &  .or. (kcase==ktt_udk)
     &  .or. (kcase==ktthWdk)
     &  .or. (kcase==kttZbbl)
     &  .or. (kcase==ktt_bbu)
     &  .or. (kcase==kttbdkl)
     &  .or. (kcase==kttbdkh)
     &  .or. (kcase==kvlchk6)) then
        n3=1
        call phi1_2m(mb,r(5),r(6),r(7),smin,p345,p5,p34,wt345,*99)
        call phi1_2m(mb,r(8),r(11),r(12),smin,p678,p6,p78,wt678,*99)
      elseif ((kcase==kWtbwdk) .or. (kcase==kW_twdk)) then
        n3=1
        call phi1_2m(mb,r(5),r(6),r(7),smin,p348,p8,p34,wt345,*99)
        call phi1_2m(zip,r(8),r(11),r(12),smin,p567,p7,p56,wt678,*99)
      elseif (kcase==ktautau) then
        n3=0
        smin=0._dp
        call phi1_2m(smin,r(5),r(6),r(7),smin,p345,p5,p34,wt345,*99)
        call phi1_2m(smin,r(8),r(11),r(12),smin,p678,p6,p78,wt678,*99)
      endif

c      write(6,*) '34',(p34(4)**2-p34(3)**2-p34(2)**2-p34(1)**2)
c      write(6,*) '78',(p78(4)**2-p78(3)**2-p78(2)**2-p78(1)**2)

      call phi3m0(r(13),r(14),p34,p3,p4,wt34,*99)
      if   ((kcase==kWtbwdk) 
     & .or. (kcase==kW_twdk)
     & ) then
      call phi3m0(r(15),r(16),p56,p5,p6,wt78,*99)
      else
      call phi3m0(r(15),r(16),p78,p7,p8,wt78,*99)
      endif
      
      wt=wt0*wt12*wt345*wt678*wt34*wt78
      
      return
 99   wt=0._dp
      zerowidth=oldzerowidth
      return 1
      end

