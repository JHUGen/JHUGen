      subroutine pdfwrap
      implicit none
      include 'types.f'
      include 'nlooprun.f'
      include 'pdlabel.f'
      include 'couple.f'
      character*40 TableFile
      character*100 gridname
      double precision cmass,bmass
      COMMON/QMASS/CMASS,BMASS

      if         (pdlabel .eq. 'mstw8lo') then
      amz=0.13939d0
      nlooprun=1
      elseif     (pdlabel .eq. 'mstw8nl') then
      amz=0.12018d0
      nlooprun=2
      elseif     (pdlabel .eq. 'mstw8nn') then
      amz=0.11707d0
      nlooprun=3
      elseif     (pdlabel .eq. 'MMHT_lo') then
      amz=0.135d0
      nlooprun=1
      elseif     (pdlabel .eq. 'MMHT_nl') then
      amz=0.120d0
      nlooprun=2
      elseif     (pdlabel .eq. 'MMHT_nn') then
      amz=0.118d0
      nlooprun=3
      elseif     (pdlabel .eq. 'mrs4nf3') then
      amz=0.1083d0
      nlooprun=2
      cmass=1000d0
      bmass=1001d0
      elseif     (pdlabel .eq. 'mrs4lf3') then
      amz=0.1186d0
      nlooprun=1
      cmass=1000d0
      bmass=1001d0
      elseif     (pdlabel .eq. 'mrs4nf4') then
      amz=0.1153d0
      nlooprun=2
      bmass=1000d0
      elseif     (pdlabel .eq. 'mrs4lf4') then
      amz=0.1251d0
      nlooprun=1
      bmass=1000d0
      elseif     (pdlabel .eq. 'mrs04nl') then
      amz=0.1205d0 
      nlooprun=2
      elseif     (pdlabel .eq. 'mrs04nn') then
      amz=0.1167d0
      nlooprun=3
      elseif     (pdlabel .eq. 'mrs02nl') then
      amz=0.1197d0 
      nlooprun=2
      elseif     (pdlabel .eq. 'mrs02nn') then
      amz=0.1154d0
      nlooprun=3
      elseif     (pdlabel .eq. 'mrs0119') then
      amz=0.119d0
      nlooprun=2
      elseif     (pdlabel .eq. 'mrs0117') then
      amz=0.117d0
      nlooprun=2
      elseif     (pdlabel .eq. 'mrs0121') then
      amz=0.121d0
      nlooprun=2
      elseif     (pdlabel .eq. 'mrs01_j') then
      amz=0.121d0
      nlooprun=2
      elseif     (pdlabel .eq. 'mrs01lo') then
      amz=0.130d0
      nlooprun=1

C MRS99
C  1     COR01  central gluon, a_s    300      0.1175   0.00537  C
C  2     COR02  higher gluon          300      0.1175   0.00497  C
C  3     COR03  lower gluon           300      0.1175   0.00398  C
C  4     COR04  lower a_s             229      0.1125   0.00585  C
C  5     COR05  higher a_s            383      0.1225   0.00384  C
C  6     COR06  quarks up             303.3    0.1178   0.00497  C
C  7     COR07  quarks down           290.3    0.1171   0.00593  C
C  8     COR08  strange up            300      0.1175   0.00524  C
C  9     COR09  strange down          300      0.1175   0.00524  C
C  10    C0R10  charm up              300      0.1175   0.00525  C
C  11    COR11  charm down            300      0.1175   0.00524  C
C  12    COR12  larger d/u            300      0.1175   0.00515  C

      elseif     (pdlabel .eq. 'mrs99_1') then
      amz=0.1175d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs99_2') then
      amz=0.1175d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs99_3') then
      amz=0.1175d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs99_4') then
      amz=0.1125d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs99_5') then
      amz=0.1225d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs99_6') then
      amz=0.1178d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs99_7') then
      amz=0.1171d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs99_8') then
      amz=0.1175d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs99_9') then
      amz=0.1175d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs9910') then
      amz=0.1175d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs9911') then
      amz=0.1175d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs9912') then
      amz=0.1175d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs98z1') then
      amz=0.1175d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs98z2') then
      amz=0.1175d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs98z3') then
      amz=0.1175d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs98z4') then
      amz=0.1125d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs98z5') then
      amz=0.1225d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs98ht') then
      amz=0.1170d0
      nlooprun=2
c      write(6,*) 'alpha_s(MZ) for mrs98ht has been modified from'
c      write(6,*) 'the inherent 0.1170 to a new value of 0.1175'    
      elseif (pdlabel .eq. 'mrs98l1') then
      amz=0.125d0
      nlooprun=1
      elseif (pdlabel .eq. 'mrs98l2') then
      amz=0.125d0
      nlooprun=1
      elseif (pdlabel .eq. 'mrs98l3') then
      amz=0.125d0
      nlooprun=1
      elseif (pdlabel .eq. 'mrs98l4') then
      amz=0.120d0
      nlooprun=1
      elseif (pdlabel .eq. 'mrs98l5') then
      amz=0.130d0
      nlooprun=1
C     TEMPORARY NAMING SCHEME:                                   C
C                                                                C
C  mode  set    comment             L(4)/MeV  a_s(M_Z)  grid#1   C
C  ----  ---    -------             --------  -------   ------   C
C                                                                C
C  1     LO05A  central gluon, a_s    174      0.1250   0.01518  C
C  2     LO09A  higher gluon          174      0.1250   0.01616  C
C  3     LO10A  lower gluon           174      0.1250   0.01533  C
C  4     LO01A  lower a_s             136      0.1200   0.01652  C
C  5     LO07A  higher a_s            216      0.1300   0.01522  C
C                                                                C
C                                                                C
C      The corresponding grid files are called lt05a.dat etc.    C
C                                                                C
      elseif (pdlabel .eq. 'mrs96r1') then
      amz=0.113d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs96r2') then
      amz=0.120d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs96r3') then
      amz=0.113d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs96r4') then
      amz=0.120d0
      nlooprun=2
      elseif (pdlabel .eq. 'hmrs90e') then
      amz=0.098382675d0
      nlooprun=2
      elseif (pdlabel .eq. 'hmrs90b') then
      amz=0.107961191d0
      amz=0.12801d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs95ap') then
      amz=0.112683043d0
      nlooprun=2
      elseif (pdlabel .eq. 'mrs95_g') then
      amz=0.114476658d0
      amz=0.13352d0
      nlooprun=2
c      amz=0.11297d0
C   1      CTEQ4M   Standard MSbar scheme   0.116        1.6      cteq4m.tbl
C   2      CTEQ4D   Standard DIS scheme     0.116        1.6      cteq4d.tbl
C   3      CTEQ4L   Leading Order           0.116        1.6      cteq4l.tbl
C   4      CTEQ4A1  Alpha_s series          0.110        1.6      cteq4a1.tbl
C   5      CTEQ4A2  Alpha_s series          0.113        1.6      cteq4a2.tbl
C   6      CTEQ4A3  same as CTEQ4M          0.116        1.6      cteq4m.tbl
C   7      CTEQ4A4  Alpha_s series          0.119        1.6      cteq4a4.tbl
C   8      CTEQ4A5  Alpha_s series          0.122        1.6      cteq4a5.tbl
C   9      CTEQ4HJ  High Jet                0.116        1.6      cteq4hj.tbl
C   10     CTEQ4LQ  Low Q0                  0.114        0.7      cteq4lq.tbl
      elseif (pdlabel .eq. 'cteq3_m') then
      amz=0.112d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq3_l') then
      amz=0.112d0
      nlooprun=1
      elseif (pdlabel .eq. 'cteq3_d') then
      amz=0.112d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq4_m') then
      amz=0.116d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq4_d') then
      amz=0.116d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq4_l') then 
      amz=0.132d0
      nlooprun=1
      elseif (pdlabel .eq. 'cteq4a1') then
      amz=0.110d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq4a2') then
      amz=0.113d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq4a3') then
      amz=0.116d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq4a4') then
      amz=0.119d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq4a5') then
      amz=0.122d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq4hj') then
      amz=0.116d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq4lq') then
      amz=0.114d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq5_m') then
      Call SetCtq5(1)
      amz=0.118d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq5_d') then
      Call SetCtq5(2)
      amz=0.118d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq5_l') then
      Call SetCtq5(3)
      amz=0.127d0
      nlooprun=1
      elseif (pdlabel .eq. 'cteq5l1') then
      amz=0.127d0
      nlooprun=1
      elseif (pdlabel .eq. 'cteq5hj') then
      Call SetCtq5(4)
      amz=0.118d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq5hq') then
      Call SetCtq5(5)
      amz=0.118d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq5f3') then
      Call SetCtq5(6)
      amz=0.106d0
      nlooprun=2
      cmass=1000d0
      bmass=1001d0
      elseif (pdlabel .eq. 'cteq5f4') then
      Call SetCtq5(7)
      amz=0.112d0
      nlooprun=2
      bmass=1000d0
      elseif (pdlabel .eq. 'cteq5m1') then
      Call SetCtq5(8)
      amz=0.118d0
      nlooprun=2
      elseif (pdlabel .eq. 'ctq5hq1') then
      Call SetCtq5(9)
      amz=0.118d0
      nlooprun=2
      elseif (pdlabel .eq. 'cteq66m') then
      amz=0.118d0
      Call SetCtq6(400)
      nlooprun=2
      elseif (pdlabel .eq. 'cteq61m') then
      amz=0.118d0
      Call SetCtq6(200)
      nlooprun=2
      elseif (pdlabel .eq. 'cteq6_m') then
      amz=0.118d0
      Call SetCtq6(1)
      nlooprun=2
      elseif (pdlabel .eq. 'cteq6_d') then
      amz=0.118d0
      Call SetCtq6(2)
      nlooprun=2
      elseif (pdlabel .eq. 'cteq6_l') then
      amz=0.118d0
      Call SetCtq6(3)
      nlooprun=1
      elseif (pdlabel .eq. 'cteq6l1') then
      amz=0.130d0
      Call SetCtq6(4)
      nlooprun=1
      elseif (pdlabel .eq. 'CT10.00') then
      amz=0.118d0
      call SetCT10(100)
      nlooprun=2
      elseif (pdlabel .eq. 'CT14.LL') then
      TableFile='CT14llo.pds'
      call SetCT14(TableFile)
      amz=0.130d0
      nlooprun=1
      elseif (pdlabel .eq. 'CT14.NL') then
      TableFile='CT14n.00.pds'
      call SetCT14(TableFile)
      amz=0.118d0
      nlooprun=2
      elseif (pdlabel .eq. 'CT14.NN') then
      TableFile='CT14nn.00.pds'
      call SetCT14(TableFile)
      amz=0.118d0
      nlooprun=3

      elseif (pdlabel .eq. 'NN2.3NL') then
      gridname='Pdfdata/NNPDF23_nlo_as_0118.LHgrid'
      amz=0.118d0
      nlooprun=2
      call NNPDFDriver(gridname)
      call NNinitPDF(0)
      elseif (pdlabel .eq. 'NN2.3NN') then
      gridname='Pdfdata/NNPDF23_nnlo_as_0118.LHgrid'
      amz=0.118d0
      nlooprun=3
      call NNPDFDriver(gridname)
      call NNinitPDF(0)

      elseif (pdlabel .eq. 'NN3.0LO') then
      gridname='Pdfdata/NNPDF30_lo_as_0118.LHgrid'
      amz=0.118d0
      nlooprun=1
      call NNPDFDriver(gridname)
      call NNinitPDF(0)
      elseif (pdlabel .eq. 'NN3.0NL') then
      gridname='Pdfdata/NNPDF30_nlo_as_0118.LHgrid'
      amz=0.118d0
      nlooprun=2
      call NNPDFDriver(gridname)
      call NNinitPDF(0)
      elseif (pdlabel .eq. 'NN3.0NN') then
      gridname='Pdfdata/NNPDF30_nnlo_as_0118.LHgrid'
      amz=0.118d0
      nlooprun=3
      call NNPDFDriver(gridname)
      call NNinitPDF(0)

c--- NEW ATTEMPT
      elseif (pdlabel .eq. 'mtungb1') then
c--- need a value here: Lambda = 200 MeV
      amz=0.109d0
      else
        write(6,*) 'Unimplemented distribution= ',pdlabel
        write(6,*) 'Recommended are: '
        write(6,*) '    mstw8lo      (LO)'
        write(6,*) '    mstw8nl      (NLO)'
        write(6,*) '    mstw8nn      (NNLO)'
        write(6,*) '    MMHT_lo      (LO)'
        write(6,*) '    MMHT_nl      (NLO)'
        write(6,*) '    MMHT_nn      (NNLO)'
        write(6,*) '    CT10.00      (NLO)'
        write(6,*) '    CT14.LL      (LO)'
        write(6,*) '    CT14.NL      (NLO)'
        write(6,*) '    CT14.NN      (NNLO)'
        write(6,*) '    NN2.3NL      (NLO)'
        write(6,*) '    NN2.3NN      (NNLO)'
        write(6,*) '    NN3.0LO      (LO)'
        write(6,*) '    NN3.0NL      (NLO)'
        write(6,*) '    NN3.0NN      (NNLO)' 
        write(6,*) 'Please refer to manual for further details.' 
!     .'mstw8lo,','mstw8nl,','mstw8nn,',
!     .'mrs4nf3,','mrs4lf3,','mrs4nf4,','mrs4lf4,',
!     .'mrs02nl,','mrs02nn,',
!     .'mrs0119,','mrs0117,','mrs0121,','mrs01_j,','mrs01lo,',
!     .'mrs99_1,','mrs99_2,','mrs99_3,','mrs99_4,','mrs99_5,','mrs99_6,',
!     .'mrs99_7,','mrs99_8,','mrs99_9,','mrs9910,','mrs9911,','mrs9912,',
!     .'mrs98z1,','mrs98z2,','mrs98z3,','mrs98z4,','mrs98z5,','mrs98ht,',
!     .'mrs98l1,','mrs98l2,','mrs98l3,','mrs98l4,','mrs98l5,',
!     .'mrs96r1,','mrs96r2,','mrs96r3,','mrs96r4,',
!     .'hmrs90e,','hmrs90b,',
!     .'mrs95ap,','mrs95_g,','mtungb1',
!     .'cteq3_m,','cteq3_l,','cteq3_d,',
!     .'cteq4_m,','cteq4_d,','cteq4_l,','cteq4a1,','cteq4a2,',
!     .'cteq4a3,','cteq4a4,','cteq4a5,','cteq4hj,','cteq4lq,',
!     .'cteq5_m,','cteq5_d,','cteq5_l,','cteq5hj,','cteq5hq,',
!     .'cteq5f3,','cteq5f4,','cteq5m1,','ctq5hq1,','cteq5l1,',
!     .'cteq6_m,','cteq6_d,','cteq6_l,','cteq6l1,',
!     .'cteq61m,','cteq66m,',
!     .'CT10.00,','CT14.LL,','CT14.NL,','CT14.NN,',
!     .'NN2.3NL,','NN2.3NN'

      stop
      endif
      return
      end
 

