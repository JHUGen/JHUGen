      block data codeversion_data
      
      include 'codeversion.f'
      data codeversion/'8.0'/      
      data      prelim/.false./      ! if true, print warning message
      end

      subroutine banner
      implicit none
      include 'types.f'
************************************************************************
*  Set the version number of MCFM and write out the banner heading     *
************************************************************************
      
      include 'codeversion.f'
      character*50 line
      integer:: vlength,lenocc

      line='**************************************************'
      vlength=lenocc(codeversion)
      vlength=vlength+17   
      line(25-vlength/2:24+(vlength+1)/2)=
     & ' MCFM - version '//codeversion//' '

      write(6,*) line

c--- warning message, if necessary
      if (prelim) then
        write(6,*) '*                                                *'
        write(6,*) '*           PRELIMINARY VERSION                  *'
        write(6,*) '*                                                *'
        write(6,*) '*  NOTE: This is a private release of the MCFM   *'
        write(6,*) '*  code that has not yet been made public on the *'
        write(6,*) '*  usual website. As such:                       *'
        write(6,*) '*                                                *'
        write(6,*) '*   + Please do not redistribute without the     *'
        write(6,*) '*     knowledge of the authors;                  *'
        write(6,*) '*                                                *'
        write(6,*) '*   + Please notify the authors of any bugs      *'
        write(6,*) '*     or problems so that they can be corrected  *'
        write(6,*) '*     before the next official release.          *'
        write(6,*) '*                                                *'
        write(6,*) line
      endif
     

      write(6,*) '*                                                *'
      write(6,*) '* MCFM, v'//codeversion//
     &                          '                  June 2nd, 2016  *'
      write(6,*) '*                                                *'
      write(6,*) '* Authors: John Campbell, Keith Ellis,           *'
      write(6,*) '*          Walter Giele, Ciaran Williams         *'
      write(6,*) '*         (johnmc@fnal.gov, ellis@fnal.gov,      *'
      write(6,*) '*          giele@fnal.gov,ciaranwi@buffalo.edu)  *'
      write(6,*) '*                                                *'
      write(6,*) '* For details see:                               *'
      write(6,*) '*                                                *'
      write(6,*) '*  Color singlet production at NNLO in MCFM      *'
      write(6,*) '*   R. Boughezal, J. Campbell, R.K. Ellis,       *'
      write(6,*) '*    C. Focke, W. Giele, X. Liu, F. Petriello,   *'
      write(6,*) '*    C. Williams,  arXiv:1605.08011              *'
      write(6,*) '*    (overview of NNLO implementation in MCFM)   *'
      write(6,*) '*                                                *'
      write(6,*) '*  arXiv:1603.02663 (diphotons at NNLO)          *'
      write(6,*) '*  arXiv:1601.00658 (VH at NNLO)                 *'
      write(6,*) '*                                                *'
      write(6,*) '*  arXiv:1502.02990 (VBF and VBS Higgs)          *'
      write(6,*) '*  arXiv:1403.2641  (Triphoton production)       *'
      write(6,*) '*  arXiv:1312.1628  (gg->WW, Higgs interference) *'
      write(6,*) '*  arXiv:1311.3589  (gg->ZZ, Higgs interference) *'
      write(6,*) '*  Phys.Rev.D87:114006, arXiv:1302.3856          *'
      write(6,*) '*  (tZ, tH -- with R. Rontsch)                   *'
      write(6,*) '*  arXiv:1211.6390 (DM, P. Fox and C. Williams)  *'
      write(6,*) '*  JHEP 1211:162 (2012), arXiv:1208.0566         *'
      write(6,*) '*  (Z+gam+jet,Z+gam+gam -- with H. Hartanto)     *'
      write(6,*) '*  arXiv:1204.1513 (top production+decay)        *'
      write(6,*) '*  JHEP 1207:052 (2012), arXiv:1204.5678 (ttW)   *'
      write(6,*) '*  JHEP 1110:005 (2011), arXiv:1107.5569         *'
      write(6,*) '*         (gg->WW,Higgs intference)              *'
      write(6,*) '*  JHEP 1107:018 (2011), arXiv:1105.0020         *'
      write(6,*) '*         (diboson update)                       *'
      write(6,*) '*  JHEP 1103:027 (2011), arXiv:1011.6647         *'
      write(6,*) '*         (Wbb for mb>0, with S. Badger)         *'
      write(6,*) '*  Phys.Rev.D81:074023, arXiv:1001.4495 (H+2jet) *'
      write(6,*) '*                                                *'
      write(6,*) '*  P.R.L. 102:142001, arXiv:0903.0005 [hep-ph]   *'
      write(6,*) '*    (t-channel single top + explicit b,         *'
      write(6,*) '*      JC, R.Frederix, F.Maltoni, F.Tramontano)  *'
      write(6,*) '*  N.P.B 726:109(2005), hep-ph/0506289 (W+t)     *'
      write(6,*) '*  Phys.Rev.D70:094012, hep-ph/0408158 (Sngl Top)*'
      write(6,*) '*       (with Francesco Tramontano)              *'
      write(6,*) '*                                                *'
      write(6,*) '*  Phys.Rev.D65:113007, hep-ph/0202176 (W,Z+2j)  *'
      write(6,*) '*  Phys.Rev.D62:114012, hep-ph/0006304 (W,Z+bb)  *'
      write(6,*) '*  Phys.Rev.D60:113006, hep-ph/9905386 (diboson) *'
      write(6,*) '*                                                *'
      write(6,*) '* On the web:  http://mcfm.fnal.gov/             *'
      write(6,*) '**************************************************'
      write(6,*) 
 
      return
      end








