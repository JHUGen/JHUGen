      subroutine writereference
      implicit none
      include 'types.f'
      include 'kpart.f'
      include 'kprocess.f'
      include 'nproc.f'
      logical writerefs
      common/writerefs/writerefs

      if (writerefs .eqv. .false.) return

        write(6,*)
        write(6,58) '****************  MCFM references   ****************'
        write(6,58) '*                                                  *'
        write(6,58) '*  An update on vector boson pair production at    *'
        write(6,58) '*    hadron colliders                              *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, PRD60 (1999) 113006 *'
        write(6,58) '*                                                  *'
        write(6,58) '*  Vector boson pair production at the LHC         *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, C. Williams,        *'
        write(6,58) '*    JHEP 1107 (2011) 018                          *'
        write(6,58) '*                                                  *'
        write(6,58) '*  A multi-threaded version of MCFM                *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, W. Giele,           *'
        write(6,58) '*    EPJC 75 (2015) 246                            *'
        write(6,58) '*                                                  *'

      if (kpart == knnlo) then
!        write(6,58) '*  NNLO implementation in MCFM                     *'
!        write(6,58) '*   and H, W, Z processes:                         *'
!        write(6,58) '*                                                  *'
        write(6,58) '*  Color singlet production at NNLO in MCFM        *'
        write(6,58) '*   R. Boughezal, J. Campbell, R.K. Ellis,         *'
        write(6,58) '*    C. Focke, W. Giele, X. Liu, F. Petriello,     *'
        write(6,58) '*    C. Williams,  to appear                       *'
        write(6,58) '*                                                  *'
        if     (nproc == 285) then
        write(6,58) '*  Predictions for diphoton production at the      *'
        write(6,58) '*    LHC through NNLO in QCD                       *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, Ye Li, C. Williams  *'
        write(6,58) '*    arXiv:1603.02663                              *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 91) .and. (nproc <= 110)) then
        write(6,58) '*  Associated production of a Higgs boson at NNLO  *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, C. Williams,        *'
        write(6,58) '*    arXiv:1601.00658                              *'
        write(6,58) '*                                                  *'
        endif
      endif

        if     (nproc == 287) then
        write(6,58) '*  Triphoton production at hadron colliders        *'
        write(6,58) '*   J.M. Campbell, C. Williams, arXiv:1403.2641    *'
        write(6,58) '*                                                  *'
        elseif (nproc == 289) then
        write(6,58) '*  Four-photon production at the LHC: an           *'
        write(6,58) '*    of 2 -> 4 analytic unitarity                  *'
        write(6,58) '*   T. Dennen, C. Williams, arXiv:1411.3237        *'
        write(6,58) '*                                                  *'
        elseif ((nproc == 20) .or. (nproc == 25)) then
        write(6,58) '*  QCD corrections to the hadronic production of a *'
        write(6,58) '*    heavy quark pair including decay correlations *'
        write(6,58) '*   S. Badger, J.M. Campbell, R.K. Ellis           *'
        write(6,58) '*    arXiv:1011.6647                               *'
        write(6,58) '*                                                  *'
        elseif ((nproc == 21) .or. (nproc == 26)) then
        write(6,58) '*  Strong radiative correction to Wbb production   *'
        write(6,58) '*    in proton-antiproton collisions               *'
        write(6,58) '*   R.K. Ellis, S. Veseli, hep-ph/9810489          *'
        write(6,58) '*                                                  *'
        elseif ((nproc == 22) .or. (nproc == 27)
     &     .or. (nproc == 44) .or. (nproc == 46)) then
        write(6,58) '*  Next-to-leading order corrections to W+2jet     *'
        write(6,58) '*    and Z+2jet production at hadron colliders     *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, hep-ph/0202176      *'
        write(6,58) '*                                                  *'
        write(6,58) '*  Next-to-leading order QCD predictions for       *'
        write(6,58) '*    W+2jet and Z+2jet production at the CERN LHC  *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, D.L. Rainwater,     *'
        write(6,58) '*    hep-ph/0308195                                *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 50) .and. (nproc <= 56)) then
        write(6,58) '*  Radiative corrections to Zbb production         *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, hep-ph/0006304      *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 123) .and. (nproc <= 127)) then
        write(6,58) '*  Gluon-gluon contributions to W+ W- production   *'
        write(6,58) '*    and Higgs interference effects                *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, C. Williams,        *'
        write(6,58) '*    arXiv:1107.5569                               *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 128) .and. (nproc <= 132)) then
        write(6,58) '*  Bounding the Higgs width at the LHC using       *'
        write(6,58) '*    full analytic results for gg -> e-e+mu-mu+    *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, C. Williams,        *'
        write(6,58) '*    arXiv:1311.3589                               *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 141) .and. (nproc <= 151)) then
        write(6,58) '*  Top-quark processes at NLO in production        *'
        write(6,58) '*     and decay                                    *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, arXiv:1204.1513     *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 161) .and. (nproc <= 177)) then
        write(6,58) '*   Single top production and decay at             *'
        write(6,58) '*     next-to-leading order                        *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, F. Tramontano,      *'
        write(6,58) '*    hep-ph/0408158                                *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 180) .and. (nproc <= 187)) then
        write(6,58) '*  Next-to-leading order corrections to Wt         *'
        write(6,58) '*    production and decay                          *'
        write(6,58) '*   J.M. Campbell, F. Tramontano, hep-ph/0506289   *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 211) .and. (nproc <= 217)) then
        write(6,58) '*  Higgs boson production in weak boson fusion     *'
        write(6,58) '*    at next-to-leading order                      *'
        write(6,58) '*   E.L. Berger, J.M. Campbell, hep-ph/0403194     *'
        write(6,58) '*                                                  *'
        elseif (((nproc >= 220) .and. (nproc <= 229)) .or.
     &          ((nproc >= 2201) .and. (nproc <= 2291))) then
        write(6,58) '*  Higgs constraints from vector boson fusion      *'
        write(6,58) '*    and scattering                                *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, arXiv:1502.02990    *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 270) .and. (nproc <= 279)) then
        write(6,58) '*  Next-to-leading order Higgs + 2 jet production  *'
        write(6,58) '*    via gluon fusion                              *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, G. Zanderighi,      *'
        write(6,58) '*    hep-ph/0608194                                *'
        write(6,58) '*                                                  *'
        write(6,58) '*  Hadronic production of a Higgs boson and        *'
        write(6,58) '*    two jets at next-to-leading order             *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, C. Williams,        *'
        write(6,58) '*    arXiv:1001.4495                               *'
        write(6,58) '*                                                  *'
        elseif ((nproc == 301) .or. (nproc == 302)
     &     .or. (nproc == 306) .or. (nproc == 307)) then
        write(6,58) '*  Next-to-leading order predictions for           *'
        write(6,58) '*    Zgam+jet and Zgamgam final states at the LHC  *'
        write(6,58) '*   J.M. Campbell, H. Hartanto, C. Williams        *'
        write(6,58) '*    arXiv:1208.0566                               *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 500) .and. (nproc <= 516)) then
        write(6,58) '*  ttW production and decay at NLO                 *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, arXiv:1204.5678     *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 540) .and. (nproc <= 569)) then
        write(6,58) '*  Single top production in association with a     *'
        write(6,58) '*    Z boson at the LHC                            *'
        write(6,58) '*   J.M. Campbell, R.K. Ellis, R. Rontsch,         *'
        write(6,58) '*    arXiv:1302.3856                               *'
        write(6,58) '*                                                  *'
        elseif ((nproc >= 800) .and. (nproc <= 848)) then
        write(6,58) '*  Next-to-leading order predictions for           *'
        write(6,58) '*    dark matter production at hadron colliders    *'
        write(6,58) '*   P.J. Fox, C. Williams, arXiv:1211.6390         *'
        write(6,58) '*                                                  *'
        endif

        write(6,58) '****************************************************'

   58 format(a53)

      return
      end

