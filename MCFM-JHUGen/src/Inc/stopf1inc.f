      complex(dp):: lc,ls,lVc,lVs,lp
     &,LsA,LsB1,LsB2,tr1Xfc,tr1Xfs,tr2fu,tr3Xc
     &,tr3c00fs,tr3c001fs,tr3c002fs,tr3Xs,tr3s00ft,tr3s001ft
     &,tr3s002ft,tr4Xc,tr4Xs,tr5Xc,tr5Xs,B0csf,B0cgsf
     &,lRc1,lRc2,lRs1,lRs2,lRcs,BfunX
      common/stopf1inc/ lc,ls,lVc,lVs,lp
     &,LsA,LsB1,LsB2,tr1Xfc,tr1Xfs,tr2fu,tr3Xc
     &,tr3c00fs,tr3c001fs,tr3c002fs,tr3Xs,tr3s00ft,tr3s001ft
     &,tr3s002ft,tr4Xc,tr4Xs,tr5Xc,tr5Xs,B0csf,B0cgsf
     &,lRc1,lRc2,lRs1,lRs2,lRcs,BfunX
!$omp threadprivate(/stopf1inc/)
      include 'sck.f'

