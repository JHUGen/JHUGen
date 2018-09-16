c--- The variables R and P provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (Q1(a,b,c,is) and leg 2 (Q2(a,b,c,is)
c--- In each case the parton labelling is using the normal QM notation
c--- of putting everything backwards
c---       emitted line after emission =    a
c---       emitter before emission     =    b
c---       spectator                   =    c
c--- There is no label for he or she who is emitted.
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

      real(dp)::
     & R1(-1:1,-1:1,-1:1,0:2,3),R2(-1:1,-1:1,-1:1,0:2,3)
      common/RP_col_new/R1,R2
!$omp threadprivate(/RP_col_new/)
