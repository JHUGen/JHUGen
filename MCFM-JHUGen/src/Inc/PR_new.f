c--- The variables R and P provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (Q1(a,b,c,is) and leg 2 (Q2(a,b,c,is)
c--- In each case the parton labelling is Using the normal QM notation of putting
c--- everything backward
c---       emitted line after emission =   a
c---       emitter before emission     =    b
c---       spectator                   =    c
c--- There is no label for he or she who is emitted.
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

      real(dp)::
     & Q1(-1:1,-1:1,-1:1,3),Q2(-1:1,-1:1,-1:1,3)
      common/RP_new/Q1,Q2
!$omp threadprivate(/RP_new/)
