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

c--- Additional label (1-8) in this array, to represent the partons
c--- in the final state, according to the look-up parameters
c--- that are also defined here

c--- code is as follows:
c---  g=0, q=1, a=-1, r=2, b=-2 --- "f" for final
      integer:: gf_gf,qf_af,qf_qf,qf_rf,bf_rf,rf_bf,qf_gf,af_gf
      parameter (gf_gf=1,qf_af=2,qf_qf=3,qf_rf=4,
     & bf_rf=5,rf_bf=6,qf_gf=7,af_gf=8)
      real(dp):: 
     & S1(-1:1,-1:1,-1:1,8,0:2,3),S2(-1:1,-1:1,-1:1,8,0:2,3)
      common/SP_twojet/S1,S2
!$omp threadprivate(/SP_twojet/)
