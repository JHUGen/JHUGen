      real(dp),parameter::
     & zeta3=1.202056903159594285399738161511449990764986292_dp,
     & zeta2=pi**2/6._dp,
     & be0 = 11/three*CA-4/three*TR*NF,
     & be1= 34/three*CA**2-(20/three*CA+4*CF)*TR*NF,
     & Ga0 = four,
     & Ga1 = 4/three*((four-pisq)*CA+5*be0),
     & Gaq0 = Ga0*CF,
     & Gaq1 = Ga1*CF,
     & Gag0 = Ga0*CA,
     & Gag1 = Ga1*CA,
     & gBq0 = 6*CF,
     & gBq1 = CF*(CA*(146/nine-80*zeta3)+CF*(three-4*pisq+48*zeta3)
     &           +be0*(121/nine+2*pisq/three)),
     & gBg0 = 2*be0,
     & gBg1 = CA*(CA*(182/nine-32*zeta3)
     & +be0*(94/nine-2/three*pisq))+2*be1,
     & gams1=CA*(-64/nine+28*zeta3)+be0*(-56/nine+2*zeta2)
