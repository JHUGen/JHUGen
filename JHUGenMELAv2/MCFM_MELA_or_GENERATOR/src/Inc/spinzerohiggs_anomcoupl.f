c---- ANOMALOUS COUPLINGS AVAILABLE FOR SPIN 0 HIGGS
c---- g1: SM, g2: Higher-order scalar, g3: Form factor dependence on g2, g4: Pseudoscalar
c---- prime: Dipole ansatz product for q1,2**2 dependence
c---- prime 2: |q1**2| + |q2**2|
c---- prime 3: |q1**2| - |q2**2|
c---- prime 4: (q1 + q2)**2
c---- prime 5: q1**4 + q2**4
c---- prime 6: q1**4 - q2**4
c---- prime 7: |q1**2| * |q2**2|
c---- NOTE: Please add new future couplings in the same order for both declaration and the common block


      integer AllowAnomalousCouplings,distinguish_HWWcouplings

      integer AnomalCouplPR
      integer AnomalCouplDK
      integer channeltoggle_stu ! 0, 1, 2 for s, t+u and s+t+u
      integer vvhvvtoggle_vbfvh ! 0, 1, 2 for VBF, VH and VB+VH

      integer cz_q1sq,cz_q2sq,cz_q12sq
      integer cw_q1sq,cw_q2sq,cw_q12sq
      integer c2z_q1sq,c2z_q2sq,c2z_q12sq
      integer c2w_q1sq,c2w_q2sq,c2w_q12sq

      double precision LambdaBSM,Lambda_Q
      double precision Lambda_z1,Lambda_z2,Lambda_z3,Lambda_z4
      double precision Lambda_z11,Lambda_z21,Lambda_z31,Lambda_z41
      double precision Lambda_z12,Lambda_z22,Lambda_z32,Lambda_z42
      double precision Lambda_z10,Lambda_z20,Lambda_z30,Lambda_z40

      double complex ghz1,ghz2,ghz3,ghz4
      double complex ghz1_prime,ghz2_prime,ghz3_prime,ghz4_prime
      double complex ghz1_prime2,ghz2_prime2,ghz3_prime2,ghz4_prime2
      double complex ghz1_prime3,ghz2_prime3,ghz3_prime3,ghz4_prime3
      double complex ghz1_prime4,ghz2_prime4,ghz3_prime4,ghz4_prime4
      double complex ghz1_prime5,ghz2_prime5,ghz3_prime5,ghz4_prime5
      double complex ghz1_prime6,ghz2_prime6,ghz3_prime6,ghz4_prime6
      double complex ghz1_prime7,ghz2_prime7,ghz3_prime7,ghz4_prime7

      double precision Lambda_zgs1
      double complex ghzgs1_prime2,ghzgs2,ghzgs3,ghzgs4
      double complex ghgsgs2,ghgsgs3,ghgsgs4

      double precision Lambda_w1,Lambda_w2,Lambda_w3,Lambda_w4
      double precision Lambda_w11,Lambda_w21,Lambda_w31,Lambda_w41
      double precision Lambda_w12,Lambda_w22,Lambda_w32,Lambda_w42
      double precision Lambda_w10,Lambda_w20,Lambda_w30,Lambda_w40

      double complex ghw1,ghw2,ghw3,ghw4
      double complex ghw1_prime,ghw2_prime,ghw3_prime,ghw4_prime
      double complex ghw1_prime2,ghw2_prime2,ghw3_prime2,ghw4_prime2
      double complex ghw1_prime3,ghw2_prime3,ghw3_prime3,ghw4_prime3
      double complex ghw1_prime4,ghw2_prime4,ghw3_prime4,ghw4_prime4
      double complex ghw1_prime5,ghw2_prime5,ghw3_prime5,ghw4_prime5
      double complex ghw1_prime6,ghw2_prime6,ghw3_prime6,ghw4_prime6
      double complex ghw1_prime7,ghw2_prime7,ghw3_prime7,ghw4_prime7

      double precision h2mass,h2width

      double precision Lambda2BSM,Lambda2_Q
      double precision Lambda2_z1,Lambda2_z2,Lambda2_z3,Lambda2_z4
      double precision Lambda2_z11,Lambda2_z21,Lambda2_z31,Lambda2_z41
      double precision Lambda2_z12,Lambda2_z22,Lambda2_z32,Lambda2_z42
      double precision Lambda2_z10,Lambda2_z20,Lambda2_z30,Lambda2_z40

      double complex gh2z1,gh2z2,gh2z3,gh2z4
      double complex gh2z1_prime,gh2z2_prime,gh2z3_prime,gh2z4_prime
      double complex gh2z1_prime2,gh2z2_prime2,gh2z3_prime2,gh2z4_prime2
      double complex gh2z1_prime3,gh2z2_prime3,gh2z3_prime3,gh2z4_prime3
      double complex gh2z1_prime4,gh2z2_prime4,gh2z3_prime4,gh2z4_prime4
      double complex gh2z1_prime5,gh2z2_prime5,gh2z3_prime5,gh2z4_prime5
      double complex gh2z1_prime6,gh2z2_prime6,gh2z3_prime6,gh2z4_prime6
      double complex gh2z1_prime7,gh2z2_prime7,gh2z3_prime7,gh2z4_prime7

      double precision Lambda2_zgs1
      double complex gh2zgs1_prime2,gh2zgs2,gh2zgs3,gh2zgs4
      double complex gh2gsgs2,gh2gsgs3,gh2gsgs4

      double precision Lambda2_w1,Lambda2_w2,Lambda2_w3,Lambda2_w4
      double precision Lambda2_w11,Lambda2_w21,Lambda2_w31,Lambda2_w41
      double precision Lambda2_w12,Lambda2_w22,Lambda2_w32,Lambda2_w42
      double precision Lambda2_w10,Lambda2_w20,Lambda2_w30,Lambda2_w40

      double complex gh2w1,gh2w2,gh2w3,gh2w4
      double complex gh2w1_prime,gh2w2_prime,gh2w3_prime,gh2w4_prime
      double complex gh2w1_prime2,gh2w2_prime2,gh2w3_prime2,gh2w4_prime2
      double complex gh2w1_prime3,gh2w2_prime3,gh2w3_prime3,gh2w4_prime3
      double complex gh2w1_prime4,gh2w2_prime4,gh2w3_prime4,gh2w4_prime4
      double complex gh2w1_prime5,gh2w2_prime5,gh2w3_prime5,gh2w4_prime5
      double complex gh2w1_prime6,gh2w2_prime6,gh2w3_prime6,gh2w4_prime6
      double complex gh2w1_prime7,gh2w2_prime7,gh2w3_prime7,gh2w4_prime7

      common/spinzerohiggs_anomcoupl/
     & AllowAnomalousCouplings,
     & distinguish_HWWcouplings,
     & AnomalCouplPR,AnomalCouplDK,
     & channeltoggle_stu,vvhvvtoggle_vbfvh,
     & cz_q1sq,cz_q2sq,cz_q12sq,
     & cw_q1sq,cw_q2sq,cw_q12sq,
     & c2z_q1sq,c2z_q2sq,c2z_q12sq,
     & c2w_q1sq,c2w_q2sq,c2w_q12sq,
     & LambdaBSM,Lambda_Q,
     & Lambda_z1,Lambda_z2,Lambda_z3,Lambda_z4,
     & Lambda_z11,Lambda_z21,Lambda_z31,Lambda_z41,
     & Lambda_z12,Lambda_z22,Lambda_z32,Lambda_z42,
     & Lambda_z10,Lambda_z20,Lambda_z30,Lambda_z40,
     & ghz1,ghz2,ghz3,ghz4,
     & ghz1_prime,ghz2_prime,ghz3_prime,ghz4_prime,
     & ghz1_prime2,ghz2_prime2,ghz3_prime2,ghz4_prime2,
     & ghz1_prime3,ghz2_prime3,ghz3_prime3,ghz4_prime3,
     & ghz1_prime4,ghz2_prime4,ghz3_prime4,ghz4_prime4,
     & ghz1_prime5,ghz2_prime5,ghz3_prime5,ghz4_prime5,
     & ghz1_prime6,ghz2_prime6,ghz3_prime6,ghz4_prime6,
     & ghz1_prime7,ghz2_prime7,ghz3_prime7,ghz4_prime7,
     & Lambda_zgs1,
     & ghzgs1_prime2,ghzgs2,ghzgs3,ghzgs4,
     & ghgsgs2,ghgsgs3,ghgsgs4,
     & Lambda_w1,Lambda_w2,Lambda_w3,Lambda_w4,
     & Lambda_w11,Lambda_w21,Lambda_w31,Lambda_w41,
     & Lambda_w12,Lambda_w22,Lambda_w32,Lambda_w42,
     & Lambda_w10,Lambda_w20,Lambda_w30,Lambda_w40,
     & ghw1,ghw2,ghw3,ghw4,
     & ghw1_prime,ghw2_prime,ghw3_prime,ghw4_prime,
     & ghw1_prime2,ghw2_prime2,ghw3_prime2,ghw4_prime2,
     & ghw1_prime3,ghw2_prime3,ghw3_prime3,ghw4_prime3,
     & ghw1_prime4,ghw2_prime4,ghw3_prime4,ghw4_prime4,
     & ghw1_prime5,ghw2_prime5,ghw3_prime5,ghw4_prime5,
     & ghw1_prime6,ghw2_prime6,ghw3_prime6,ghw4_prime6,
     & ghw1_prime7,ghw2_prime7,ghw3_prime7,ghw4_prime7,
     & h2mass,h2width,
     & Lambda2BSM,Lambda2_Q,
     & Lambda2_z1,Lambda2_z2,Lambda2_z3,Lambda2_z4,
     & Lambda2_z11,Lambda2_z21,Lambda2_z31,Lambda2_z41,
     & Lambda2_z12,Lambda2_z22,Lambda2_z32,Lambda2_z42,
     & Lambda2_z10,Lambda2_z20,Lambda2_z30,Lambda2_z40,
     & gh2z1,gh2z2,gh2z3,gh2z4,
     & gh2z1_prime,gh2z2_prime,gh2z3_prime,gh2z4_prime,
     & gh2z1_prime2,gh2z2_prime2,gh2z3_prime2,gh2z4_prime2,
     & gh2z1_prime3,gh2z2_prime3,gh2z3_prime3,gh2z4_prime3,
     & gh2z1_prime4,gh2z2_prime4,gh2z3_prime4,gh2z4_prime4,
     & gh2z1_prime5,gh2z2_prime5,gh2z3_prime5,gh2z4_prime5,
     & gh2z1_prime6,gh2z2_prime6,gh2z3_prime6,gh2z4_prime6,
     & gh2z1_prime7,gh2z2_prime7,gh2z3_prime7,gh2z4_prime7,
     & Lambda2_zgs1,
     & gh2zgs1_prime2,gh2zgs2,gh2zgs3,gh2zgs4,
     & gh2gsgs2,gh2gsgs3,gh2gsgs4,
     & Lambda2_w1,Lambda2_w2,Lambda2_w3,Lambda2_w4,
     & Lambda2_w11,Lambda2_w21,Lambda2_w31,Lambda2_w41,
     & Lambda2_w12,Lambda2_w22,Lambda2_w32,Lambda2_w42,
     & Lambda2_w10,Lambda2_w20,Lambda2_w30,Lambda2_w40,
     & gh2w1,gh2w2,gh2w3,gh2w4,
     & gh2w1_prime,gh2w2_prime,gh2w3_prime,gh2w4_prime,
     & gh2w1_prime2,gh2w2_prime2,gh2w3_prime2,gh2w4_prime2,
     & gh2w1_prime3,gh2w2_prime3,gh2w3_prime3,gh2w4_prime3,
     & gh2w1_prime4,gh2w2_prime4,gh2w3_prime4,gh2w4_prime4,
     & gh2w1_prime5,gh2w2_prime5,gh2w3_prime5,gh2w4_prime5,
     & gh2w1_prime6,gh2w2_prime6,gh2w3_prime6,gh2w4_prime6,
     & gh2w1_prime7,gh2w2_prime7,gh2w3_prime7,gh2w4_prime7
