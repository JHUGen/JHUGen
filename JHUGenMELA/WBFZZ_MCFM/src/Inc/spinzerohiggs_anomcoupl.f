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


      double precision LambdaBSM,Lambda_z1,Lambda_z2,Lambda_z3,Lambda_z4
      double precision Lambda_Q
      double complex ghz1,ghz2,ghz3,ghz4
      double complex ghz1_prime,ghz2_prime,ghz3_prime,ghz4_prime
      double complex ghz1_prime2,ghz2_prime2,ghz3_prime2,ghz4_prime2
      double complex ghz1_prime3,ghz2_prime3,ghz3_prime3,ghz4_prime3
      double complex ghz1_prime4,ghz2_prime4,ghz3_prime4,ghz4_prime4
      double complex ghz1_prime5,ghz2_prime5,ghz3_prime5,ghz4_prime5
      double complex ghz1_prime6,ghz2_prime6,ghz3_prime6,ghz4_prime6
      double complex ghz1_prime7,ghz2_prime7,ghz3_prime7,ghz4_prime7
      double precision Lambda_w1,Lambda_w2,Lambda_w3,Lambda_w4
      double complex ghw1,ghw2,ghw3,ghw4
      double complex ghw1_prime,ghw2_prime,ghw3_prime,ghw4_prime
      double complex ghw1_prime2,ghw2_prime2,ghw3_prime2,ghw4_prime2
      double complex ghw1_prime3,ghw2_prime3,ghw3_prime3,ghw4_prime3
      double complex ghw1_prime4,ghw2_prime4,ghw3_prime4,ghw4_prime4
      double complex ghw1_prime5,ghw2_prime5,ghw3_prime5,ghw4_prime5
      double complex ghw1_prime6,ghw2_prime6,ghw3_prime6,ghw4_prime6
      double complex ghw1_prime7,ghw2_prime7,ghw3_prime7,ghw4_prime7

      logical AllowAnomalousCouplings

      common/spinzerohiggs_anomcoupl/
     & LambdaBSM,Lambda_z1,Lambda_z2,Lambda_z3,Lambda_z4,Lambda_Q,
     & ghz1,ghz2,ghz3,ghz4,
     & ghz1_prime,ghz2_prime,ghz3_prime,ghz4_prime,
     & ghz1_prime2,ghz2_prime2,ghz3_prime2,ghz4_prime2,
     & ghz1_prime3,ghz2_prime3,ghz3_prime3,ghz4_prime3,
     & ghz1_prime4,ghz2_prime4,ghz3_prime4,ghz4_prime4,
     & ghz1_prime5,ghz2_prime5,ghz3_prime5,ghz4_prime5,
     & ghz1_prime6,ghz2_prime6,ghz3_prime6,ghz4_prime6,
     & ghz1_prime7,ghz2_prime7,ghz3_prime7,ghz4_prime7,
     & Lambda_w1,Lambda_w2,Lambda_w3,Lambda_w4,
     & ghw1,ghw2,ghw3,ghw4,
     & ghw1_prime,ghw2_prime,ghw3_prime,ghw4_prime,
     & ghw1_prime2,ghw2_prime2,ghw3_prime2,ghw4_prime2,
     & ghw1_prime3,ghw2_prime3,ghw3_prime3,ghw4_prime3,
     & ghw1_prime4,ghw2_prime4,ghw3_prime4,ghw4_prime4,
     & ghw1_prime5,ghw2_prime5,ghw3_prime5,ghw4_prime5,
     & ghw1_prime6,ghw2_prime6,ghw3_prime6,ghw4_prime6,
     & ghw1_prime7,ghw2_prime7,ghw3_prime7,ghw4_prime7,
     & AllowAnomalousCouplings
