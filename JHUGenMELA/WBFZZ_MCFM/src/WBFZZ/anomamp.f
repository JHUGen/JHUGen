      !-- H -> [l^-(i3) lb(i4)] [l^-(i5) lb(i6)]
      function anomzzamp(i3,i4,i5,i6,s3456,s34,s56,za,zb)
      implicit none
      include 'mxpart.f'
      include 'masses.f'
      include 'anom_higgs.f' 
      include 'spinzerohiggs_anomcoupl.f'
      include 'zprods_decl.f'
      double complex anomzzamp
      integer i3,i4,i5,i6
      double precision s3456,s34,s56
      double complex ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn
      double complex aa1,aa2,aa3
      double precision shat,q3_q3,q4_q4

      shat = abs(s3456)
      q3_q3 = abs(s34)
      q4_q4 = abs(s56)

c--- MARKUS: define q^2 dependent couplings

      ghz1_dyn = ghz1 
     &         + ghz1_prime * Lambda_z1**4/( Lambda_z1**2 
     &         + abs(q3_q3) )/(Lambda_z1**2 + abs(q4_q4))
     &         + ghz1_prime2*(abs(q3_q3)+abs(q4_q4))/Lambda_z1**2                                   
     &         + ghz1_prime3*(abs(q3_q3)-abs(q4_q4))/Lambda_z1**2                                                                                   
     &         + ghz1_prime4*(        shat         )/Lambda_Q**2                                                                                    
     &         + ghz1_prime5*(abs(q3_q3)**2+abs(q4_q4)**2)/Lambda_z1**4
     &         + ghz1_prime6*(abs(q3_q3)**2-abs(q4_q4)**2)/Lambda_z1**4
     &         + ghz1_prime7*(abs(q3_q3)*abs(q4_q4))/Lambda_z1**4



      ghz2_dyn = ghz2 
     &         + ghz2_prime * Lambda_z2**4/( Lambda_z2**2 
     &         + abs(q3_q3) )/(Lambda_z2**2 + abs(q4_q4))
     &         + ghz2_prime2*(abs(q3_q3)+abs(q4_q4))/Lambda_z2**2                                   
     &         + ghz2_prime3*(abs(q3_q3)-abs(q4_q4))/Lambda_z2**2                                                                                   
     &         + ghz2_prime4*(        shat         )/Lambda_Q**2                                                                                    
     &         + ghz2_prime5*(abs(q3_q3)**2+abs(q4_q4)**2)/Lambda_z2**4
     &         + ghz2_prime6*(abs(q3_q3)**2-abs(q4_q4)**2)/Lambda_z2**4
     &         + ghz2_prime7*(abs(q3_q3)*abs(q4_q4))/Lambda_z2**4


      ghz3_dyn = ghz3 
     &         + ghz3_prime * Lambda_z3**4/( Lambda_z3**2 
     &         + abs(q3_q3) )/(Lambda_z3**2 + abs(q4_q4))
     &         + ghz3_prime2*(abs(q3_q3)+abs(q4_q4))/Lambda_z3**2                                   
     &         + ghz3_prime3*(abs(q3_q3)-abs(q4_q4))/Lambda_z3**2                                                                                   
     &         + ghz3_prime4*(        shat         )/Lambda_Q**2                                                                                    
     &         + ghz3_prime5*(abs(q3_q3)**2+abs(q4_q4)**2)/Lambda_z3**4
     &         + ghz3_prime6*(abs(q3_q3)**2-abs(q4_q4)**2)/Lambda_z3**4
     &         + ghz3_prime7*(abs(q3_q3)*abs(q4_q4))/Lambda_z3**4


      ghz4_dyn = ghz4 
     &         + ghz4_prime * Lambda_z4**4/( Lambda_z4**2 
     &         + abs(q3_q3) )/(Lambda_z4**2 + abs(q4_q4))
     &         + ghz4_prime2*(abs(q3_q3)+abs(q4_q4))/Lambda_z4**2                                   
     &         + ghz4_prime3*(abs(q3_q3)-abs(q4_q4))/Lambda_z4**2                                                                                   
     &         + ghz4_prime4*(        shat         )/Lambda_Q**2                                                                                    
     &         + ghz4_prime5*(abs(q3_q3)**2+abs(q4_q4)**2)/Lambda_z4**4
     &         + ghz4_prime6*(abs(q3_q3)**2-abs(q4_q4)**2)/Lambda_z4**4
     &         + ghz4_prime7*(abs(q3_q3)*abs(q4_q4))/Lambda_z4**4


      aa1 =ghz1_dyn*zmass**2/shat
     &     + (shat-abs(q3_q3)-abs(q4_q4))/shat*
     &       (ghz2_dyn
     &       +ghz3_dyn*(shat-abs(q3_q3)-abs(q4_q4))/4d0/LambdaBSM**2)
      aa2 =-2d0*ghz2_dyn
     &     -ghz3_dyn*(shat-abs(q3_q3)-abs(q4_q4))/2d0/LambdaBSM**2
      aa3 =-2d0*ghz4_dyn


      aa1 = aa1 / zmass**2 !-- F
      aa2 = aa2 / zmass**2 !-- F
      aa3 = aa3 / zmass**2 !-- F

      anomzzamp = s3456 * aa1 * za(i3,i5)*zb(i6,i4) + 
     & aa2 * 0.5d0 * (za(i3,i5)**2*zb(i5,i4)*zb(i6,i3) + 
     &  za(i3,i5)*za(i4,i5)*zb(i5,i4)*zb(i6,i4) + 
     &  za(i3,i5)*za(i3,i6)*zb(i6,i3)*zb(i6,i4) + 
     &  za(i3,i6)*za(i4,i5)*zb(i6,i4)**2) +
     & aa3 * (0d0,0.5d0) * (za(i3,i4)*za(i5,i6)*zb(i6,i4)**2 -
     &  za(i3,i5)**2*zb(i4,i3)*zb(i6,i5))

      return

      end
      

      !-- H -> [l^-(i3) lb(i4)] [l^-(i5) lb(i6)]
      function anomwwamp(i3,i4,i5,i6,s3456,s34,s56,za,zb)
      implicit none
      include 'mxpart.f'
      include 'masses.f'
      include 'anom_higgs.f' 
      include 'spinzerohiggs_anomcoupl.f'
      include 'zprods_decl.f'
      double complex anomwwamp
      integer i3,i4,i5,i6
      double precision s3456,s34,s56
      double complex ghw1_dyn,ghw2_dyn,ghw3_dyn,ghw4_dyn
      double complex aa1,aa2,aa3
      double precision shat,q3_q3,q4_q4

      shat = abs(s3456)
      q3_q3 = abs(s34)
      q4_q4 = abs(s56)

c--- MARKUS: define q^2 dependent couplings

      ghw1_dyn = ghw1 
     &         + ghw1_prime * Lambda_z1**4/( Lambda_z1**2 
     &         + abs(q3_q3) )/(Lambda_z1**2 + abs(q4_q4))
     &         + ghw1_prime2*(abs(q3_q3)+abs(q4_q4))/Lambda_z1**2                                   
     &         + ghw1_prime3*(abs(q3_q3)-abs(q4_q4))/Lambda_z1**2                                                                                   
     &         + ghw1_prime4*(        shat         )/Lambda_Q**2                                                                                    
     &         + ghw1_prime5*(abs(q3_q3)**2+abs(q4_q4)**2)/Lambda_z1**4
     &         + ghw1_prime6*(abs(q3_q3)**2-abs(q4_q4)**2)/Lambda_z1**4
     &         + ghw1_prime7*(abs(q3_q3)*abs(q4_q4))/Lambda_z1**4



      ghw2_dyn = ghw2 
     &         + ghw2_prime * Lambda_z2**4/( Lambda_z2**2 
     &         + abs(q3_q3) )/(Lambda_z2**2 + abs(q4_q4))
     &         + ghw2_prime2*(abs(q3_q3)+abs(q4_q4))/Lambda_z2**2                                   
     &         + ghw2_prime3*(abs(q3_q3)-abs(q4_q4))/Lambda_z2**2                                                                                   
     &         + ghw2_prime4*(        shat         )/Lambda_Q**2                                                                                    
     &         + ghw2_prime5*(abs(q3_q3)**2+abs(q4_q4)**2)/Lambda_z2**4
     &         + ghw2_prime6*(abs(q3_q3)**2-abs(q4_q4)**2)/Lambda_z2**4
     &         + ghw2_prime7*(abs(q3_q3)*abs(q4_q4))/Lambda_z2**4


      ghw3_dyn = ghw3 
     &         + ghw3_prime * Lambda_z3**4/( Lambda_z3**2 
     &         + abs(q3_q3) )/(Lambda_z3**2 + abs(q4_q4))
     &         + ghw3_prime2*(abs(q3_q3)+abs(q4_q4))/Lambda_z3**2                                   
     &         + ghw3_prime3*(abs(q3_q3)-abs(q4_q4))/Lambda_z3**2                                                                                   
     &         + ghw3_prime4*(        shat         )/Lambda_Q**2                                                                                    
     &         + ghw3_prime5*(abs(q3_q3)**2+abs(q4_q4)**2)/Lambda_z3**4
     &         + ghw3_prime6*(abs(q3_q3)**2-abs(q4_q4)**2)/Lambda_z3**4
     &         + ghw3_prime7*(abs(q3_q3)*abs(q4_q4))/Lambda_z3**4


      ghw4_dyn = ghw4 
     &         + ghw4_prime * Lambda_z4**4/( Lambda_z4**2 
     &         + abs(q3_q3) )/(Lambda_z4**2 + abs(q4_q4))
     &         + ghw4_prime2*(abs(q3_q3)+abs(q4_q4))/Lambda_z4**2                                   
     &         + ghw4_prime3*(abs(q3_q3)-abs(q4_q4))/Lambda_z4**2                                                                                   
     &         + ghw4_prime4*(        shat         )/Lambda_Q**2                                                                                    
     &         + ghw4_prime5*(abs(q3_q3)**2+abs(q4_q4)**2)/Lambda_z4**4
     &         + ghw4_prime6*(abs(q3_q3)**2-abs(q4_q4)**2)/Lambda_z4**4
     &         + ghw4_prime7*(abs(q3_q3)*abs(q4_q4))/Lambda_z4**4


      aa1 =ghw1_dyn*wmass**2/shat
     &     + (shat-abs(q3_q3)-abs(q4_q4))/shat*
     &       (ghw2_dyn
     &       +ghw3_dyn*(shat-abs(q3_q3)-abs(q4_q4))/4d0/LambdaBSM**2)
      aa2 =-2d0*ghw2_dyn
     &     -ghw3_dyn*(shat-abs(q3_q3)-abs(q4_q4))/2d0/LambdaBSM**2
      aa3 =-2d0*ghw4_dyn


      aa1 = aa1 / wmass**2 !-- F
      aa2 = aa2 / wmass**2 !-- F
      aa3 = aa3 / wmass**2 !-- F

     
      anomwwamp = s3456 * aa1 * za(i3,i5)*zb(i6,i4) + 
     & aa2 * 0.5d0 * (za(i3,i5)**2*zb(i5,i4)*zb(i6,i3) + 
     &  za(i3,i5)*za(i4,i5)*zb(i5,i4)*zb(i6,i4) + 
     &  za(i3,i5)*za(i3,i6)*zb(i6,i3)*zb(i6,i4) + 
     &  za(i3,i6)*za(i4,i5)*zb(i6,i4)**2) +
     & aa3 * (0d0,0.5d0) * (za(i3,i4)*za(i5,i6)*zb(i6,i4)**2 -
     &  za(i3,i5)**2*zb(i4,i3)*zb(i6,i5))

      return

      end
      
