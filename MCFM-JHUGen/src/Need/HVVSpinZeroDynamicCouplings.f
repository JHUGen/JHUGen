      subroutine HVVSpinZeroDynCoupl(res,ic,jh,sWp,sWm,sWW,tryWW)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'spinzerohiggs_anomcoupl.f'

      integer ic,jh
      double precision sWp,sWm,sWW
      double complex res,restmp
      double precision sWp_signed,sWm_signed
      double precision sWW_signed,Q2_012Factor
      double precision lambda_v,lambda_v120(1:3)
      double complex vvcoupl(1:8)
      logical tryWW
      logical forceZZcoupl
      logical doQ2_012_Coupl

      forceZZcoupl = .not.tryWW
      forceZZcoupl = forceZZcoupl.or.(distinguish_HWWcouplings.ne.1)
      forceZZcoupl = forceZZcoupl.or.(ic.gt.4)

      doQ2_012_Coupl = .false.
      sWp_signed=zip
      sWm_signed=zip
      sWW_signed=zip
      vvcoupl(:)=czip
      res=czip
      if(jh.eq.1) then

      if( forceZZcoupl ) then
         if(cz_q1sq.ne.0) then
            sWp_signed=abs(sWp)*dble(sign(1,cz_q1sq))
         endif
         if(cz_q2sq.ne.0) then
            sWm_signed=abs(sWm)*dble(sign(1,cz_q2sq))
         endif
         if(cz_q12sq.ne.0) then
            sWW_signed=abs(sWW)*dble(sign(1,cz_q12sq))
         endif
         doQ2_012_Coupl = doQ2_012_Coupl .or. cz_q1sq.ne.0
         doQ2_012_Coupl = doQ2_012_Coupl .or. cz_q2sq.ne.0
         doQ2_012_Coupl = doQ2_012_Coupl .or. cz_q12sq.ne.0
         if(ic.eq.1) then
            vvcoupl(1)=ghz1
            vvcoupl(2)=ghz1_prime
            vvcoupl(3)=ghz1_prime2
            vvcoupl(4)=ghz1_prime3
            vvcoupl(5)=ghz1_prime4
            vvcoupl(6)=ghz1_prime5
            vvcoupl(7)=ghz1_prime6
            vvcoupl(8)=ghz1_prime7
            lambda_v = Lambda_z1
            lambda_v120 = (/ Lambda_z11, Lambda_z12, Lambda_z10 /)
         elseif(ic.eq.2) then
            vvcoupl(1)=ghz2
            vvcoupl(2)=ghz2_prime
            vvcoupl(3)=ghz2_prime2
            vvcoupl(4)=ghz2_prime3
            vvcoupl(5)=ghz2_prime4
            vvcoupl(6)=ghz2_prime5
            vvcoupl(7)=ghz2_prime6
            vvcoupl(8)=ghz2_prime7
            lambda_v = Lambda_z2
            lambda_v120 = (/ Lambda_z21, Lambda_z22, Lambda_z20 /)
         elseif(ic.eq.3) then
            vvcoupl(1)=ghz3
            vvcoupl(2)=ghz3_prime
            vvcoupl(3)=ghz3_prime2
            vvcoupl(4)=ghz3_prime3
            vvcoupl(5)=ghz3_prime4
            vvcoupl(6)=ghz3_prime5
            vvcoupl(7)=ghz3_prime6
            vvcoupl(8)=ghz3_prime7
            lambda_v = Lambda_z3
            lambda_v120 = (/ Lambda_z31, Lambda_z32, Lambda_z30 /)
         elseif(ic.eq.4) then
            vvcoupl(1)=ghz4
            vvcoupl(2)=ghz4_prime
            vvcoupl(3)=ghz4_prime2
            vvcoupl(4)=ghz4_prime3
            vvcoupl(5)=ghz4_prime4
            vvcoupl(6)=ghz4_prime5
            vvcoupl(7)=ghz4_prime6
            vvcoupl(8)=ghz4_prime7
            lambda_v = Lambda_z4
            lambda_v120 = (/ Lambda_z41, Lambda_z42, Lambda_z40 /)
         elseif(ic.eq.5) then ! Zgs 1
            vvcoupl(3) = ghzgs1_prime2
            lambda_v = Lambda_zgs1
            lambda_v120 = (/ Lambda_z11, Lambda_z12, Lambda_z10 /)
         elseif(ic.eq.6) then ! Zgs 2-4
            vvcoupl(1) = ghzgs2
            lambda_v = 1d0 ! Not present
            lambda_v120 = (/ Lambda_z21, Lambda_z22, Lambda_z20 /)
         elseif(ic.eq.7) then
            vvcoupl(1) = ghzgs3
            lambda_v = 1d0 ! Not present
            lambda_v120 = (/ Lambda_z31, Lambda_z32, Lambda_z30 /)
         elseif(ic.eq.8) then
            vvcoupl(1) = ghzgs4
            lambda_v = 1d0 ! Not present
            lambda_v120 = (/ Lambda_z41, Lambda_z42, Lambda_z40 /)
         elseif(ic.eq.9) then ! gsgs 2-4
            vvcoupl(1) = ghgsgs2
            lambda_v = 1d0 ! Not present
            lambda_v120 = (/ Lambda_z21, Lambda_z22, Lambda_z20 /)
         elseif(ic.eq.10) then
            vvcoupl(1) = ghgsgs3
            lambda_v = 1d0 ! Not present
            lambda_v120 = (/ Lambda_z31, Lambda_z32, Lambda_z30 /)
         elseif(ic.eq.11) then
            vvcoupl(1) = ghgsgs4
            lambda_v = 1d0 ! Not present
            lambda_v120 = (/ Lambda_z41, Lambda_z42, Lambda_z40 /)
         endif
      else
         if(cw_q1sq.ne.0) then
            sWp_signed=abs(sWp)*dble(sign(1,cw_q1sq))
         endif
         if(cw_q2sq.ne.0) then
            sWm_signed=abs(sWm)*dble(sign(1,cw_q2sq))
         endif
         if(cw_q12sq.ne.0) then
            sWW_signed=abs(sWW)*dble(sign(1,cw_q12sq))
         endif
         doQ2_012_Coupl = doQ2_012_Coupl .or. cw_q1sq.ne.0
         doQ2_012_Coupl = doQ2_012_Coupl .or. cw_q2sq.ne.0
         doQ2_012_Coupl = doQ2_012_Coupl .or. cw_q12sq.ne.0
         if(ic.eq.1) then
            vvcoupl(1)=ghw1
            vvcoupl(2)=ghw1_prime
            vvcoupl(3)=ghw1_prime2
            vvcoupl(4)=ghw1_prime3
            vvcoupl(5)=ghw1_prime4
            vvcoupl(6)=ghw1_prime5
            vvcoupl(7)=ghw1_prime6
            vvcoupl(8)=ghw1_prime7
            lambda_v = Lambda_w1
            lambda_v120 = (/ Lambda_w11, Lambda_w12, Lambda_w10 /)
         elseif(ic.eq.2) then
            vvcoupl(1)=ghw2
            vvcoupl(2)=ghw2_prime
            vvcoupl(3)=ghw2_prime2
            vvcoupl(4)=ghw2_prime3
            vvcoupl(5)=ghw2_prime4
            vvcoupl(6)=ghw2_prime5
            vvcoupl(7)=ghw2_prime6
            vvcoupl(8)=ghw2_prime7
            lambda_v = Lambda_w2
            lambda_v120 = (/ Lambda_w21, Lambda_w22, Lambda_w20 /)
         elseif(ic.eq.3) then
            vvcoupl(1)=ghw3
            vvcoupl(2)=ghw3_prime
            vvcoupl(3)=ghw3_prime2
            vvcoupl(4)=ghw3_prime3
            vvcoupl(5)=ghw3_prime4
            vvcoupl(6)=ghw3_prime5
            vvcoupl(7)=ghw3_prime6
            vvcoupl(8)=ghw3_prime7
            lambda_v = Lambda_w3
            lambda_v120 = (/ Lambda_w31, Lambda_w32, Lambda_w30 /)
         elseif(ic.eq.4) then
            vvcoupl(1)=ghw4
            vvcoupl(2)=ghw4_prime
            vvcoupl(3)=ghw4_prime2
            vvcoupl(4)=ghw4_prime3
            vvcoupl(5)=ghw4_prime4
            vvcoupl(6)=ghw4_prime5
            vvcoupl(7)=ghw4_prime6
            vvcoupl(8)=ghw4_prime7
            lambda_v = Lambda_w4
            lambda_v120 = (/ Lambda_w41, Lambda_w42, Lambda_w40 /)
         endif
      endif

      if(vvcoupl(2).ne.czip) then
         restmp = vvcoupl(2)*lambda_v**4
         restmp = restmp/(lambda_v**2+abs(sWp))
         restmp = restmp/(lambda_v**2+abs(sWm))
         res = res+restmp
      endif
      if(vvcoupl(3).ne.czip) then
         res = res + vvcoupl(3) * ( sWp + sWm )/lambda_v**2
      endif
      if(vvcoupl(4).ne.czip) then
         res = res + vvcoupl(4) * ( sWp - sWm )/lambda_v**2
      endif
      if(vvcoupl(5).ne.czip) then
         res = res + vvcoupl(5) * ( sWW )/Lambda_Q**2
      endif
      if(vvcoupl(6).ne.czip) then
         res = res + vvcoupl(6) * ( sWp**2 + sWm**2 )/lambda_v**4
      endif
      if(vvcoupl(7).ne.czip) then
        res = res + vvcoupl(7) * ( sWp**2 - sWm**2 )/lambda_v**4
      endif
      if(vvcoupl(8).ne.czip) then
         res = res + vvcoupl(8) * ( sWp * sWm )/lambda_v**4
      endif

      if(ic.eq.1) then
         if(doQ2_012_Coupl) then
      Q2_012Factor = (lambda_v120(1)**2+sWp_signed)
      Q2_012Factor = Q2_012Factor*(lambda_v120(2)**2+sWm_signed)
      Q2_012Factor = Q2_012Factor*(lambda_v120(3)**2+sWW_signed)
            if(Q2_012Factor.ne.zip) then
      Q2_012Factor = 1d0/Q2_012Factor
      Q2_012Factor = Q2_012Factor*(lambda_v120(1))**2
      Q2_012Factor = Q2_012Factor*(lambda_v120(2))**2
      Q2_012Factor = Q2_012Factor*(lambda_v120(3))**2
            endif
            res = res * Q2_012Factor
         endif
         if(vvcoupl(1).ne.czip) then
            res = res + vvcoupl(1)
         endif
      else
         if(vvcoupl(1).ne.czip) then
            res = res + vvcoupl(1)
         endif
         if(doQ2_012_Coupl) then
      Q2_012Factor = (lambda_v120(1)**2+sWp_signed)
      Q2_012Factor = Q2_012Factor*(lambda_v120(2)**2+sWm_signed)
      Q2_012Factor = Q2_012Factor*(lambda_v120(3)**2+sWW_signed)
            if(Q2_012Factor.ne.zip) then
      Q2_012Factor = 1d0/Q2_012Factor
      Q2_012Factor = Q2_012Factor*(lambda_v120(1))**2
      Q2_012Factor = Q2_012Factor*(lambda_v120(2))**2
      Q2_012Factor = Q2_012Factor*(lambda_v120(3))**2
            endif
            res = res * Q2_012Factor
         endif
      endif

      elseif(jh.eq.2) then !-- HM Higgs

      if( forceZZcoupl ) then
         if(c2z_q1sq.ne.0) then
            sWp_signed=abs(sWp)*dble(sign(1,c2z_q1sq))
         endif
         if(c2z_q2sq.ne.0) then
            sWm_signed=abs(sWm)*dble(sign(1,c2z_q2sq))
         endif
         if(c2z_q12sq.ne.0) then
            sWW_signed=abs(sWW)*dble(sign(1,c2z_q12sq))
         endif
         doQ2_012_Coupl = doQ2_012_Coupl .or. c2z_q1sq.ne.0
         doQ2_012_Coupl = doQ2_012_Coupl .or. c2z_q2sq.ne.0
         doQ2_012_Coupl = doQ2_012_Coupl .or. c2z_q12sq.ne.0
         if(ic.eq.1) then
            vvcoupl(1)=gh2z1
            vvcoupl(2)=gh2z1_prime
            vvcoupl(3)=gh2z1_prime2
            vvcoupl(4)=gh2z1_prime3
            vvcoupl(5)=gh2z1_prime4
            vvcoupl(6)=gh2z1_prime5
            vvcoupl(7)=gh2z1_prime6
            vvcoupl(8)=gh2z1_prime7
            lambda_v = Lambda2_z1
            lambda_v120 = (/ Lambda2_z11, Lambda2_z12, Lambda2_z10 /)
         elseif(ic.eq.2) then
            vvcoupl(1)=gh2z2
            vvcoupl(2)=gh2z2_prime
            vvcoupl(3)=gh2z2_prime2
            vvcoupl(4)=gh2z2_prime3
            vvcoupl(5)=gh2z2_prime4
            vvcoupl(6)=gh2z2_prime5
            vvcoupl(7)=gh2z2_prime6
            vvcoupl(8)=gh2z2_prime7
            lambda_v = Lambda2_z2
            lambda_v120 = (/ Lambda2_z21, Lambda2_z22, Lambda2_z20 /)
         elseif(ic.eq.3) then
            vvcoupl(1)=gh2z3
            vvcoupl(2)=gh2z3_prime
            vvcoupl(3)=gh2z3_prime2
            vvcoupl(4)=gh2z3_prime3
            vvcoupl(5)=gh2z3_prime4
            vvcoupl(6)=gh2z3_prime5
            vvcoupl(7)=gh2z3_prime6
            vvcoupl(8)=gh2z3_prime7
            lambda_v = Lambda2_z3
            lambda_v120 = (/ Lambda2_z31, Lambda2_z32, Lambda2_z30 /)
         elseif(ic.eq.4) then
            vvcoupl(1)=gh2z4
            vvcoupl(2)=gh2z4_prime
            vvcoupl(3)=gh2z4_prime2
            vvcoupl(4)=gh2z4_prime3
            vvcoupl(5)=gh2z4_prime4
            vvcoupl(6)=gh2z4_prime5
            vvcoupl(7)=gh2z4_prime6
            vvcoupl(8)=gh2z4_prime7
            lambda_v = Lambda2_z4
            lambda_v120 = (/ Lambda2_z41, Lambda2_z42, Lambda2_z40 /)
         elseif(ic.eq.5) then ! Zgs 1
            vvcoupl(3) = gh2zgs1_prime2
            lambda_v = Lambda2_zgs1
            lambda_v120 = (/ Lambda2_z11, Lambda2_z12, Lambda2_z10 /)
         elseif(ic.eq.6) then ! Zgs 2-4
            vvcoupl(1) = gh2zgs2
            lambda_v = 1d0 ! Not present
            lambda_v120 = (/ Lambda2_z21, Lambda2_z22, Lambda2_z20 /)
         elseif(ic.eq.7) then
            vvcoupl(1) = gh2zgs3
            lambda_v = 1d0 ! Not present
            lambda_v120 = (/ Lambda2_z31, Lambda2_z32, Lambda2_z30 /)
         elseif(ic.eq.8) then
            vvcoupl(1) = gh2zgs4
            lambda_v = 1d0 ! Not present
            lambda_v120 = (/ Lambda2_z41, Lambda2_z42, Lambda2_z40 /)
         elseif(ic.eq.9) then ! gsgs 2-4
            vvcoupl(1) = gh2gsgs2
            lambda_v = 1d0 ! Not present
            lambda_v120 = (/ Lambda2_z21, Lambda2_z22, Lambda2_z20 /)
         elseif(ic.eq.10) then
            vvcoupl(1) = gh2gsgs3
            lambda_v = 1d0 ! Not present
            lambda_v120 = (/ Lambda2_z31, Lambda2_z32, Lambda2_z30 /)
         elseif(ic.eq.11) then
            vvcoupl(1) = gh2gsgs4
            lambda_v = 1d0 ! Not present
            lambda_v120 = (/ Lambda2_z41, Lambda2_z42, Lambda2_z40 /)
         endif
      else
         if(c2w_q1sq.ne.0) then
            sWp_signed=abs(sWp)*dble(sign(1,c2w_q1sq))
         endif
         if(c2w_q2sq.ne.0) then
            sWm_signed=abs(sWm)*dble(sign(1,c2w_q2sq))
         endif
         if(c2w_q12sq.ne.0) then
            sWW_signed=abs(sWW)*dble(sign(1,c2w_q12sq))
         endif
         doQ2_012_Coupl = doQ2_012_Coupl .or. c2w_q1sq.ne.0
         doQ2_012_Coupl = doQ2_012_Coupl .or. c2w_q2sq.ne.0
         doQ2_012_Coupl = doQ2_012_Coupl .or. c2w_q12sq.ne.0
         if(ic.eq.1) then
            vvcoupl(1)=gh2w1
            vvcoupl(2)=gh2w1_prime
            vvcoupl(3)=gh2w1_prime2
            vvcoupl(4)=gh2w1_prime3
            vvcoupl(5)=gh2w1_prime4
            vvcoupl(6)=gh2w1_prime5
            vvcoupl(7)=gh2w1_prime6
            vvcoupl(8)=gh2w1_prime7
            lambda_v = Lambda2_w1
            lambda_v120 = (/ Lambda2_w11, Lambda2_w12, Lambda2_w10 /)
         elseif(ic.eq.2) then
            vvcoupl(1)=gh2w2
            vvcoupl(2)=gh2w2_prime
            vvcoupl(3)=gh2w2_prime2
            vvcoupl(4)=gh2w2_prime3
            vvcoupl(5)=gh2w2_prime4
            vvcoupl(6)=gh2w2_prime5
            vvcoupl(7)=gh2w2_prime6
            vvcoupl(8)=gh2w2_prime7
            lambda_v = Lambda2_w2
            lambda_v120 = (/ Lambda2_w21, Lambda2_w22, Lambda2_w20 /)
         elseif(ic.eq.3) then
            vvcoupl(1)=gh2w3
            vvcoupl(2)=gh2w3_prime
            vvcoupl(3)=gh2w3_prime2
            vvcoupl(4)=gh2w3_prime3
            vvcoupl(5)=gh2w3_prime4
            vvcoupl(6)=gh2w3_prime5
            vvcoupl(7)=gh2w3_prime6
            vvcoupl(8)=gh2w3_prime7
            lambda_v = Lambda2_w3
            lambda_v120 = (/ Lambda2_w31, Lambda2_w32, Lambda2_w30 /)
         elseif(ic.eq.4) then
            vvcoupl(1)=gh2w4
            vvcoupl(2)=gh2w4_prime
            vvcoupl(3)=gh2w4_prime2
            vvcoupl(4)=gh2w4_prime3
            vvcoupl(5)=gh2w4_prime4
            vvcoupl(6)=gh2w4_prime5
            vvcoupl(7)=gh2w4_prime6
            vvcoupl(8)=gh2w4_prime7
            lambda_v = Lambda2_w4
            lambda_v120 = (/ Lambda2_w41, Lambda2_w42, Lambda2_w40 /)
         endif
      endif

      if(vvcoupl(2).ne.czip) then
         restmp = vvcoupl(2)*lambda_v**4
         restmp = restmp/(lambda_v**2+abs(sWp))
         restmp = restmp/(lambda_v**2+abs(sWm))
         res = res+restmp
      endif
      if(vvcoupl(3).ne.czip) then
         res = res + vvcoupl(3) * ( sWp + sWm )/lambda_v**2
      endif
      if(vvcoupl(4).ne.czip) then
         res = res + vvcoupl(4) * ( sWp - sWm )/lambda_v**2
      endif
      if(vvcoupl(5).ne.czip) then
         res = res + vvcoupl(5) * ( sWW )/Lambda2_Q**2
      endif
      if(vvcoupl(6).ne.czip) then
         res = res + vvcoupl(6) * ( sWp**2 + sWm**2 )/lambda_v**4
      endif
      if(vvcoupl(7).ne.czip) then
        res = res + vvcoupl(7) * ( sWp**2 - sWm**2 )/lambda_v**4
      endif
      if(vvcoupl(8).ne.czip) then
         res = res + vvcoupl(8) * ( sWp * sWm )/lambda_v**4
      endif

      if(ic.eq.1) then
         if(doQ2_012_Coupl) then
      Q2_012Factor = (lambda_v120(1)**2+sWp_signed)
      Q2_012Factor = Q2_012Factor*(lambda_v120(2)**2+sWm_signed)
      Q2_012Factor = Q2_012Factor*(lambda_v120(3)**2+sWW_signed)
            if(Q2_012Factor.ne.zip) then
      Q2_012Factor = 1d0/Q2_012Factor
      Q2_012Factor = Q2_012Factor*(lambda_v120(1))**2
      Q2_012Factor = Q2_012Factor*(lambda_v120(2))**2
      Q2_012Factor = Q2_012Factor*(lambda_v120(3))**2
            endif
            res = res * Q2_012Factor
         endif
         if(vvcoupl(1).ne.czip) then
            res = res + vvcoupl(1)
         endif
      else
         if(vvcoupl(1).ne.czip) then
            res = res + vvcoupl(1)
         endif
         if(doQ2_012_Coupl) then
      Q2_012Factor = (lambda_v120(1)**2+sWp_signed)
      Q2_012Factor = Q2_012Factor*(lambda_v120(2)**2+sWm_signed)
      Q2_012Factor = Q2_012Factor*(lambda_v120(3)**2+sWW_signed)
            if(Q2_012Factor.ne.zip) then
      Q2_012Factor = 1d0/Q2_012Factor
      Q2_012Factor = Q2_012Factor*(lambda_v120(1))**2
      Q2_012Factor = Q2_012Factor*(lambda_v120(2))**2
      Q2_012Factor = Q2_012Factor*(lambda_v120(3))**2
            endif
            res = res * Q2_012Factor
         endif
      endif

      else !-- No other High mass Higgses
         res = czip
      endif !-- End HM Higgs

      return
      end

      function checkHZAcouplings(jh)
      implicit none
      include 'constants.f'
      include 'spinzerohiggs_anomcoupl.f'
      logical checkHZAcouplings
      integer jh
      double complex sumchecked
      sumchecked = czip
      if(jh.eq.1) then
      sumchecked=sumchecked
     & + ghzgs1_prime2
     & + ghzgs2
     & + ghzgs3
     & + ghzgs4
      else if(jh.eq.2) then
      sumchecked=sumchecked
     & + gh2zgs1_prime2
     & + gh2zgs2
     & + gh2zgs3
     & + gh2zgs4
      endif
      checkHZAcouplings=(sumchecked.ne.czip)
      return
      end

      function checkHAAcouplings(jh)
      implicit none
      include 'constants.f'
      include 'spinzerohiggs_anomcoupl.f'
      logical checkHAAcouplings
      integer jh
      double complex sumchecked
      sumchecked = czip
      if(jh.eq.1) then
      sumchecked=sumchecked
     & + ghgsgs2
     & + ghgsgs3
     & + ghgsgs4
      else if(jh.eq.2) then
      sumchecked=sumchecked
     & + gh2gsgs2
     & + gh2gsgs3
     & + gh2gsgs4
      endif
      checkHAAcouplings=(sumchecked.ne.czip)
      return
      end
