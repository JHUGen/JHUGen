      function fun4(tcs,qsq,msq)
      implicit none
      include 'types.f'
      complex(dp):: fun4
      
      real(dp):: tcs,qsq,msq,qsqhat,taucs,qcs
      complex(dp):: lnrat,C0fa2m

      qsqhat=qsq-msq
      taucs=tcs-msq
      qcs=qsq-tcs

      fun4=  + qcs**(-1) + lnrat( - taucs,msq)*taucs*qcs**(-2) -
     &    lnrat( - qsqhat,msq)*tcs*qsq**(-1)*qsqhat*qcs**(-2) + C0fa2m(
     &    tcs,qsq,msq)*msq*qcs**(-1)

      return
      end


      function fun3(tcg,qsq,msq)
      implicit none
      include 'types.f'
      complex(dp):: fun3
      
      real(dp):: tcg,qsq,msq,qsqhat,taucg,qcg
      complex(dp):: lnrat

      qsqhat=qsq-msq
      taucg=tcg-msq
      qcg=qsq-tcg

      fun3=  + taucg*qcg**(-1) + lnrat( - taucg,msq)*taucg**2*
     &    qcg**(-2) - lnrat( - qsqhat,msq)*taucg*tcg*qsq**(-1)*qsqhat*
     &    qcg**(-2) - lnrat( - qsqhat,msq)*taucg*qsq**(-1)*qsqhat*
     &    qcg**(-1) + lnrat( - qsqhat,msq)*tcg*qsq**(-1)*qsqhat*
     &    qcg**(-1)

      return
      end


      function fun2(tx,ty,msq)
      implicit none
      include 'types.f'
      complex(dp):: fun2
      
      real(dp):: tx,ty,msq,taux,tauy,tymtx
      complex(dp):: lnrat

      taux=tx-msq
      tauy=ty-msq
      tymtx=ty-tx

      fun2=taux*lnrat(-taux,msq)/tymtx-tx*tauy*lnrat(-tauy,msq)/ty/tymtx

      return
      end


      function fun1(tcg,qsq,msq)
      implicit none
      include 'types.f'
      complex(dp):: fun1
      
      real(dp):: tcg,qsq,msq,qsqhat,taucg,qcg
      complex(dp):: lnrat

      qsqhat=qsq-msq
      taucg=tcg-msq
      qcg=qsq-tcg

      fun1=  + qsq*qcg**(-1) + lnrat( - taucg,msq)*taucg*tcg**(-1)*
     &    qsq**2*qcg**(-2) + lnrat( - qsqhat,msq)*taucg**(-1)*qsq*
     &    qsqhat*qcg**(-1) - lnrat( - qsqhat,msq)*qsq*qsqhat*qcg**(-2)
     &     - lnrat( - qsqhat,msq)*qsqhat*qcg**(-1)

      return
      end


      function L6m1(taucg,taucs,taugs,msq,qsq)
      implicit none
      include 'types.f'
      complex(dp):: L6m1
      
      real(dp):: taucg,taucs,taugs,msq,qsq
      complex(dp):: I3me,Lsm1_2m

      L6m1=  + Lsm1_2m(taugs,taucs,qsq,msq) + 2.D0*I3me(msq,taugs,qsq)*
     & taucs**(-1)*taugs*msq - I3me(msq,taugs,qsq)*taucs -
     & I3me(msq,taugs,qsq)*taucg

      return
      end


      function L6m2(taucg,taucs,taugs,msq,qsq)
      implicit none
      include 'types.f'
      complex(dp):: L6m2
      
      real(dp):: taucg,taucs,taugs,msq,qsq,tcg,tcs
      complex(dp):: C0fa2m,C0fb2m,Lsm2_2m
      tcs=taucs+msq
      tcg=taucg+msq

      L6m2=  + Lsm2_2m(taucs,taucg,qsq,msq) + 2.D0*C0fb2m(tcg,msq)*
     & taucs**(-1)*taugs*msq - C0fb2m(tcg,msq)*taucg - 2.D0*
     & C0fa2m(tcs,qsq,msq)*taucs**(-1)*taucg**(-1)*taugs**2*msq -
     & 2.D0*C0fa2m(tcs,qsq,msq)*taucs**(-1)*taugs*msq +
     & C0fa2m(tcs,qsq,msq)*taucg + C0fa2m(tcs,qsq,msq)*taugs

      return
      end


      function L6m3(taucg,taucs,taugs,msq,qsq)
      implicit none
      include 'types.f'
      complex(dp):: L6m3
      
      real(dp):: taucg,taucs,taugs,msq,qsq
      complex(dp):: I3me,Lsm1_2m

      L6m3=  + Lsm1_2m(taugs,taucg,qsq,msq) - I3me(msq,taugs,qsq)*
     & taucs + 2.D0*I3me(msq,taugs,qsq)*taucg**(-1)*taugs*msq -
     & I3me(msq,taugs,qsq)*taucg

      return
      end
