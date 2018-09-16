      SUBROUTINE TAUTAUWW(I1,I2,RPLPL,RPLMN,RMNPL,RMNMN,FAC)
* THE MATRIX ELEMENT FOR THE PROCESS
*    QBAR(I1) Q(I2) ---> TAU TAUBAR
* FOLLOWED BY
*    TAUBAR ---> NUBAR(3) NU(4) E^+(5)
*    TAU    ---> NU(6)    E^-(7)  NUBAR(8)
*
      implicit none
      include 'types.f'
      include 'cplx.h'
c      IMPLICIT real(dp):: (A-H,O-Z)
      complex(dp):: RPLPL,RPLMN,RMNPL,RMNMN,DEN,SPROD
      complex(dp):: SPL(10,10),SMN(10,10)
      real(dp):: FAC,DOTKS,CT,CTBAR,D12,D34,D35,D45,D67,D68,D78,
     & XT,XTBAR,RMW2,RMGW,RMT2,RMB2,RMGT,SCALE
      real(dp):: RMT,RGT,RMW,RGW,RMB,RMTLO,RMTUP,GW,GS,PLAB(4,10)
      integer:: I1,I2,K
      COMMON/COUPS/GW,GS
      COMMON/CSTD/SPL,SMN
      COMMON/MOM/PLAB
      COMMON/PARS/RMT,RGT,RMW,RGW,RMB,RMTLO,RMTUP
      integer, save::init=0
      SAVE RMW2,RMGW,RMT2,RMB2,RMGT,SCALE
!$omp threadprivate(RMW2,RMGW,RMT2,RMB2,RMGT,SCALE,INIT)
!$omp threadprivate(/COUPS/,/PARS/,/MOM/,/CSTD/)  


      IF(INIT.NE.0) GOTO 1
      INIT=1
      SCALE=1D1
      RMW2=RMW**2/SCALE
      RMGW=RMW*RGW/SCALE
      RMT2=RMT**2/SCALE
      RMB2=RMB**2/SCALE
      RMGT=RMT*RGT/SCALE
  1   CONTINUE
*
* THE DENOMINATOR
 
      D12=DOTKS(1,2)/SCALE
      D34=DOTKS(3,4)/SCALE
      D35=DOTKS(3,5)/SCALE
      D45=DOTKS(4,5)/SCALE
      D67=DOTKS(6,7)/SCALE
      D68=DOTKS(6,8)/SCALE
      D78=DOTKS(7,8)/SCALE
      XT=   2._dp*(D67+D68+D78) + RMB2
      XTBAR=2._dp*(D34+D35+D45) + RMB2
      DEN=2._dp*D12*
     & cplx2(XT   -RMT2,RMGT)*
     & cplx2(XTBAR-RMT2,RMGT)*
     & cplx2(2._dp*D45-RMW2,RMGW)*
     & cplx2(2._dp*D78-RMW2,RMGW)

*
* THE AUXILIARY VECTORS FROM THE T AND TBAR MOMENTA
      CT   =XT   /(2._dp*(D68+D78))
      CTBAR=XTBAR/(2._dp*(D34+D45))
      DO 11 K=1,4
      PLAB(K,9) =PLAB(K,3)+(1._dp-CTBAR)*PLAB(K,4)+PLAB(K,5)
      PLAB(K,10)=PLAB(K,6)+PLAB(K,7)+(1._dp-CT   )*PLAB(K,8)
   11 CONTINUE
*
      CALL STD
* THE SPINOR PART OF THE NUMERATOR
      SPROD=-SPL(8,10)*SMN(9,4)/cplx1(SCALE)
      RPLMN=SPROD*SMN(10,I2)*SPL(I1,9)/SCALE
      RPLPL=RMT2*SPL(8,I1)*SMN(I2,4)/SCALE
      RMNMN=SPROD*SMN(10,I1)*SPL(I2,9)/SCALE
      RMNPL=RMT2*SPL(8,I2)*SMN(I1,4)/SCALE

      FAC=GW**8*(1._dp/SCALE)**4*
     & D67*D35*128._dp**2/ABS(DEN)**2

* ADD NU,NUBAR SPIN SUM, Q,QBAR AVERAGE, COLOUR SUM &AVERAGE
      FAC=FAC *4._dp         /4._dp           *3._dp       /9._dp

      RETURN
      END
