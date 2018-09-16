      subroutine gagaOneloopsq(ss,tt,uu,onelsq)
!==== returns the MS-SCET renormalized one-loop sq amplitude for gaga
!==== C. Williams August 2015
      
!==== return answer in units of N_C similar to other routines 
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'scale.f' 
      real(dp), intent(in) :: ss,tt,uu
      real(dp), intent(out) ::onelsq
      real(dp)  xtag_ipi
      common/xtag_ipi/xtag_ipi
!$omp threadprivate(/xtag_ipi/)

      
!======AGTY function for one-loop squared finite pieces
      real(dp) :: AGTYG1s,BigX,BigY

!=====finite bits squared (from AGTY)

      BigX=log(-tt/ss)
      BigY=log(-uu/ss)

      
      onelsq=Cf**2*AGTYG1s(tt,uu,Bigx,Bigy)

!==== xtag_ipi is sign of i(pi) in Log[mu^2/-s]-> Log[mu^2/s]+xtag_ipi*(im*pi)
      
      
!========first renormalization term (linear in (I(eps)+Z(eps))
      onelsq=onelsq+
     &  (2*CF**2*(-6*(3*tt**2 + 4*tt*uu + 3*uu**2)*
     -       Log(musq/ss)**3 - 
     -      6*Log(musq/ss)**2*
     -       (2*tt**2 + 2*BigX**2*tt**2 + 
     -         BigY**2*tt**2 + 12*tt*uu + 
     -         2*BigX**2*tt*uu + 2*BigY**2*tt*uu + 
     -         2*uu**2 + BigX**2*uu**2 + 
     -         2*BigY**2*uu**2 - 
     -         uu*(2*tt + 3*uu)*Log(-(musq/tt)) - 
     -         tt*(3*tt + 2*uu)*Log(-(musq/uu))) + 
     -      Log(musq/ss)*
     -       (126*tt**2 - 36*BigX**2*tt**2 - 
     -         18*BigY**2*tt**2 + 3*Pisq*tt**2 - 
     -         36*BigX**2*tt*uu - 36*BigY**2*tt*uu + 
     -         4*Pisq*tt*uu + 126*uu**2 - 
     -         18*BigX**2*uu**2 - 36*BigY**2*uu**2 + 
     -         3*Pisq*uu**2 + 
     -         48*BigX*Pisq*tt**2*xtag_ipi + 
     -         24*BigY*Pisq*tt**2*xtag_ipi + 
     -         48*BigX*Pisq*tt*uu*xtag_ipi + 
     -         48*BigY*Pisq*tt*uu*xtag_ipi + 
     -         24*BigX*Pisq*uu**2*xtag_ipi + 
     -         48*BigY*Pisq*uu**2*xtag_ipi + 
     -         54*Pisq*tt**2*xtag_ipi**2 + 
     -         72*Pisq*tt*uu*xtag_ipi**2 + 
     -         54*Pisq*uu**2*xtag_ipi**2 + 
     -         18*uu*(2*tt + 3*uu)*Log(-(musq/tt)) + 
     -         18*tt*(3*tt + 2*uu)*Log(-(musq/uu)))
     -       + Pisq*(-7*tt**2 + 2*BigX**2*tt**2 + 
     -         BigY**2*tt**2 + 2*BigX**2*tt*uu + 
     -         2*BigY**2*tt*uu - 7*uu**2 + 
     -         BigX**2*uu**2 + 2*BigY**2*uu**2 + 
     -         72*BigX*tt**2*xtag_ipi + 
     -         36*BigY*tt**2*xtag_ipi + 
     -         72*BigX*tt*uu*xtag_ipi + 
     -         72*BigY*tt*uu*xtag_ipi + 
     -         36*BigX*uu**2*xtag_ipi + 
     -         72*BigY*uu**2*xtag_ipi + 
     -         12*tt**2*xtag_ipi**2 + 
     -         12*BigX**2*tt**2*xtag_ipi**2 + 
     -         6*BigY**2*tt**2*xtag_ipi**2 + 
     -         72*tt*uu*xtag_ipi**2 + 
     -         12*BigX**2*tt*uu*xtag_ipi**2 + 
     -         12*BigY**2*tt*uu*xtag_ipi**2 + 
     -         12*uu**2*xtag_ipi**2 + 
     -         6*BigX**2*uu**2*xtag_ipi**2 + 
     -         12*BigY**2*uu**2*xtag_ipi**2 - 
     -         uu*(2*tt + 3*uu)*(1 + 6*xtag_ipi**2)*
     -          Log(-(musq/tt)) - 
     -         tt*(3*tt + 2*uu)*(1 + 6*xtag_ipi**2)*
     -          Log(-(musq/uu)))))/(3.*tt*uu)
!=====second renormalization bit ( (I(eps)+Z(eps)*conj(I(eps+Z(eps)))

      onelsq=onelsq+
     & (CF**2*(tt/uu + uu/tt)*
     -    (324*Pi**2*xtag_ipi**2 + Pi**4*(1 + 6*xtag_ipi**2)**2 - 
     -      36*Pi**2*(1 - 6*xtag_ipi**2)*Log(musq/ss) - 
     -      12*(-27 + Pi**2*(1 - 6*xtag_ipi**2))*Log(musq/ss)**2 + 
     -      216*Log(musq/ss)**3 + 36*Log(musq/ss)**4))/18._dp
      return
      end
