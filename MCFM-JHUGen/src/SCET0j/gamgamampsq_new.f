!---  C.W Jan 16
!---- routine to fill hard term for pp=>Gam Gam (re-written)
      subroutine gamgamampsq_new(order,p,i1,i2,i3,i4,qqb,hard,tlrem)
!======tlrem is two-loop remainder and is the bit which goes like the
!======sum of the quark charges and doesnt talk to the rest of the calc.

!=====note this is now written entirely in terms of |<M|M>|^2 so there are no
!=====longer any complex variables.

!===== this one has Nigel et als formula for I2 as is with no (-mu^2/s)^2e added
!----- in H2 prefactor
      implicit none 
      include 'types.f' 
      include 'constants.f' 
      include 'nf.f' 
      include 'mxpart.f' 
      include 'sprods_com.f' 
      include 'scale.f' 
      include 'ewcouple.f' 
      include 'qcdcouple.f'
      include 'scet_const.f'
      real(dp):: p(mxpart,4)
      integer:: order,i1,i2,i3,i4
      real(dp):: qqb,hard(2),tlrem,Qpie
      real(dp):: ss,tt,uu
      real(dp):: BigS,BigT,BigU,BigX,BigY,xx,yy
!      complex(dp):: lnrat,coeff(2),onelsq
!      complex(dp):: coeff1fin,coeff1C0,Eq13curlybr,twoloopcorr
!     complex(dp):: m2gagafin,renfac_gaga2loop,renfac_gaga2loop_2
      real(dp) Fin1x1,Fin2x0
      
      hard(:)=zip
      call dotem(4,p,s)

!=====Mandelstam invariants 
      ss=s(i1,i2)
      tt=s(i1,i3)
      uu=s(i1,i4)
      
!----- Big X and Big Y defined for physical invariants, so real(dp) is fine 
      BigX=log(-tt/ss)
      BigY=log(-uu/ss)
      BigS=log(ss/musq) 
      BigU=log(-uu/musq)
      BigT=log(-tt/musq)
      xx=-tt/ss
      yy=-uu/ss 
      tlrem=zip

!======LO (minus an overall esq**2*xn factor)
      qqb=8._dp*(tt/uu+uu/tt)

!===== this is 2Re(<0|1>) matched 
      hard(1)= (4*CF*(-42*ss**2 + 12*BigX**2*ss**2 + 6*BigY**2*ss**2 + 
     -      7*Pi**2*ss**2 - 84*ss*uu + 12*BigX**2*ss*uu + 
     -      14*Pi**2*ss*uu - 84*uu**2 + 6*BigX**2*uu**2 + 
     -      6*BigY**2*uu**2 + 14*Pi**2*uu**2 + 
     -      24*tt*uu*Log(musq/ss) - 
     -      6*(ss**2 + 2*ss*uu + 2*uu**2)*Log(musq/ss)**2 + 
     -      12*ss*uu*Log(-(musq/tt)) - 6*uu**2*Log(-(musq/tt)) - 
     -      6*(3*ss**2 + 4*ss*uu + uu**2)*Log(-(musq/uu))))/
     -     (3.*tt*uu)
!===== extract factor of tree to recycle some DY code later
      hard(1)=hard(1)/qqb

!===== call AGTY functions      
      call AGTYassemble_new(ss,tt,uu,Fin1x1,Fin2x0,tlrem)
      
!======this is 2Re(<0|2>) + <1|1> matched (in units of as/2pi**2
      hard(2)= -(-2764*CA*CF*Pi**2*ss**2 + 3528*CF**2*Pi**2*ss**2 - 
     -     1008*BigX**2*CF**2*Pi**2*ss**2 - 
     -     504*BigY**2*CF**2*Pi**2*ss**2 + 448*CF*NF*Pi**2*ss**2 + 
     -     75*CA*CF*Pi**4*ss**2 - 294*CF**2*Pi**4*ss**2 + 
     -     54*Fin1x1*ss*uu + 54*Fin2x0*ss*uu - 5528*CA*CF*Pi**2*ss*uu + 
     -     7056*CF**2*Pi**2*ss*uu - 1008*BigX**2*CF**2*Pi**2*ss*uu + 
     -     896*CF*NF*Pi**2*ss*uu + 150*CA*CF*Pi**4*ss*uu - 
     -     588*CF**2*Pi**4*ss*uu + 54*Fin1x1*uu**2 + 54*Fin2x0*uu**2 - 
     -     5528*CA*CF*Pi**2*uu**2 + 7056*CF**2*Pi**2*uu**2 - 
     -     504*BigX**2*CF**2*Pi**2*uu**2 - 
     -     504*BigY**2*CF**2*Pi**2*uu**2 + 896*CF*NF*Pi**2*uu**2 + 
     -     150*CA*CF*Pi**4*uu**2 - 588*CF**2*Pi**4*uu**2 + 
     -     396*CA*CF*ss**2*Zeta3 - 72*CF*NF*ss**2*Zeta3 + 
     -     792*CA*CF*ss*uu*Zeta3 - 144*CF*NF*ss*uu*Zeta3 + 
     -     792*CA*CF*uu**2*Zeta3 - 144*CF*NF*uu**2*Zeta3 - 
     -     24*CF*(-72*CF*tt*uu - 11*CA*(ss**2 + 2*ss*uu + 2*uu**2) + 
     -        2*NF*(ss**2 + 2*ss*uu + 2*uu**2))*Log(musq/ss)**3 - 
     -     216*CF**2*(ss**2 + 2*ss*uu + 2*uu**2)*Log(musq/ss)**4 - 
     -     504*CF**2*Pi**2*ss*uu*Log(-(musq/tt)) - 
     -     252*CF**2*Pi**2*(2*ss - uu)*uu*Log(-(musq/tt)) + 
     -     252*CF**2*Pi**2*uu**2*Log(-(musq/tt)) + 
     -     1512*CF**2*Pi**2*ss**2*Log(-(musq/uu)) + 
     -     2016*CF**2*Pi**2*ss*uu*Log(-(musq/uu)) + 
     -     504*CF**2*Pi**2*uu**2*Log(-(musq/uu)) - 
     -     12*CF*Log(musq/ss)**2*
     -      (-233*CA*ss**2 + 90*CF*ss**2 - 72*BigX**2*CF*ss**2 - 
     -        36*BigY**2*CF*ss**2 + 38*NF*ss**2 + 6*CA*Pi**2*ss**2 - 
     -        42*CF*Pi**2*ss**2 - 466*CA*ss*uu + 612*CF*ss*uu - 
     -        72*BigX**2*CF*ss*uu + 76*NF*ss*uu + 12*CA*Pi**2*ss*uu - 
     -        84*CF*Pi**2*ss*uu - 466*CA*uu**2 + 612*CF*uu**2 - 
     -        36*BigX**2*CF*uu**2 - 36*BigY**2*CF*uu**2 + 76*NF*uu**2 + 
     -        12*CA*Pi**2*uu**2 - 84*CF*Pi**2*uu**2 - 
     -        36*CF*ss*uu*Log(-(musq/tt)) + 
     -        18*CF*uu**2*Log(-(musq/tt)) + 
     -        18*CF*uu*(-2*ss + uu)*Log(-(musq/tt)) + 
     -        36*CF*(3*ss**2 + 4*ss*uu + uu**2)*Log(-(musq/uu))) - 
     -     6*CF*Log(musq/ss)*(-804*CA*ss**2 + 1512*CF*ss**2 - 
     -        432*BigX**2*CF*ss**2 - 216*BigY**2*CF*ss**2 + 
     -        120*NF*ss**2 + 157*CA*Pi**2*ss**2 - 22*NF*Pi**2*ss**2 - 
     -        1608*CA*ss*uu + 3024*CF*ss*uu - 432*BigX**2*CF*ss*uu + 
     -        240*NF*ss*uu + 314*CA*Pi**2*ss*uu - 336*CF*Pi**2*ss*uu - 
     -        44*NF*Pi**2*ss*uu - 1608*CA*uu**2 + 3024*CF*uu**2 - 
     -        216*BigX**2*CF*uu**2 - 216*BigY**2*CF*uu**2 + 
     -        240*NF*uu**2 + 314*CA*Pi**2*uu**2 - 336*CF*Pi**2*uu**2 - 
     -        44*NF*Pi**2*uu**2 - 216*CF*ss*uu*Log(-(musq/tt)) - 
     -        108*CF*(2*ss - uu)*uu*Log(-(musq/tt)) + 
     -        108*CF*uu**2*Log(-(musq/tt)) + 
     -        216*CF*(3*ss**2 + 4*ss*uu + uu**2)*Log(-(musq/uu))))/
     -  (54.*tt*uu)
!=======factor out Born 
      hard(2)=hard(2)/qqb

      return
      end
