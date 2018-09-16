      subroutine gamgamampsq(order,p,i1,i2,i3,i4,qqb,hard,tlrem) 
!+==== routine to fill LO, and hard M.E's for pp=>Gam gam 
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
      real(dp):: BigS,BigT,BigU,BigX,BigY,xx,yy,onelsq
      complex(dp):: lnrat,coeff(2)
      complex(dp):: coeff1fin,coeff1C0,Eq13curlybr,twoloopcorr
      real(dp)  xtag_ipi
      common/xtag_ipi/xtag_ipi
      complex(dp):: m2gagafin,renfac_gaga2loop,renfac_gaga2loop_2
      include 'cplx.h'
!$omp threadprivate(/xtag_ipi/)

      hard(:)=zip
      coeff(:)=czip
      call dotem(4,p,s)

!==== xtag_ipi is sign of i(pi) in Log[mu^2/-s]-> Log[mu^2/s]+xtag_ipi*(im*pi)
      xtag_ipi=1._dp

!=====Mandelstam invariants 
      ss=s(i1,i2)
      tt=s(i1,i3)
      uu=s(i1,i4)
          
!----- Big X and Big Y defined for physical invariants, so this is ok  
      BigX=log(-tt/ss)
      BigY=log(-uu/ss)
      BigS=log(ss/musq) 
      BigU=log(-uu/musq)
      BigT=log(-tt/musq)
      xx=-tt/ss
      yy=-uu/ss 
      tlrem=zip

!===== LO (minus overall esq**2 * xn) 

!--- Tree-level from Eq.(2.37)
!--- overall factor of N from Eq.(2.35) applied later
      qqb=8._dp*(tt/uu+uu/tt)     

!--- renormalized interference of 1-loop and tree-level amplitude
!--- <M0|M1> in notation of Eq.(2.38)
!--- overall factor of 2*N*CF applied later
!--- 2 from Eq.(2.38) and N*CF from Eq.(3.8)
      coeff(1)=   -(  -48 - (12*ss)/tt - (12*ss)/uu + 
     -  (12*ss**2)/(tt*uu) - (4*BigX**2*ss**2)/(tt*uu) - 
     -  (4*BigY**2*ss**2)/(tt*uu) + 
     -  (8*im*Pi*ss**2)/(tt*uu) - 
     -  (8*BigX*im*Pi*ss**2)/(tt*uu) - 
     -  (8*BigY*im*Pi*ss**2)/(tt*uu) + (4*tt)/uu - 
     -  (4*BigX**2*tt)/uu - (8*im*Pi*tt)/uu - 
     -  (8*BigX*im*Pi*tt)/uu - (14*Pi**2*tt)/(3.*uu) + 
     -  (4*uu)/tt - (4*BigY**2*uu)/tt - (8*im*Pi*uu)/tt - 
     -  (8*BigY*im*Pi*uu)/tt - (14*Pi**2*uu)/(3.*tt) - 
     -  (8*ss**2*Log(musq/ss))/(tt*uu) + 
     -  (8*tt*Log(musq/ss))/uu - 
     -  (8*im*Pi*tt*Log(musq/ss))/uu + 
     -  (8*uu*Log(musq/ss))/tt - 
     -  (8*im*Pi*uu*Log(musq/ss))/tt + 
     -  (4*tt*Log(musq/ss)**2)/uu + 
     -  (4*uu*Log(musq/ss)**2)/tt - 4*Log(-(musq/tt)) - 
     -  (12*ss*Log(-(musq/tt)))/tt - 4*Log(-(musq/uu)) - 
     -  (12*ss*Log(-(musq/uu)))/uu)

c      write(6,*) 'old coeff(1)',coeff(1)

      coeff1fin =
     &  + im*pi * ( 16.D0 + 12.D0*xx**(-1)*yy + 12.D0*xx*yy**(-1) + 16.D
     &    0*BigY*xx**(-1)*yy + 16.D0*BigY + 8.D0*BigY*xx*yy**(-1) + 8.D0
     &    *BigX*xx**(-1)*yy + 16.D0*BigX + 16.D0*BigX*xx*yy**(-1) )
      coeff1fin = coeff1fin - 28.D0*xx**(-1)*yy - 28.D0*xx*yy**(-1) + 8.
     &    D0*BigU + 12.D0*BigU*xx*yy**(-1) + 12.D0*BigT*xx**(-1)*yy + 8.
     &    D0*BigT - 12.D0*BigS*xx**(-1)*yy - 16.D0*BigS - 12.D0*BigS*xx
     &    *yy**(-1) + 8.D0*BigY**2*xx**(-1)*yy + 8.D0*BigY**2 + 4.D0*
     &    BigY**2*xx*yy**(-1) + 4.D0*BigX**2*xx**(-1)*yy + 8.D0*BigX**2
     &     + 8.D0*BigX**2*xx*yy**(-1)

      coeff1C0 =
     &  + im*pi * (  - 3.D0/2.D0 + BigS )
     &  + 3.D0/2.D0*BigS - 1.D0/2.D0*BigS**2 + 7.D0/12.D0*pisq

! JC result is commented out below
!      coeff(1)=coeff1fin+coeff1C0*qqb

c      write(6,*) 'new coeff(1)',coeff(1)

      Eq13curlybr =
     &  + im*pi * ( 5.D0/3.D0*CF*TR*nf - 19.D0/9.D0*CF*TR*nf*BigS + 1.D0
     &    /3.D0*CF*TR*nf*BigS**2 - 1.D0/12.D0*CF*TR*nf*pisq - 67.D0/12.D
     &    0*CF*CA + 233.D0/36.D0*CF*CA*BigS - 11.D0/12.D0*CF*CA*BigS**2
     &     + 23.D0/48.D0*CF*CA*pisq - 1.D0/6.D0*CF*CA*pisq*BigS - 9.D0/
     &    4.D0*CF**2*BigS + 9.D0/4.D0*CF**2*BigS**2 - 1.D0/2.D0*CF**2*
     &    BigS**3 - 7.D0/8.D0*CF**2*pisq + 7.D0/12.D0*CF**2*pisq*BigS )
      Eq13curlybr = Eq13curlybr - 5.D0/3.D0*CF*TR*nf*BigS + 19.D0/18.D0
     &    *CF*TR*nf*BigS**2 - 1.D0/9.D0*CF*TR*nf*BigS**3 + 1.D0/6.D0*CF
     &    *TR*nf*zeta3 - 28.D0/27.D0*CF*TR*nf*pisq + 11.D0/36.D0*CF*TR*
     &    nf*pisq*BigS + 67.D0/12.D0*CF*CA*BigS - 233.D0/72.D0*CF*CA*
     &    BigS**2 + 11.D0/36.D0*CF*CA*BigS**3 - 11.D0/24.D0*CF*CA*zeta3
     &     + 691.D0/216.D0*CF*CA*pisq - 157.D0/144.D0*CF*CA*pisq*BigS
     &     + 1.D0/12.D0*CF*CA*pisq*BigS**2 - 25.D0/288.D0*CF*CA*pisq**2
     &     + 9.D0/8.D0*CF**2*BigS**2 - 3.D0/4.D0*CF**2*BigS**3 + 1.D0/8.
     &    D0*CF**2*BigS**4 - 9.D0/8.D0*CF**2*pisq + 19.D0/8.D0*CF**2*
     &    pisq*BigS - 19.D0/24.D0*CF**2*pisq*BigS**2 + 49.D0/288.D0*
     &    CF**2*pisq**2

! JC result is commented out below
      twoloopcorr=CF**2*coeff1C0*coeff1fin+Eq13curlybr*qqb
c      write(6,*) 'twoloopcorr',twoloopcorr
c      write(6,*) 'ciaran',real(renfac_gaga2loop(ss,tt,uu,BigX,BigY),kind=dp)
c      write(6,*)
c      pause
      
!==== call two loop 2*Re(|<0|2>|) from AGTY
!==== note tlrem is the piece which does not go like Q(j)^4 
!==== and is thus treated as a separate finite correction
      call gagaAGTYassemble(ss,tt,uu,tlrem,Qpie) 

      coeff(2)=cplx1(Qpie)
      
c      write(6,*) 'ss,tt,uu,musq',ss,tt,uu,musq
      
!==== next we need the one-loop squared piece,
c      call gagaOneloopsq(ss,tt,uu,onelsq)       
c      write(6,*) 'old onelsq',onelsq

! JC result is commented out below
      call gagaAGTYassemble1loop(ss,tt,uu,coeff1fin,coeff1C0,qqb,onelsq)
      
c      write(6,*) 'new onelsq',onelsq
c      write(6,*)
c      pause
     
!==== overall factor =
!====                 x cf (color factor of 3.8 in AGT-Y (with no xn)
!====                 x 2_dp because we take 2*Re(intf) 
!====
!+===== overall factor is then ason2pi 
!===== hard function dressed with couplings 

      hard(1)=2._dp*cf*real(coeff(1)/qqb,kind=dp)

!==== overall factor = 1 for AGTY piece, 2 for ren fac (its written for <0|2>)  
!===== and one for the one-loop squared piece
!+===== overall factor is then ason2pi**2 

c      hard(2)=real(coeff(2)/qqb,kind=dp)
c     &    +2._dp*real(renfac_gaga2loop(ss,tt,uu,BigX,BigY)/qqb,kind=dp)
c     &       +real(onelsq/qqb,kind=dp)
  
      hard(2)=real(coeff(2)/qqb,kind=dp)
     &    +2._dp*real(twoloopcorr/qqb,kind=dp)
     &       +real(onelsq/qqb,kind=dp)
  
c      write(6,*) 'qqb        ',qqb
c      write(6,*) 'hard(1)*qqb',hard(1)*qqb
c      write(6,*) 'hard(2)*qqb',hard(2)*qqb
c      write(6,*) 'tlrem  ',tlrem
c      write(6,*) 'Qpie   ',Qpie
c      write(6,*) 'onelsq ',onelsq
  
      return 
      end


       function renfac_gaga2loop(ss,tt,uu,BigX,BigY) 
      implicit none 
!====== MS-bar renormalziation for ga ga process 
!====== following is the piece which is added to the finite 
!===== function of AGT-Y 
!===== in units of alpha_S/2pi and this is <M_0 | M _2 >  
      
      include 'types.f' 
      complex(dp):: renfac_gaga2loop
      include 'constants.f' 
      include 'nf.f' 
      include 'scale.f' 
      include 'scet_const.f'
      real(dp):: ss,tt,uu,BigX,BigY
      real(dp)  xtag_ipi
      common/xtag_ipi/xtag_ipi
!$omp threadprivate(/xtag_ipi/)
!==== xtag_ipi is sign of i(pi) in Log[mu^2/-s]-> Log[mu^2/s]+xtag_ipi*(im*pi)

      renfac_gaga2loop=
     &   (CF*((tt/uu + uu/tt)*
     -       (Pi**2*(-32*CA + 8*NF + 3*(-CA + CF)*Pi**2) + 
     -         36*(-11*CA + 2*NF)*zeta3 + 
     -         6*(im*Pi*xtag_ipi + Log(musq/ss))*
     -          (-18*CF*Pi**2 + 2*NF*(60 + Pi**2) + 
     -            CA*(-804 + 25*Pi**2) + 
     -            2*(im*Pi*xtag_ipi + Log(musq/ss))*
     -             (-233*CA + 81*CF + 38*NF - 3*(-2*CA + CF)*Pi**2 + 
     -               (im*Pi*xtag_ipi + Log(musq/ss))*
     -                (-22*CA + 54*CF + 4*NF + 
     -                  9*CF*(im*Pi*xtag_ipi + Log(musq/ss)))))) + 
     -      (36*CF*(Pi**2 - 18*(im*Pi*xtag_ipi + Log(musq/ss)) - 
     -           6*(im*Pi*xtag_ipi + Log(musq/ss))**2)*
     -         (-7*tt**2 + 2*BigX**2*tt**2 + BigY**2*tt**2 + 
     -           4*BigX*im*Pi*tt**2 + 2*BigY*im*Pi*tt**2 + 
     -           2*BigX**2*tt*uu + 2*BigY**2*tt*uu + 
     -           4*BigX*im*Pi*tt*uu + 4*BigY*im*Pi*tt*uu - 7*uu**2 + 
     -           BigX**2*uu**2 + 2*BigY**2*uu**2 + 
     -           2*BigX*im*Pi*uu**2 + 4*BigY*im*Pi*uu**2 + 
     -           (3*tt**2 + 4*tt*uu + 3*uu**2)*
     -            (im*Pi*xtag_ipi + Log(musq/ss)) - 
     -           uu*(2*tt + 3*uu)*Log(-(musq/tt)) - 
     -           3*tt**2*Log(-(musq/uu)) - 2*tt*uu*Log(-(musq/uu))))/
     -       (tt*uu)))/108. 
      
      return 
      end

!===== following is old code for checking
c$$$
c$$$      function renfac_gaga2loop_2(ss,tt,uu,BigX,BigY) 
c$$$      implicit none 
c$$$!====== MS-bar renormalziation for ga ga process 
c$$$!====== following is the piece which is added to the finite 
c$$$!===== function of AGT-Y 
c$$$!===== in units of alpha_S/2pi and this is <M_0 | M _2 >  
c$$$      
c$$$      include 'types.f' 
c$$$      complex(dp):: renfac_gaga2loop_2
c$$$      include 'constants.f' 
c$$$      include 'nf.f' 
c$$$      include 'scale.f' 
c$$$      include 'scet_const.f'
c$$$      real(dp):: ss,tt,uu,BigX,BigY
c$$$
c$$$      renfac_gaga2loop_2=
c$$$     & (CF*(36*CF*(Pi**2 + 
c$$$     -         18*(im*Pi - Log(musq/ss)) - 
c$$$     -         6*(-(im*Pi) + Log(musq/ss))**2)*
c$$$     -       (-7*tt**2 + 2*BigX**2*tt**2 + 
c$$$     -         BigY**2*tt**2 + 4*BigX*im*Pi*tt**2 + 
c$$$     -         2*BigY*im*Pi*tt**2 + 2*BigX**2*tt*uu + 
c$$$     -         2*BigY**2*tt*uu + 4*BigX*im*Pi*tt*uu + 
c$$$     -         4*BigY*im*Pi*tt*uu - 7*uu**2 + 
c$$$     -         BigX**2*uu**2 + 2*BigY**2*uu**2 + 
c$$$     -         2*BigX*im*Pi*uu**2 + 
c$$$     -         4*BigY*im*Pi*uu**2 + 
c$$$     -         (3*tt**2 + 4*tt*uu + 3*uu**2)*
c$$$     -          (-(im*Pi) + Log(musq/ss)) - 
c$$$     -         uu*(2*tt + 3*uu)*Log(-(musq/tt)) - 
c$$$     -         3*tt**2*Log(-(musq/uu)) - 
c$$$     -         2*tt*uu*Log(-(musq/uu))) + 
c$$$     -      (tt**2 + uu**2)*
c$$$     -       (Pi**2*(-32*CA + 8*NF + 
c$$$     -            3*(-CA + CF)*Pi**2) + 
c$$$     -         6*(-(im*Pi) + Log(musq/ss))*
c$$$     -          (-18*CF*Pi**2 + 2*NF*(60 + Pi**2) + 
c$$$     -            CA*(-804 + 25*Pi**2) + 
c$$$     -            2*(-(im*Pi) + Log(musq/ss))*
c$$$     -             (-233*CA + 81*CF + 38*NF - 
c$$$     -               3*(-2*CA + CF)*Pi**2 + 
c$$$     -               (-(im*Pi) + Log(musq/ss))*
c$$$     -                (-22*CA + 54*CF + 4*NF - 
c$$$     -                  9*CF*im*Pi + 9*CF*Log(musq/ss))
c$$$     -               )) + 36*(-11*CA + 2*NF)*Zeta3)))
c$$$     -   /(108.*tt*uu)
c$$$      
c$$$      return 
c$$$      end


     
