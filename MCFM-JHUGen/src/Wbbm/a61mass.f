      subroutine a61mass(k1,k2,k3,k4,k5,k6,mqsq,a61mm,a61mp,a61pm,a61pp,
     &                   a6treemm,a6treemp,a6treepm,a6treepp)
      implicit none
      include 'types.f'
c--- Routine to compute the 1-loop primitive amplitude A6;1 for the Wbb
c--- process, where the mass of the b-quark is kept non-zero.
c--- The labels on this routine refer to the momenta in "mom" (passed
c--- via common block), in which the massive momenta have been made massless
c--- a la Rodrigo and appear in positions 5 and 6.
c--- The tree-level amplitudes are also returned
c---
c---     0 -> q(k1) + qb(k4) + W(->e(k6)+nubar(k5)) + Q(k3) + Qbar(k2)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'masses.f'
      include 'momwbbm.f'
      include 'nflav.f'
      include 'scale.f'
      include 'Wbbmlabels.f'
      include 'zprods_com.f'
      include 'swapxz.f'
      integer:: k1,k2,k3,k4,k5,k6,iep,nu
      logical:: checkcoeffs,numcheck,writescalars
      real(dp):: mqsq,ren,p2(4),p3(4),p1Dp2,p3Dp4,p1Dp3,p2Dp4,
     & pole(-2:-1),mhloopsq,mq
      complex(dp):: a61mm,a61mp,a61pm,a61pp
      complex(dp):: a6treemm,a6treemp,a6treepm,a6treepp
      complex(dp):: a61mmep(-2:0),a61mpep(-2:0),
     & a61pmep(-2:0),a61ppep(-2:0)
      complex(dp):: ampALC(2,2),ampBLC(2,2)
      complex(dp):: ALCmm(-2:0),ALCmp(-2:0),ALCpm(-2:0),ALCpp(-2:0)
      complex(dp):: BLCmm(-2:0),BLCmp(-2:0),BLCpm(-2:0),BLCpp(-2:0)
      complex(dp):: ASLmp(-2:0),ASLmm(-2:0),ASLpm(-2:0),ASLpp(-2:0)
      complex(dp):: Afmm(-2:0),Afmp(-2:0),Afpm(-2:0),Afpp(-2:0)
      complex(dp):: Ahmm(-2:0),Ahmp(-2:0),Ahpm(-2:0),Ahpp(-2:0)
      complex(dp):: coeff2,coeff1,atree,atreesl,
     & Lnrat,xl12,xl34,xl13,xl24
      parameter (checkcoeffs=.false.)
      common/numcheck/numcheck
      common/writescalars/writescalars
!$omp threadprivate(/numcheck/,/writescalars/)

      writescalars=.false.
      swapxz=.true.
      
      mq=sqrt(mqsq)
      
c--- first calculate the tree-level amplitudes
      call a6treemass(k1,k2,k3,k4,k5,k6,mq,
     & a6treemm,a6treemp,a6treepm,a6treepp)

c--- compute logarithms appearing in pole contributions      
      do nu=1,4
      p2(nu)=bp*mom(k2,nu)+bm*mom(k3,nu)
      p3(nu)=bp*mom(k3,nu)+bm*mom(k2,nu)
      enddo
      p1Dp2=mom(k1,4)*p2(4)-mom(k1,1)*p2(1)
     &     -mom(k1,2)*p2(2)-mom(k1,3)*p2(3)
      p1Dp3=mom(k1,4)*p3(4)-mom(k1,1)*p3(1)
     &     -mom(k1,2)*p3(2)-mom(k1,3)*p3(3)
      p2Dp4=mom(k4,4)*p2(4)-mom(k4,1)*p2(1)
     &     -mom(k4,2)*p2(2)-mom(k4,3)*p2(3)
      p3Dp4=mom(k4,4)*p3(4)-mom(k4,1)*p3(1)
     &     -mom(k4,2)*p3(2)-mom(k4,3)*p3(3)
      xl12=lnrat(-two*p1Dp2,musq)
      xl13=lnrat(-two*p1Dp3,musq)
      xl24=lnrat(-two*p2Dp4,musq)
      xl34=lnrat(-two*p3Dp4,musq)
      
c--- compute the necessary scalar integrals for leading colour A
      call computescalars(k1,k2,k3,k4,k5,k6,scints)
c      call computescalarsnew(k1,k2,k3,k4,k5,k6,scints)
      
c--- Leading colour A
c--- calculate (-,+) and (+,-) amplitudes separately
      call ALC_mp(k1,k2,k3,k4,k5,k6,coeff,ampALC)
      call sumamp(coeff,scints,ALCmp,'ALC mp')
      if (checkcoeffs) then
        write(6,*) 'ALC: MINUS-PLUS'
        call writecoeffs(k1,k2,k3,k4,k5,k6,coeff)
      endif
      
      call ALC_pm(k1,k2,k3,k4,k5,k6,coeff,ampALC)
      call sumamp(coeff,scints,ALCpm,'ALC pm')
      if (checkcoeffs) then
        write(6,*) 'ALC: PLUS-MINUS'
        call writecoeffs(k1,k2,k3,k4,k5,k6,coeff)
      endif

c--- calculate (-,-) and (+,+) amplitudes, the latter obtained by symmetry
      call ALC_mm(k1,k2,k3,k4,k5,k6,coeff,ampALC)
      call sumamp(coeff,scints,ALCmm,'ALC mm')
      if (checkcoeffs) then
        write(6,*) 'ALC: MINUS-MINUS'
        call writecoeffs(k1,k2,k3,k4,k5,k6,coeff)
      endif

      call ALC_pp(k1,k2,k3,k4,k5,k6,coeff,ampALC)
      call sumamp(coeff,scints,ALCpp,'ALC pp')
      if (checkcoeffs) then
        write(6,*) 'ALC: PLUS-PLUS'
        call writecoeffs(k1,k2,k3,k4,k5,k6,coeff)
      endif

c--- explicit form for poles in ALC
      pole(-2)=-1._dp
      pole(-1)=xl12+xl34+8._dp/3._dp-log(mqsq/musq)
      do iep=-2,-1
      ALCmp(iep)=pole(iep)*a6treemp
      ALCpm(iep)=pole(iep)*a6treepm
      ALCmm(iep)=pole(iep)*a6treemm
      ALCpp(iep)=pole(iep)*a6treepp
      enddo
      
      if (numcheck) then
        call catani(mqsq,k1,k2,k3,k4,k5,k6,coeff2,coeff1)
        write(6,*) 'Catani',coeff2
        write(6,*) 'Catani',coeff1
        write(6,*)
        write(6,*) 'ALCmp(-2)/tree/Catani',ALCmp(-2)/a6treemp/coeff2
        write(6,*) 'ALCmp(-1)/tree/Catani',ALCmp(-1)/a6treemp/coeff1
        write(6,*)
        write(6,*) 'ALCpm(-2)/tree/Catani',ALCpm(-2)/a6treepm/coeff2
        write(6,*) 'ALCpm(-1)/tree/Catani',ALCpm(-1)/a6treepm/coeff1
        write(6,*)
        write(6,*) 'ALCmm(-2)/tree/Catani',ALCmm(-2)/a6treemm/coeff2
        write(6,*) 'ALCmm(-1)/tree/Catani',ALCmm(-1)/a6treemm/coeff1
        write(6,*)
        write(6,*) 'ALCpp(-2)/tree/Catani',ALCpp(-2)/a6treepp/coeff2
        write(6,*) 'ALCpp(-1)/tree/Catani',ALCpp(-1)/a6treepp/coeff1
        write(6,*)
      endif
      
c      return
      
c--- compute the necessary scalar integrals for leading colour B
      call computescalars(k1,k3,k2,k4,k5,k6,scints)
     
c--- Leading colour B
c--- calculate (-,+) and (+,-) amplitudes separately
      call BLC_mp(k1,k2,k3,k4,k5,k6,coeff,ampBLC)
      call sumamp(coeff,scints,BLCmp,'BLC mp')
      if (checkcoeffs) then
        write(6,*) 'BLC: MINUS-PLUS'
        call writecoeffs(k1,k3,k2,k4,k5,k6,coeff)
      endif

      call BLC_pm(k1,k2,k3,k4,k5,k6,coeff,ampBLC)
      call sumamp(coeff,scints,BLCpm,'BLC pm')
      if (checkcoeffs) then
        write(6,*) 'BLC: PLUS-MINUS'
        call writecoeffs(k1,k3,k2,k4,k5,k6,coeff)
      endif

c--- calculate (-,-) and (+,+) amplitudes, the latter obtained by symmetry
      call BLC_mm(k1,k2,k3,k4,k5,k6,coeff,ampBLC)
      call sumamp(coeff,scints,BLCmm,'BLC mm')
      if (checkcoeffs) then
        write(6,*) 'BLC: MINUS-MINUS'
        call writecoeffs(k1,k3,k2,k4,k5,k6,coeff)
      endif

      call BLC_pp(k1,k2,k3,k4,k5,k6,coeff,ampBLC)
      call sumamp(coeff,scints,BLCpp,'BLC pp')
      if (checkcoeffs) then
        write(6,*) 'BLC: PLUS-PLUS'
        call writecoeffs(k1,k3,k2,k4,k5,k6,coeff)
      endif

c--- explicit form for poles in BLC
      pole(-2)=-1._dp
      pole(-1)=xl13+xl24+8._dp/3._dp-log(mqsq/musq)
      do iep=-2,-1
      BLCmp(iep)=-pole(iep)*a6treemp
      BLCpm(iep)=-pole(iep)*a6treepm
      BLCmm(iep)=-pole(iep)*a6treemm
      BLCpp(iep)=-pole(iep)*a6treepp
      enddo
      
      if (numcheck) then
        call catani(mqsq,k1,k3,k2,k4,k5,k6,coeff2,coeff1)
        write(6,*) 'Catani',coeff2
        write(6,*) 'Catani',coeff1
        write(6,*)
        write(6,*) 'BLCmp(-2)/tree/Catani',BLCmp(-2)/a6treemp/coeff2
        write(6,*) 'BLCmp(-1)/tree/Catani',BLCmp(-1)/a6treemp/coeff1
        write(6,*)
        write(6,*) 'BLCpm(-2)/tree/Catani',BLCpm(-2)/a6treepm/coeff2
        write(6,*) 'BLCpm(-1)/tree/Catani',BLCpm(-1)/a6treepm/coeff1
        write(6,*)
        write(6,*) 'BLCmm(-2)/tree/Catani',BLCmm(-2)/a6treemm/coeff2
        write(6,*) 'BLCmm(-1)/tree/Catani',BLCmm(-1)/a6treemm/coeff1
        write(6,*)
        write(6,*) 'BLCpp(-2)/tree/Catani',BLCpp(-2)/a6treepp/coeff2
        write(6,*) 'BLCpp(-1)/tree/Catani',BLCpp(-1)/a6treepp/coeff1
        write(6,*)
      endif
      
c--- subleading colour Asl
c--- calculate (-,+) and (+,-) amplitudes separately

      call BDKfillmm(k1,k2,k3,k4,k5,k6,mq,za,zb,ASLmm,ASLpp)
      call BDKfillmp(k1,k2,k3,k4,k5,k6,mq,za,zb,ASLmp)
      call BDKfillmp(k1,k3,k2,k4,k5,k6,mq,za,zb,ASLpm)

      if (numcheck) then
c        write(6,89) 'ASL: MINUS-PLUS'
        write(6,89) 'ASL mp(-2) =',ASLmp(-2)
        write(6,89) 'ASL mp(-1) =',ASLmp(-1)
        write(6,89) 'ASL mp( 0) =',ASLmp( 0)
        write(6,89) 
c        write(6,89) 'ASL: PLUS-MINUS'
        write(6,89) 'ASL pm(-2) =',ASLpm(-2)
        write(6,89) 'ASL pm(-1) =',ASLpm(-1)
        write(6,89) 'ASL pm( 0) =',ASLpm( 0)
        write(6,89) 
c        write(6,89) 'ASL: MINUS-MINUS'
        write(6,89) 'ASL mm(-2) =',ASLmm(-2)
        write(6,89) 'ASL mm(-1) =',ASLmm(-1)
        write(6,89) 'ASL mm( 0) =',ASLmm( 0)
        write(6,89) 
c        write(6,89) 'ASL: PLUS-PLUS'
        write(6,89) 'ASL pp(-2) =',ASLpp(-2)
        write(6,89) 'ASL pp(-1) =',ASLpp(-1)
        write(6,89) 'ASL pp( 0) =',ASLpp( 0)
        write(6,89)       
      endif

      if (numcheck) then
        atreesl=atree('sl',k3,k2,k1,k4,k5,k6,zb,za)
        write(6,*) 'a6treemp',a6treemp
        write(6,*) 'atreesl ',atreesl
        call catanisl(mqsq,k1,k2,k3,k4,k5,k6,coeff2,coeff1)
        write(6,*) 'Catani(-2)*tree',coeff2*atreesl
        write(6,*) 'Catani(-1)*tree',coeff1*atreesl
        write(6,*)
        write(6,*) 'ASLmp(-2)/tree/Catani',ASLmp(-2)/atreesl/coeff2
        write(6,*) 'ASLmp(-1)/tree/Catani',ASLmp(-1)/atreesl/coeff1
        write(6,*)
        write(6,*) 'ASLpm(-2)/tree/Catani',ASLpm(-2)/a6treepm/coeff2
        write(6,*) 'ASLpm(-1)/tree/Catani',ASLpm(-1)/a6treepm/coeff1
        write(6,*)
        write(6,*) 'ASLmm(-2)/tree/Catani',ASLmm(-2)/a6treemm/coeff2
        write(6,*) 'ASLmm(-1)/tree/Catani',ASLmm(-1)/a6treemm/coeff1
        write(6,*)
        write(6,*) 'ASLpp(-2)/tree/Catani',ASLpp(-2)/a6treepp/coeff2
        write(6,*) 'ASLpp(-1)/tree/Catani',ASLpp(-1)/a6treepp/coeff1
        write(6,*)
      endif
      
      mhloopsq=mt**2
c--- light fermion loop contribution Af and heavy fermion contribution Ah
c--- (computed for a fermion of mass-squared mhloopsq, passed into the routine)
      call Afh(k1,k2,k3,k4,k5,k6,mq,Afmm,Afmp,Afpm,Afpp,
     &                              Ahmm,Ahmp,Ahpm,Ahpp,mhloopsq)

      if (numcheck) then
c        write(6,89) ' Af: MINUS-MINUS' 
        write(6,89) ' Af mm(-2) =',Afmm(-2)
        write(6,89) ' Af mm(-1) =',Afmm(-1)
        write(6,89) ' Af mm( 0) =',Afmm( 0)
        write(6,89) 
      
c        write(6,89) ' Af: MINUS-PLUS' 
        write(6,89) ' Af mp(-2) =',Afmp(-2)
        write(6,89) ' Af mp(-1) =',Afmp(-1)
        write(6,89) ' Af mp( 0) =',Afmp( 0)
        write(6,89) 
      
c        write(6,89) ' Af: PLUS-MINUS' 
        write(6,89) ' Af pm(-2) =',Afpm(-2)
        write(6,89) ' Af pm(-1) =',Afpm(-1)
        write(6,89) ' Af pm( 0) =',Afpm( 0)
        write(6,89) 
      
c        write(6,89) ' Af: PLUS-PLUS' 
        write(6,89) ' Af pp(-2) =',Afpp(-2)
        write(6,89) ' Af pp(-1) =',Afpp(-1)
        write(6,89) ' Af pp( 0) =',Afpp( 0)
        write(6,89) 
      
c        write(6,89) ' Ah: MINUS-MINUS' 
        write(6,89) ' Ah mm(-2) =',Ahmm(-2)
        write(6,89) ' Ah mm(-1) =',Ahmm(-1)
        write(6,89) ' Ah mm( 0) =',Ahmm( 0)
        write(6,89) 
      
c        write(6,89) ' Ah: MINUS-PLUS' 
        write(6,89) ' Ah mp(-2) =',Ahmp(-2)
        write(6,89) ' Ah mp(-1) =',Ahmp(-1)
        write(6,89) ' Ah mp( 0) =',Ahmp( 0)
        write(6,89) 
      
c        write(6,89) ' Ah: PLUS-MINUS' 
        write(6,89) ' Ah pm(-2) =',Ahpm(-2)
        write(6,89) ' Ah pm(-1) =',Ahpm(-1)
        write(6,89) ' Ah pm( 0) =',Ahpm( 0)
        write(6,89) 
      
c        write(6,89) ' Ah: PLUS-PLUS' 
        write(6,89) ' Ah pp(-2) =',Ahpp(-2)
        write(6,89) ' Ah pp(-1) =',Ahpp(-1)
        write(6,89) ' Ah pp( 0) =',Ahpp( 0)
        write(6,89) 
      endif
      
c      call writetable(
c     & A6treemp,A6treepm,A6treemm,A6treepp,
c     & ALCmp,ALCpm,ALCmm,ALCpp,
c     & ASLmp,ASLpm,ASLmm,ASLpp,
c     & Afmp,Afpm,Afmm,Afpp)
      
c--- we have now computed all the ingredients necessary for computing A6;1,
c--- so now add up according to the color sum
c--- Note: no heavy quark contribution included for now
      do iep=-2,0
      a61mmep(iep)=
     & ALCmm(iep)-two/xn**2*(BLCmm(iep)+ALCmm(iep))
     &  -1._dp/xn**2*ASLmm(iep)-real(nflav,dp)/xn*Afmm(iep)
     &  -1._dp/xn*Ahmm(iep)

      a61mpep(iep)=
     & ALCmp(iep)-two/xn**2*(BLCmp(iep)+ALCmp(iep))
     &  -1._dp/xn**2*ASLmp(iep)-real(nflav,dp)/xn*Afmp(iep)
     &  -1._dp/xn*Ahmp(iep)

      a61pmep(iep)=
     & ALCpm(iep)-two/xn**2*(BLCpm(iep)+ALCpm(iep))
     &  -1._dp/xn**2*ASLpm(iep)-real(nflav,dp)/xn*Afpm(iep)
     &  -1._dp/xn*Ahpm(iep)

      a61ppep(iep)=
     & ALCpp(iep)-two/xn**2*(BLCpp(iep)+ALCpp(iep))
     &  -1._dp/xn**2*ASLpp(iep)-real(nflav,dp)/xn*Afpp(iep)
     &  -1._dp/xn*Ahpp(iep)
      enddo

c--- add up coefficients with appropriate powers of epinv      
      a61mm=a61mmep(-2)*epinv**2+a61mmep(-1)*epinv+a61mmep(0)
      a61mp=a61mpep(-2)*epinv**2+a61mpep(-1)*epinv+a61mpep(0)
      a61pm=a61pmep(-2)*epinv**2+a61pmep(-1)*epinv+a61pmep(0)
      a61pp=a61ppep(-2)*epinv**2+a61ppep(-1)*epinv+a61ppep(0)

c--- overall wave function and couping constant renormalization
      ren=
     & +two*((real(nflav,dp)/3._dp-11._dp/6._dp*xn)*epinv+xn/6._dp)
     & +two/3._dp*(epinv+log(musq/mhloopsq))
     & -(xn**2-1._dp)/two/xn*(3._dp*(epinv+log(musq/mqsq))+5._dp)
      

c--- overall normalization of amplitude wrt. t-tbar process
      ren=ren/xn
      
      a61mm=a61mm+ren*a6treemm
      a61mp=a61mp+ren*a6treemp
      a61pm=a61pm+ren*a6treepm
      a61pp=a61pp+ren*a6treepp

c      write(6,*) 'a61mmep(-1),real(nflav,dp)*two/3._dp',
c     & a61mmep(-1),real(nflav,dp)*two/3._dp*a6treemm/xn
c      write(6,*) 'a61mpep(-1),real(nflav,dp)*two/3._dp',
c     & a61mpep(-1),real(nflav,dp)*two/3._dp*a6treemp/xn
c      write(6,*) 'a61pmep(-1),real(nflav,dp)*two/3._dp',
c     & a61pmep(-1),real(nflav,dp)*two/3._dp*a6treepm/xn
c      write(6,*) 'a61ppep(-1),real(nflav,dp)*two/3._dp',
c     & a61ppep(-1),real(nflav,dp)*two/3._dp*a6treepp/xn
c      pause

      return

   89 format(SP,a13,2e20.11,' i')   
      
      end
      
      
