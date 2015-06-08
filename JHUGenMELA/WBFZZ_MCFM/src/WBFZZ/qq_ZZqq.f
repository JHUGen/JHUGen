      subroutine qq_ZZqq(p,msq,zzcoupl,wwcoupl,lambda_bsm,lambdaq,
     . lambda)
      implicit none
c--- Author: R.K. Ellis, October 2014
c--- q(-p1)+q(-p2)->Z(p3,p4)+Z(p5,p6)+q(p7)+q(p8);
      include 'spinzerohiggs_anomcoupl.f' !--F
      include 'constants.f'
      include 'cmplxmass.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'runstring.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f'
      include 'first.f'
      include 'WWbits.f'
      double complex zzcoupl(32),wwcoupl(32) !--F
      double precision lambdaq,lambda_bsm,lambda(4)!--F
      integer nmax,jmax
      parameter(jmax=12,nmax=10)
      integer j,k,ll,
     & uqcq_uqcq,uquq_uquq,dqsq_dqsq,
     & dqdq_dqdq,uqbq_uqbq,dqcq_dqcq,
     & dquq_dquq,dqcq_uqsq,uqsq_dqcq
      parameter(
     & uqcq_uqcq=1,uquq_uquq=2,dqsq_dqsq=3,
     & dqdq_dqdq=4,uqbq_uqbq=5,dqcq_dqcq=6,
     & dquq_dquq=7,dqcq_uqsq=8,uqsq_dqcq=9)
      integer h1,h2,h3,h5
      double precision p(mxpart,4),msq(fn:nf,fn:nf),temp(fn:nf,fn:nf),
     & tempw(fn:nf,fn:nf),stat,spinavge,mult
      double complex zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),cdotpr,
     & propw71,propw81,propw72,propw82,
     & propw7341,propw7561,propw7342,propw7562,
     & ZZ7341(2,2,2,2,2,2),ZZ7561(2,2,2,2,2,2),
     & ZZ8341(2,2,2,2,2,2),ZZ8561(2,2,2,2,2,2),
     & WWp7341(2),WWm7341(2),WWp7561(2),WWm7561(2),
     & WWp8341(2),WWm8341(2),WWp8561(2),WWm8561(2),
     & j7_1(4,2),j7_2(4,2),j8_1(4,2),j8_2(4,2),
     & j7_34_1(4,2,2,2),j7_34_2(4,2,2,2),
     & j7_56_1(4,2,2,2),j7_56_2(4,2,2,2),
     & j8_34_1(4,2,2,2),j8_34_2(4,2,2,2),
     & j8_56_1(4,2,2,2),j8_56_2(4,2,2,2),
     & jl7_34_1(4,2,2,2),jl7_34_2(4,2,2,2),
     & jl7_56_1(4,2,2,2),jl7_56_2(4,2,2,2),
     & jl8_34_1(4,2,2,2),jl8_34_2(4,2,2,2),
     & jl8_56_1(4,2,2,2),jl8_56_2(4,2,2,2),
     & jw7_34_1(4,2,2),jw7_34_2(4,2,2),
     & jw7_56_1(4,2,2),jw7_56_2(4,2,2),
     & jw8_34_1(4,2,2),jw8_34_2(4,2,2),
     & jw8_56_1(4,2,2),jw8_56_2(4,2,2),
     & j7_3456_1(4,2,2,2,2),j8_3456_2(4,2,2,2,2),
     & j7_3456_2(4,2,2,2,2),j8_3456_1(4,2,2,2,2),
     & jw7_3456_1(4,2,2,2),jw8_3456_2(4,2,2,2),
     & jw7_3456_2(4,2,2,2),jw8_3456_1(4,2,2,2),
     & amp(nmax,2,2,2,2),ampa(nmax,2,2,2,2),ampb(nmax,2,2,2,2),
     & ZZHamp71_82(2,2,2,2,2,2),ZZHamp81_72(2,2,2,2,2,2),
     & WWZZ71_82amp(2,2),WWZZ81_72amp(2,2),
     & gmZ7341(2,2,2,2),gmZ7561(2,2,2,2),gmZ71(2,2,2,2),gmZ82(2,2,2,2),
     & gmZ7342(2,2,2,2),gmZ7562(2,2,2,2),gmZ72(2,2,2,2),gmZ81(2,2,2,2),
     & ll7341(2,2),ll7561(2,2),ll7342(2,2),ll7562(2,2),
     & gmZl7341(2,2,2),gmZl7561(2,2,2),gmZl7342(2,2,2),gmZl7562(2,2,2),
     & gmZl8562(2,2,2),gmZl8342(2,2,2),gmZl8561(2,2,2),gmZl8341(2,2,2),
     & k7341(4),k1567(4),k7342(4),k7562(4),ggWW(2,2),
     & srWWZZ71_82amp(2,2),srWWZZ81_72amp(2,2)
c     ,j3_4(4,2),j5_6(4,2),
      logical doHO,doBO
      parameter(spinavge=0.25d0,stat=0.5d0)
      integer,parameter:: j1(jmax)=(/1,2,8,8,7,2,7,1,1,7,2,7/)
      integer,parameter:: j2(jmax)=(/2,1,7,7,2,7,1,7,7,1,7,2/)
      integer,parameter:: j7(jmax)=(/7,7,2,1,1,8,2,8,2,8,1,8/)
      integer,parameter:: j8(jmax)=(/8,8,1,2,8,1,8,2,8,2,8,1/)
      save doHO,doBO,mult

      data Q(-5)/+0.333333333333333d0/
      data Q(-4)/-0.666666666666667d0/
      data Q(-3)/+0.333333333333333d0/
      data Q(-2)/-0.666666666666667d0/
      data Q(-1)/+0.333333333333333d0/
      data Q(0)/+0d0/
      data Q(+1)/-0.333333333333333d0/
      data Q(+2)/+0.666666666666667d0/
      data Q(+3)/-0.333333333333333d0/
      data Q(+4)/+0.666666666666667d0/
      data Q(+5)/-0.333333333333333d0/
      data tau/1d0,-1d0,1d0,-1d0,1d0,0d0,-1d0,1d0,-1d0,1d0,-1d0/
      data mt,twidth/173.2d0,2.5d0/
      data hmass,hwidth/125d0,0.00415d0/
      data wmass,wwidth/80.39d0,2.085d0/
      data zmass,zwidth/91.19d0,2.4952d0/
      data Gf,vevsq/1.16639d-5,246d0/
      data gw,xw,gwsq,esq/0.42464d0,0.23119d0,0.1d0,0.3133285d0/

      lambdaBSM = lambda_BSM
      lambda_Q = lambdaQ

      lambda_z1 = lambda(1)
      lambda_z2 = lambda(2)
      lambda_z3 = lambda(3)
      lambda_z4 = lambda(4)

      !-- coupling 1-4 and lambdas are missing
      ghz1 = zzcoupl(1) !-- check that these are correct according to Markus' definition
      ghz2 = zzcoupl(2)
      ghz3 = zzcoupl(3)
      ghz4 = zzcoupl(4)

      ghz1_prime = zzcoupl(5) !-- from here on taken from Markus
      ghz1_prime2= zzcoupl(6) 
      ghz1_prime3= zzcoupl(7) 
      ghz1_prime4= zzcoupl(8)
      ghz1_prime5= zzcoupl(9)
      
      ghz2_prime = zzcoupl(10) 
      ghz2_prime2= zzcoupl(11)
      ghz2_prime3= zzcoupl(12)
      ghz2_prime4= zzcoupl(13)
      ghz2_prime5= zzcoupl(14)
      
      ghz3_prime = zzcoupl(15)
      ghz3_prime2= zzcoupl(16)
      ghz3_prime3= zzcoupl(17)
      ghz3_prime4= zzcoupl(18)
      ghz3_prime5= zzcoupl(19)
      
      ghz4_prime = zzcoupl(20)
      ghz4_prime2= zzcoupl(21)
      ghz4_prime3= zzcoupl(22)
      ghz4_prime4= zzcoupl(23)
      ghz4_prime5= zzcoupl(24)
      
      ghz1_prime6= zzcoupl(25)
      ghz1_prime7= zzcoupl(26)
      
      ghz2_prime6= zzcoupl(27)
      ghz2_prime7= zzcoupl(28)
      
      ghz3_prime6= zzcoupl(29)
      ghz3_prime7= zzcoupl(30)
      
      ghz4_prime6= zzcoupl(31)
      ghz4_prime7= zzcoupl(32)
      
      ghz1_prime = zzcoupl(5) 
      ghz1_prime2= zzcoupl(6) 
      ghz1_prime3= zzcoupl(7) 
      ghz1_prime4= zzcoupl(8)
      ghz1_prime5= zzcoupl(9)
      
      ghz2_prime = zzcoupl(10) 
      ghz2_prime2= zzcoupl(11)
      ghz2_prime3= zzcoupl(12)
      ghz2_prime4= zzcoupl(13)
      ghz2_prime5= zzcoupl(14)
      
      ghz3_prime = zzcoupl(15)
      ghz3_prime2= zzcoupl(16)
      ghz3_prime3= zzcoupl(17)
      ghz3_prime4= zzcoupl(18)
      ghz3_prime5= zzcoupl(19)
      
      ghz4_prime = zzcoupl(20)
      ghz4_prime2= zzcoupl(21)
      ghz4_prime3= zzcoupl(22)
      ghz4_prime4= zzcoupl(23)
      ghz4_prime5= zzcoupl(24)
      
      ghz1_prime6= zzcoupl(25)
      ghz1_prime7= zzcoupl(26)
      
      ghz2_prime6= zzcoupl(27)
      ghz2_prime7= zzcoupl(28)
      
      ghz3_prime6= zzcoupl(29)
      ghz3_prime7= zzcoupl(30)
      
      ghz4_prime6= zzcoupl(31)
      ghz4_prime7= zzcoupl(32)


      ghw1 = wwcoupl(1) !-- check that these are correct according to Markus' definition
      ghw2 = wwcoupl(2)
      ghw3 = wwcoupl(3)
      ghw4 = wwcoupl(4)

      ghw1_prime = wwcoupl(5) 
      ghw1_prime2= wwcoupl(6) 
      ghw1_prime3= wwcoupl(7) 
      ghw1_prime4= wwcoupl(8)
      ghw1_prime5= wwcoupl(9)
      
      ghw2_prime = wwcoupl(10) 
      ghw2_prime2= wwcoupl(11)
      ghw2_prime3= wwcoupl(12)
      ghw2_prime4= wwcoupl(13)
      ghw2_prime5= wwcoupl(14)
      
      ghw3_prime = wwcoupl(15)
      ghw3_prime2= wwcoupl(16)
      ghw3_prime3= wwcoupl(17)
      ghw3_prime4= wwcoupl(18)
      ghw3_prime5= wwcoupl(19)
      
      ghw4_prime = wwcoupl(20)
      ghw4_prime2= wwcoupl(21)
      ghw4_prime3= wwcoupl(22)
      ghw4_prime4= wwcoupl(23)
      ghw4_prime5= wwcoupl(24)
      
      ghw1_prime6= wwcoupl(25)
      ghw1_prime7= wwcoupl(26)
      
      ghw2_prime6= wwcoupl(27)
      ghw2_prime7= wwcoupl(28)
      
      ghw3_prime6= wwcoupl(29)
      ghw3_prime7= wwcoupl(30)
      
      ghw4_prime6= wwcoupl(31)
      ghw4_prime7= wwcoupl(32)
      
      ghw1_prime = wwcoupl(5) 
      ghw1_prime2= wwcoupl(6) 
      ghw1_prime3= wwcoupl(7) 
      ghw1_prime4= wwcoupl(8)
      ghw1_prime5= wwcoupl(9)
      
      ghw2_prime = wwcoupl(10) 
      ghw2_prime2= wwcoupl(11)
      ghw2_prime3= wwcoupl(12)
      ghw2_prime4= wwcoupl(13)
      ghw2_prime5= wwcoupl(14)
      
      ghw3_prime = wwcoupl(15)
      ghw3_prime2= wwcoupl(16)
      ghw3_prime3= wwcoupl(17)
      ghw3_prime4= wwcoupl(18)
      ghw3_prime5= wwcoupl(19)
      
      ghw4_prime = wwcoupl(20)
      ghw4_prime2= wwcoupl(21)
      ghw4_prime3= wwcoupl(22)
      ghw4_prime4= wwcoupl(23)
      ghw4_prime5= wwcoupl(24)
      
      ghw1_prime6= wwcoupl(25)
      ghw1_prime7= wwcoupl(26)
      
      ghw2_prime6= wwcoupl(27)
      ghw2_prime7= wwcoupl(28)
      
      ghw3_prime6= wwcoupl(29)
      ghw3_prime7= wwcoupl(30)
      
      ghw4_prime6= wwcoupl(31)
      ghw4_prime7= wwcoupl(32)
      
      
!$omp threadprivate(doHO,doBO,mult)      
      msq(:,:)=0d0

c--- This calculation uses the complex-mass scheme (c.f. arXiv:hep-ph/0605312)
c--- and the following lines set up the appropriate masses and sin^2(theta_w)

      if (first) then
       cwmass2=dcmplx(wmass**2,-wmass*wwidth)
       czmass2=dcmplx(zmass**2,-zmass*zwidth)
       cxw=cone-cwmass2/czmass2
       call couplz(xw)
       
       ! MARKUS: choosing ZZ final state
       l1=le
       l2=le
       r1=le
       r2=le
       
c       cxw=dcmplx(xw,0d0) ! DEBUG: Madgraph comparison
       write(6,*)
       write(6,*) '**************** Complex-mass scheme ***************'
       write(6,*) '*                                                  *'
       write(6,77) cwmass2
       write(6,78) czmass2
       write(6,79) cxw
       write(6,*) '*                                                  *'
       write(6,*) '****************************************************'
       write(6,*)
       doHO=.false.
       doBO=.false.
                 runstring(4:5)="HO"
       if     (runstring(4:5) .eq. 'HO') then
         doHO=.true.
       write(6,*) '>>>>>>>>>>>>>> Higgs contribution only <<<<<<<<<<<<<'
       write(6,*)
       elseif (runstring(4:5) .eq. 'BO') then
         doBO=.true.
       write(6,*)
       write(6,*) '>>>>>>>>>>> Background contribution only <<<<<<<<<<<'
       write(6,*)
       endif
       mult=1d0
c--- rescaling factor for Higgs amplitudes, if anomalous Higgs width
       if (anom_Higgs) then
         mult=chi_higgs**2
       endif
       first=.false.
       call flush(6)
      endif

      if (doHO) then
        Hbit=mult*cone
        Bbit=czip
      elseif (doBO) then
        Hbit=czip
        Bbit=cone
      else
        Hbit=mult*cone
        Bbit=cone
      endif

C---setup spinors and spinorvector products
      call spinorcurr(8,p,za,zb,zab,zba)

      do j=1,jmax
      temp(:,:)=0d0
      tempw(:,:)=0d0
      amp(:,:,:,:,:)=czip
      ampa(:,:,:,:,:)=czip
      ampb(:,:,:,:,:)=czip

c--- propagators and currents are not used in calculation of Higgs contribution
      if (doHO .eqv. .false.) then
      call setupzprops(j1(j),j2(j),3,4,5,6,j7(j),j8(j),
     & gmZ7341,gmZ7561,gmZ71,gmZ82,
     & gmZ7342,gmZ7562,gmZ72,gmZ81,
     & ggWW,propw71,propw81,propw72,propw82,
     & propw7341,propw7561,propw7342,propw7562,
     & ll7341,ll7561,ll7342,ll7562,
     & gmZl7341,gmZl7561,gmZl7342,gmZl7562,
     & gmZl8562,gmZl8342,gmZl8561,gmZl8341)

      call jzero(j7(j),j1(j),zab,zba,j7_1)
      call jzero(j7(j),j2(j),zab,zba,j7_2)
      call jzero(j8(j),j1(j),zab,zba,j8_1)
      call jzero(j8(j),j2(j),zab,zba,j8_2)
 
      call jone(j7(j),3,4,j1(j),za,zb,zab,zba,j7_34_1,jw7_34_1,jl7_34_1)
      call jone(j7(j),3,4,j2(j),za,zb,zab,zba,j7_34_2,jw7_34_2,jl7_34_2)
      call jone(j7(j),5,6,j1(j),za,zb,zab,zba,j7_56_1,jw7_56_1,jl7_56_1)
      call jone(j7(j),5,6,j2(j),za,zb,zab,zba,j7_56_2,jw7_56_2,jl7_56_2)
      call jone(j8(j),3,4,j1(j),za,zb,zab,zba,j8_34_1,jw8_34_1,jl8_34_1)
      call jone(j8(j),3,4,j2(j),za,zb,zab,zba,j8_34_2,jw8_34_2,jl8_34_2)
      call jone(j8(j),5,6,j1(j),za,zb,zab,zba,j8_56_1,jw8_56_1,jl8_56_1)
      call jone(j8(j),5,6,j2(j),za,zb,zab,zba,j8_56_2,jw8_56_2,jl8_56_2)

      call jtwo(j7(j),3,4,5,6,j1(j),za,zb,zab,zba,j7_3456_1,jw7_3456_1)
      call jtwo(j7(j),3,4,5,6,j2(j),za,zb,zab,zba,j7_3456_2,jw7_3456_2)
      call jtwo(j8(j),3,4,5,6,j1(j),za,zb,zab,zba,j8_3456_1,jw8_3456_1)
      call jtwo(j8(j),3,4,5,6,j2(j),za,zb,zab,zba,j8_3456_2,jw8_3456_2)

      k7341(:)=0.5d0*(zab(j1(j),:,j1(j))+zab(3,:,3)
     & +zab(4,:,4)+zab(j7(j),:,j7(j)))
      k1567(:)=0.5d0*(zab(j1(j),:,j1(j))+zab(5,:,5)
     & +zab(6,:,6)+zab(j7(j),:,j7(j)))
      k7342(:)=0.5d0*(zab(j2(j),:,j2(j))+zab(3,:,3)
     & +zab(4,:,4)+zab(j7(j),:,j7(j)))
      k7562(:)=0.5d0*(zab(j2(j),:,j2(j))+zab(5,:,5)
     & +zab(6,:,6)+zab(j7(j),:,j7(j)))

C-----Singly resonant production in VBF style diagrams
      call ZZSingleres(j1(j),j2(j),3,4,5,6,j7(j),j8(j),za,zb,
     & ZZ7341,WWp7341,WWm7341)
      call ZZSingleres(j1(j),j2(j),5,6,3,4,j7(j),j8(j),za,zb,
     & ZZ7561,WWp7561,WWm7561)
      call ZZSingleres(j1(j),j2(j),3,4,5,6,j8(j),j7(j),za,zb,
     & ZZ8341,WWp8341,WWm8341)
      call ZZSingleres(j1(j),j2(j),5,6,3,4,j8(j),j7(j),za,zb,
     & ZZ8561,WWp8561,WWm8561)
      endif
      
      
      
C----ZZ->ZZ scattering with the exchange of a H
      call ZZHZZamp(j1(j),j2(j),3,4,5,6,j7(j),j8(j),
     & za,zb,ZZHamp71_82)
      call ZZHZZamp(j1(j),j2(j),3,4,5,6,j8(j),j7(j),
     & za,zb,ZZHamp81_72)
C----Four boson vertex + WW->Higgs diagram 
      call WWZZ(j1(j),j2(j),3,4,5,6,j7(j),j8(j),
     & za,zb,WWZZ71_82amp,srWWZZ71_82amp) 
      call WWZZ(j1(j),j2(j),3,4,5,6,j8(j),j7(j),
     & za,zb,WWZZ81_72amp,srWWZZ81_72amp) 

C-----setup for (uqbq_uqbq) (2,5)->(2,5)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(uqbq_uqbq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,2,h1,h3),j8_56_2(:,1,h2,h5))*gmZ7341(2,1,h1,h2)
     & +cdotpr(jl7_34_1(:,2,h1,h3),jl8_56_2(:,1,h2,h5))*ll7341(h3,h5)
     & +cdotpr(j7_34_1(:,2,h1,h3),jl8_56_2(:,1,h2,h5))*gmZl8562(2,h1,h5)
     & +cdotpr(jl7_34_1(:,2,h1,h3),j8_56_2(:,1,h2,h5))*gmZl7341(1,h2,h3)
     
      amp(uqbq_uqbq,h1,h2,h3,h5)=amp(uqbq_uqbq,h1,h2,h3,h5)
     & +cdotpr(j7_56_1(:,2,h1,h5),j8_34_2(:,1,h2,h3))*gmZ7561(2,1,h1,h2)
     & +cdotpr(jl7_56_1(:,2,h1,h5),jl8_34_2(:,1,h2,h3))*ll7561(h3,h5)
     & +cdotpr(j7_56_1(:,2,h1,h5),jl8_34_2(:,1,h2,h3))*gmZl8342(2,h1,h3)
     & +cdotpr(jl7_56_1(:,2,h1,h5),j8_34_2(:,1,h2,h3))*gmZl7561(1,h2,h5)

      amp(uqbq_uqbq,h1,h2,h3,h5)=amp(uqbq_uqbq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,2,h1,h3,h5),j8_2(:,h2))*gmZ82(2,1,h1,h2)
     & +cdotpr(j7_1(:,h1),j8_3456_2(:,1,h2,h3,h5))*gmZ71(2,1,h1,h2)

      amp(uqbq_uqbq,h1,h2,h3,h5)=amp(uqbq_uqbq,h1,h2,h3,h5)
     & +ZZ7341(2,1,h1,h2,h3,h5)+ZZ7561(2,1,h1,h2,h5,h3)
      endif

      amp(uqbq_uqbq,h1,h2,h3,h5)=Bbit*amp(uqbq_uqbq,h1,h2,h3,h5)
     & +Hbit*ZZHamp71_82(2,1,h1,h2,h3,h5)

      temp(2,5)=temp(2,5)+esq**6*spinavge
     &   *dble(amp(uqbq_uqbq,h1,h2,h3,h5)
     & *dconjg(amp(uqbq_uqbq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      temp(4,5)=temp(2,5)
C--------------------------------------------------------------------------
C-----setup for (uqcq_uqcq) (2,4)->(2,4)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2

      if (doHO .eqv. .false.) then
      amp(uqcq_uqcq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,2,h1,h3),j8_56_2(:,2,h2,h5))*gmZ7341(2,2,h1,h2)
     & +cdotpr(jl7_34_1(:,2,h1,h3),jl8_56_2(:,2,h2,h5))*ll7341(h3,h5)
     & +cdotpr(j7_34_1(:,2,h1,h3),jl8_56_2(:,2,h2,h5))*gmZl8562(2,h1,h5)
     & +cdotpr(jl7_34_1(:,2,h1,h3),j8_56_2(:,2,h2,h5))*gmZl7341(2,h2,h3)

      amp(uqcq_uqcq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)
     & +cdotpr(j7_56_1(:,2,h1,h5),j8_34_2(:,2,h2,h3))*gmZ7561(2,2,h1,h2)
     & +cdotpr(jl7_56_1(:,2,h1,h5),jl8_34_2(:,2,h2,h3))*ll7561(h3,h5)
     & +cdotpr(j7_56_1(:,2,h1,h5),jl8_34_2(:,2,h2,h3))*gmZl8342(2,h1,h3)
     & +cdotpr(jl7_56_1(:,2,h1,h5),j8_34_2(:,2,h2,h3))*gmZl7561(2,h2,h5)
      
      amp(uqcq_uqcq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,2,h1,h3,h5),j8_2(:,h2))*gmZ82(2,2,h1,h2)
     & +cdotpr(j7_1(:,h1),j8_3456_2(:,2,h2,h3,h5))*gmZ71(2,2,h1,h2)

      amp(uqcq_uqcq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)
     & +ZZ7341(2,2,h1,h2,h3,h5)+ZZ7561(2,2,h1,h2,h5,h3)
      endif

      amp(uqcq_uqcq,h1,h2,h3,h5)=Bbit*amp(uqcq_uqcq,h1,h2,h3,h5)
     & +Hbit*ZZHamp71_82(2,2,h1,h2,h3,h5)

      temp(2,4)=temp(2,4)+esq**6*spinavge
     &   *dble(amp(uqcq_uqcq,h1,h2,h3,h5)
     & *dconjg(amp(uqcq_uqcq,h1,h2,h3,h5)))

c      if (j .eq. 1) then
c      write(6,*) h1,h2,h3,h5,esq**6*spinavge
c     &   *dble(amp(uqcq_uqcq,h1,h2,h3,h5)
c     & *dconjg(amp(uqcq_uqcq,h1,h2,h3,h5)))
c      endif
      
      enddo
      enddo
      enddo
      enddo


C-----setup for uqsq_dqcq W diagrams (2,3)->(1,4)
      do h1=1,1
      do h2=1,1
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(uqsq_dqcq,h1,h2,h3,h5)=
     & +cdotpr(jw7_3456_1(:,2,h3,h5),j8_2(:,h2))*0.5d0/propw82/cxw

      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +cdotpr(j7_1(:,h1),jw8_3456_2(:,1,h3,h5))*0.5d0/propw71/cxw

      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +(cdotpr(jw7_34_1(:,2,h3),jw8_56_2(:,1,h5))
     &  -cdotpr(jw7_34_1(:,2,h3),k7341(:))
     &  *cdotpr(k7341(:),jw8_56_2(:,1,h5))/cwmass2)/propw7341

      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +(cdotpr(jw7_56_1(:,2,h5),jw8_34_2(:,1,h3))
     &  -cdotpr(jw7_56_1(:,2,h5),k1567(:))
     &  *cdotpr(k1567(:),jw8_34_2(:,1,h3))/cwmass2)/propw7561

      if (h3 .eq. 1)
     & amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +WWp7341(h5)
      if (h5 .eq. 1)
     & amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +WWp7561(h3)

      amp(uqsq_dqcq,h1,h2,h3,h5)=amp(uqsq_dqcq,h1,h2,h3,h5)
     & +srWWZZ71_82amp(h3,h5)
      endif

      amp(uqsq_dqcq,h1,h2,h3,h5)=Bbit*amp(uqsq_dqcq,h1,h2,h3,h5)
     & +WWZZ71_82amp(h3,h5)

      tempw(2,3)=tempw(2,3)+esq**6*spinavge
     &   *dble(amp(uqsq_dqcq,h1,h2,h3,h5)
     & *dconjg(amp(uqsq_dqcq,h1,h2,h3,h5)))
      enddo     
      enddo     
      enddo     
      enddo     
      temp(2,3)=temp(2,5)

C-----setup for dqcq_uqsq (1,4)-->(2,3)
      do h1=1,1
      do h2=1,1
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(dqcq_uqsq,h1,h2,h3,h5)=
     & +cdotpr(jw7_3456_1(:,1,h3,h5),j8_2(:,h2))*0.5d0/propw82/cxw

      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +cdotpr(j7_1(:,h1),jw8_3456_2(:,2,h3,h5))*0.5d0/propw71/cxw

      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +(cdotpr(jw7_34_1(:,1,h3),jw8_56_2(:,2,h5))
     &  -cdotpr(jw7_34_1(:,1,h3),k7341(:))
     &  *cdotpr(k7341(:),jw8_56_2(:,2,h5))/cwmass2)/propw7341

      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +(cdotpr(jw7_56_1(:,1,h5),jw8_34_2(:,2,h3))
     &  -cdotpr(jw7_56_1(:,1,h5),k1567(:))
     &  *cdotpr(k1567(:),jw8_34_2(:,2,h3))/cwmass2)/propw7561

      if (h3 .eq. 1)
     & amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +WWm7341(h5)
      if (h5 .eq. 1)
     & amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & +WWm7561(h3)

      amp(dqcq_uqsq,h1,h2,h3,h5)=amp(dqcq_uqsq,h1,h2,h3,h5)
     & -srWWZZ71_82amp(h3,h5) ! note minus sign instead of exchanging 1<->7,2<->8
      endif

      amp(dqcq_uqsq,h1,h2,h3,h5)=Bbit*amp(dqcq_uqsq,h1,h2,h3,h5)
     & +WWZZ71_82amp(h3,h5)

      tempw(1,4)=tempw(1,4)+esq**6*spinavge
     &   *dble(amp(dqcq_uqsq,h1,h2,h3,h5)
     & *dconjg(amp(dqcq_uqsq,h1,h2,h3,h5)))
      enddo     
      enddo     
      enddo     
      enddo     

C-----setup for (dqcq_dqcq) (1,4)-->(1,4)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(dqcq_dqcq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,1,h1,h3),j8_56_2(:,2,h2,h5))*gmZ7341(1,2,h1,h2)
     & +cdotpr(jl7_34_1(:,1,h1,h3),jl8_56_2(:,2,h2,h5))*ll7341(h3,h5)
     & +cdotpr(j7_34_1(:,1,h1,h3),jl8_56_2(:,2,h2,h5))*gmZl8562(1,h1,h5)
     & +cdotpr(jl7_34_1(:,1,h1,h3),j8_56_2(:,2,h2,h5))*gmZl7341(2,h2,h3)
      amp(dqcq_dqcq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)
     & +cdotpr(j7_56_1(:,1,h1,h5),j8_34_2(:,2,h2,h3))*gmZ7561(1,2,h1,h2)
     & +cdotpr(jl7_56_1(:,1,h1,h5),jl8_34_2(:,2,h2,h3))*ll7561(h3,h5)
     & +cdotpr(j7_56_1(:,1,h1,h5),jl8_34_2(:,2,h2,h3))*gmZl8342(1,h1,h3)
     & +cdotpr(jl7_56_1(:,1,h1,h5),j8_34_2(:,2,h2,h3))*gmZl7561(2,h2,h5)

      amp(dqcq_dqcq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,1,h1,h3,h5),j8_2(:,h2))*gmZ82(1,2,h1,h2)

      amp(dqcq_dqcq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)
     & +cdotpr(j7_1(:,h1),j8_3456_2(:,2,h2,h3,h5))*gmZ71(1,2,h1,h2)

      amp(dqcq_dqcq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)
     & +ZZ7341(1,2,h1,h2,h3,h5)+ZZ7561(1,2,h1,h2,h5,h3)
      endif

      amp(dqcq_dqcq,h1,h2,h3,h5)=Bbit*amp(dqcq_dqcq,h1,h2,h3,h5)
     & +Hbit*ZZHamp71_82(1,2,h1,h2,h3,h5)

      temp(1,4)=temp(1,4)+esq**6*spinavge
     &   *dble(amp(dqcq_dqcq,h1,h2,h3,h5)
     & *dconjg(amp(dqcq_dqcq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo

C--------------------------------------------------------
C-----setup for dquq_dquq W diagrams (1,2)-->(1,2)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      if ((h1.eq.1) .and. (h2.eq.1)) then
      if (doHO .eqv. .false.) then
      ampa(dquq_dquq,h1,h2,h3,h5)=
     & +cdotpr(jw8_3456_1(:,1,h3,h5),j7_2(:,h2))*0.5d0/propw72/cxw

      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +cdotpr(j8_1(:,h1),jw7_3456_2(:,2,h3,h5))*0.5d0/propw81/cxw

      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +(cdotpr(jw8_34_1(:,1,h3),jw7_56_2(:,2,h5))
     &  -cdotpr(jw8_34_1(:,1,h3),k7562(:))
     &  *cdotpr(k7562(:),jw7_56_2(:,2,h5))/cwmass2)/propw7562

      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +(cdotpr(jw8_56_1(:,1,h5),jw7_34_2(:,2,h3))
     &  -cdotpr(jw8_56_1(:,1,h5),k7342(:))
     &  *cdotpr(k7342(:),jw7_34_2(:,2,h3))/cwmass2)/propw7342

      if (h3 .eq. 1)
     & ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +WWm8341(h5)
      if (h5 .eq. 1)
     & ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & +WWm8561(h3)

      ampa(dquq_dquq,h1,h2,h3,h5)=ampa(dquq_dquq,h1,h2,h3,h5)
     & -srWWZZ81_72amp(h3,h5) ! note minus sign instead of exchanging 1<->7,2<->8
      endif
      
      ampa(dquq_dquq,h1,h2,h3,h5)=Bbit*ampa(dquq_dquq,h1,h2,h3,h5)
     & +WWZZ81_72amp(h3,h5)
      endif

C--Fill Z exchange diagrams
      ampb(dquq_dquq,h1,h2,h3,h5)=amp(dqcq_dqcq,h1,h2,h3,h5)
       
      temp(1,2)=temp(1,2)+esq**6*spinavge
     &   *dble(ampa(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampa(dquq_dquq,h1,h2,h3,h5)))
      temp(1,2)=temp(1,2)+esq**6*spinavge
     &   *dble(ampb(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampb(dquq_dquq,h1,h2,h3,h5)))
      temp(1,2)=temp(1,2)-2d0/xn*esq**6*spinavge
     &   *dble(ampa(dquq_dquq,h1,h2,h3,h5)
     & *dconjg(ampb(dquq_dquq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      temp(3,4)=temp(1,2)
C-----------------------------------------------------------------

C-----setup for (dqsq_dqsq) (1,3)-->(1,3)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
      if (doHO .eqv. .false.) then
      amp(dqsq_dqsq,h1,h2,h3,h5)=
     & +cdotpr(j7_34_1(:,1,h1,h3),j8_56_2(:,1,h2,h5))*gmZ7341(1,1,h1,h2)
     & +cdotpr(jl7_34_1(:,1,h1,h3),jl8_56_2(:,1,h2,h5))*ll7341(h3,h5)
     & +cdotpr(j7_34_1(:,1,h1,h3),jl8_56_2(:,1,h2,h5))*gmZl8562(1,h1,h5)
     & +cdotpr(jl7_34_1(:,1,h1,h3),j8_56_2(:,1,h2,h5))*gmZl7341(1,h2,h3)

      amp(dqsq_dqsq,h1,h2,h3,h5)=amp(dqsq_dqsq,h1,h2,h3,h5)
     & +cdotpr(j7_56_1(:,1,h1,h5),j8_34_2(:,1,h2,h3))*gmZ7561(1,1,h1,h2)
     & +cdotpr(jl7_56_1(:,1,h1,h5),jl8_34_2(:,1,h2,h3))*ll7561(h3,h5)
     & +cdotpr(j7_56_1(:,1,h1,h5),jl8_34_2(:,1,h2,h3))*gmZl8342(1,h1,h3)
     & +cdotpr(jl7_56_1(:,1,h1,h5),j8_34_2(:,1,h2,h3))*gmZl7561(1,h2,h5)

      amp(dqsq_dqsq,h1,h2,h3,h5)=amp(dqsq_dqsq,h1,h2,h3,h5)
     & +cdotpr(j7_3456_1(:,1,h1,h3,h5),j8_2(:,h2))*gmZ82(1,1,h2,h1)
     & +cdotpr(j8_3456_2(:,1,h2,h3,h5),j7_1(:,h1))*gmZ71(1,1,h1,h2)

      amp(dqsq_dqsq,h1,h2,h3,h5)=amp(dqsq_dqsq,h1,h2,h3,h5)
     & +ZZ7341(1,1,h1,h2,h3,h5)+ZZ7561(1,1,h1,h2,h5,h3)
      endif
     
      amp(dqsq_dqsq,h1,h2,h3,h5)=Bbit*amp(dqsq_dqsq,h1,h2,h3,h5)
     & +Hbit*ZZHamp71_82(1,1,h1,h2,h3,h5)
     
      temp(1,3)=temp(1,3)+esq**6*spinavge
     &   *dble(amp(dqsq_dqsq,h1,h2,h3,h5)
     & *dconjg(amp(dqsq_dqsq,h1,h2,h3,h5)))
      enddo
      enddo
      enddo
      enddo
      temp(1,5)=temp(1,3)
      temp(3,5)=temp(1,3)

      if ((j.eq.2).or.(j.eq.4).or.(j.eq.6).or.(j.eq.8)
     & .or.(j.eq.10).or.(j.eq.12)) go to 100
C-----setup for ((uquq_uquq)  (2,2)-->(2,2)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
C-----------------ampa
      ampa(uquq_uquq,h1,h2,h3,h5)=amp(uqcq_uqcq,h1,h2,h3,h5)

C-----------------ampb
      if (doHO .eqv. .false.) then
      ampb(uquq_uquq,h1,h2,h3,h5)=
     & +cdotpr(j8_34_1(:,2,h1,h3),j7_56_2(:,2,h2,h5))*gmZ7562(2,2,h1,h2)
     & +cdotpr(jl8_34_1(:,2,h1,h3),jl7_56_2(:,2,h2,h5))*ll7562(h3,h5)
     & +cdotpr(j8_34_1(:,2,h1,h3),jl7_56_2(:,2,h2,h5))*gmZl7562(2,h1,h5)
     & +cdotpr(jl8_34_1(:,2,h1,h3),j7_56_2(:,2,h2,h5))*gmZl8341(2,h2,h3)
      ampb(uquq_uquq,h1,h2,h3,h5)=ampb(uquq_uquq,h1,h2,h3,h5)
     & +cdotpr(j8_56_1(:,2,h1,h5),j7_34_2(:,2,h2,h3))*gmZ7342(2,2,h1,h2)
     & +cdotpr(jl8_56_1(:,2,h1,h5),jl7_34_2(:,2,h2,h3))*ll7342(h3,h5)
     & +cdotpr(j8_56_1(:,2,h1,h5),jl7_34_2(:,2,h2,h3))*gmZl7342(2,h1,h3)
     & +cdotpr(jl8_56_1(:,2,h1,h5),j7_34_2(:,2,h2,h3))*gmZl8561(2,h2,h5)

      ampb(uquq_uquq,h1,h2,h3,h5)=ampb(uquq_uquq,h1,h2,h3,h5)
     & +cdotpr(j8_3456_1(:,2,h1,h3,h5),j7_2(:,h2))*gmZ72(2,2,h1,h2)
     & +cdotpr(j8_1(:,h1),j7_3456_2(:,2,h2,h3,h5))*gmZ81(2,2,h1,h2)

      ampb(uquq_uquq,h1,h2,h3,h5)=ampb(uquq_uquq,h1,h2,h3,h5)
     & +ZZ8341(2,2,h1,h2,h3,h5)+ZZ8561(2,2,h1,h2,h5,h3)
      endif

      ampb(uquq_uquq,h1,h2,h3,h5)=Bbit*ampb(uquq_uquq,h1,h2,h3,h5)
     & +Hbit*ZZHamp81_72(2,2,h1,h2,h3,h5)

      temp(2,2)=temp(2,2)+esq**6*spinavge
     & *dble(ampa(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampa(uquq_uquq,h1,h2,h3,h5)))
      temp(2,2)=temp(2,2)+esq**6*spinavge
     & *dble(ampb(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampb(uquq_uquq,h1,h2,h3,h5)))
      if (h1 .eq. h2) then
      temp(2,2)=temp(2,2)-2d0/xn*esq**6*spinavge
     & *dble(ampa(uquq_uquq,h1,h2,h3,h5)
     & *dconjg(ampb(uquq_uquq,h1,h2,h3,h5)))
      endif

      enddo
      enddo
      enddo
      enddo
      temp(4,4)=temp(2,2)


C-----setup for ((dqdq_dqdq)  (1,1)-->(1,1)
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h5=1,2
     
C-----------------ampa
      ampa(dqdq_dqdq,h1,h2,h3,h5)=amp(dqsq_dqsq,h1,h2,h3,h5)

C-----------------ampb
      if (doHO .eqv. .false.) then
      ampb(dqdq_dqdq,h1,h2,h3,h5)=
     & +cdotpr(j8_34_1(:,1,h1,h3),j7_56_2(:,1,h2,h5))*gmZ7562(1,1,h1,h2)
     & +cdotpr(jl8_34_1(:,1,h1,h3),jl7_56_2(:,1,h2,h5))*ll7562(h3,h5)
     & +cdotpr(j8_34_1(:,1,h1,h3),jl7_56_2(:,1,h2,h5))*gmZl7562(1,h1,h5)
     & +cdotpr(jl8_34_1(:,1,h1,h3),j7_56_2(:,1,h2,h5))*gmZl8341(1,h2,h3)
      ampb(dqdq_dqdq,h1,h2,h3,h5)=ampb(dqdq_dqdq,h1,h2,h3,h5)
     & +cdotpr(j8_56_1(:,1,h1,h5),j7_34_2(:,1,h2,h3))*gmZ7342(1,1,h1,h2)
     & +cdotpr(jl8_56_1(:,1,h1,h5),jl7_34_2(:,1,h2,h3))*ll7342(h3,h5)
     & +cdotpr(j8_56_1(:,1,h1,h5),jl7_34_2(:,1,h2,h3))*gmZl7342(1,h1,h3)
     & +cdotpr(jl8_56_1(:,1,h1,h5),j7_34_2(:,1,h2,h3))*gmZl8561(1,h2,h5)

      ampb(dqdq_dqdq,h1,h2,h3,h5)=ampb(dqdq_dqdq,h1,h2,h3,h5)
     & +cdotpr(j8_3456_1(:,1,h1,h3,h5),j7_2(:,h2))*gmZ72(1,1,h1,h2)
     & +cdotpr(j8_1(:,h1),j7_3456_2(:,1,h2,h3,h5))*gmZ81(1,1,h1,h2)

      ampb(dqdq_dqdq,h1,h2,h3,h5)=ampb(dqdq_dqdq,h1,h2,h3,h5)
     & +ZZ8341(1,1,h1,h2,h3,h5)+ZZ8561(1,1,h1,h2,h5,h3)
      endif

      ampb(dqdq_dqdq,h1,h2,h3,h5)=Bbit*ampb(dqdq_dqdq,h1,h2,h3,h5)
     & +Hbit*ZZHamp81_72(1,1,h1,h2,h3,h5)

      temp(1,1)=temp(1,1)+esq**6*spinavge
     & *dble(ampa(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampa(dqdq_dqdq,h1,h2,h3,h5)))
      temp(1,1)=temp(1,1)+esq**6*spinavge
     & *dble(ampb(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampb(dqdq_dqdq,h1,h2,h3,h5)))
      if (h1 .eq. h2) then
      temp(1,1)=temp(1,1)-2d0/xn*esq**6*spinavge
     & *dble(ampa(dqdq_dqdq,h1,h2,h3,h5)
     & *dconjg(ampb(dqdq_dqdq,h1,h2,h3,h5)))
      endif

      enddo
      enddo
      enddo
      enddo
      temp(3,3)=temp(1,1)
      temp(5,5)=temp(1,1)
 100  continue

      if (j.eq.1) then
      do k=1,nf
      msq(k,k)=temp(k,k)*stat
      do ll=k+1,nf
      msq(k,ll)=temp(k,ll)
      enddo
      enddo
      msq(2,3)=msq(2,3)+tempw(2,3)
      msq(1,4)=msq(1,4)+tempw(1,4)

      elseif (j.eq.2) then
      do k=1,nf
      do ll=k+1,nf
      msq(ll,k)=temp(k,ll)
      enddo
      enddo
      msq(3,2)=msq(3,2)+tempw(2,3)
      msq(4,1)=msq(4,1)+tempw(1,4)

      elseif (j.eq.3) then
      do k=-nf,-1
      msq(k,k)=temp(-k,-k)*stat
      do ll=k+1,-1
      msq(k,ll)=temp(-ll,-k)
      enddo
      enddo
      msq(-3,-2)=msq(-3,-2)+tempw(1,4)
      msq(-4,-1)=msq(-4,-1)+tempw(2,3)

      elseif (j.eq.4) then
      do k=-nf,-1
      do ll=k+1,-1
      msq(ll,k)=temp(-ll,-k)
      enddo
      enddo
      msq(-2,-3)=msq(-2,-3)+tempw(1,4)
      msq(-1,-4)=msq(-1,-4)+tempw(2,3)

c--- qbar-q
      elseif (j.eq.5) then
      do k=-nf,-1
      msq(k,-k)=temp(-k,-k)
      do ll=1,nf
      if (abs(k) .lt. abs(ll)) then
      msq(k,ll)=temp(-k,ll)
      endif
      enddo
      enddo
      msq(-1,3)=msq(-1,3)+tempw(2,3)
      msq(-2,4)=msq(-2,4)+tempw(1,4)
      
c--- qbar-q
      elseif (j.eq.6) then
      do k=-nf,-1
      do ll=1,nf
      if (abs(k) .gt. abs(ll)) then
      msq(k,ll)=temp(ll,-k)
      endif
      enddo
      enddo
      msq(-3,1)=msq(-3,1)+tempw(1,4)
      msq(-4,2)=msq(-4,2)+tempw(2,3)

c--- q-qbar
      elseif (j.eq.7) then
      do k=-nf,-1
      msq(-k,k)=temp(-k,-k)
      do ll=1,nf
      if (abs(k) .lt. abs(ll)) then
      msq(ll,k)=temp(-k,ll)
      endif
      enddo
      enddo
      msq(3,-1)=msq(3,-1)+tempw(2,3)
      msq(4,-2)=msq(4,-2)+tempw(1,4)

c--- q-qbar
      elseif (j.eq.8) then
      do k=-nf,-1
      do ll=-nf,-1
      if (abs(k) .lt. abs(ll)) then
      msq(-k,ll)=temp(-k,-ll)
      endif
      enddo
      enddo
      msq(1,-3)=msq(1,-3)+tempw(1,4)
      msq(2,-4)=msq(2,-4)+tempw(2,3)
      
c--- q-qbar extra pieces
      elseif (j.eq.9) then
      do k=1,nf
      do ll=1,nf
      if (k .lt. ll) then
      msq(k,-k)=msq(k,-k)+temp(k,ll)
      endif
      enddo
      enddo
      msq(1,-2)=msq(1,-2)+tempw(1,4) ! d u~ -> c~ s
      msq(3,-4)=msq(1,-2)
      msq(2,-1)=msq(2,-1)+tempw(2,3)
      msq(4,-3)=msq(2,-1)

c--- q-qbar extra pieces
      elseif (j.eq.10) then
      do k=1,nf
      do ll=1,nf
      if (k .gt. ll) then
      msq(k,-k)=msq(k,-k)+temp(ll,k)
      endif
      enddo
      enddo
 
c--- qbar-q extra pieces
      elseif (j.eq.11) then
      do k=1,nf
      do ll=1,nf
      if (k .lt. ll) then
      msq(-k,k)=msq(-k,k)+temp(k,ll)
      endif
      enddo
      enddo
      msq(-2,1)=msq(-2,1)+tempw(1,4) ! u~ d -> c~ s
      msq(-4,3)=msq(-2,1)
      msq(-1,2)=msq(-1,2)+tempw(2,3) ! d~ u -> s~ c
      msq(-3,4)=msq(-1,2)

c--- qbar-q extra pieces
      elseif (j.eq.12) then
      do k=1,nf
      do ll=1,nf
      if (k .gt. ll) then
      msq(-k,k)=msq(-k,k)+temp(ll,k)
      endif
      enddo
      enddo
  
      endif

      enddo

      return

   77 format(' *      W-mass^2     (',f11.5,',',f11.5,')      *')
   78 format(' *      Z-mass^2     (',f11.5,',',f11.5,')      *')
   79 format(' *  sin^2(theta_w)   (',f11.5,',',f11.5,')      *')

      end
      


      subroutine couplz(xw)
      implicit none       
      include 'constants.f'
      include 'zcouple.f'
      include 'ewcharge.f'
c---calculate the couplings as given in Kunszt and Gunion
c---Modified to notation of DKS (ie divided by 2*sw*cw)
c---xw=sin^2 theta_w
      integer j
      double precision xw
      sin2w=two*sqrt(xw*(1d0-xw))
      do j=1,nf
      l(j)=(tau(j)-two*Q(j)*xw)/sin2w
      r(j)=      (-two*Q(j)*xw)/sin2w
      enddo

      le=(-1d0-two*(-1d0)*xw)/sin2w
      re=(-two*(-1d0)*xw)/sin2w

      ln=(+1d0-two*(+0d0)*xw)/sin2w
      rn=0d0
             
      return
      end
