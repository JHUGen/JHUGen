      subroutine gg_dm_top(p,msq)
      implicit none
      include 'types.f'
      
c---Matrix element squared averaged over initial colors and spins
c     f(-p1) + f(-p2) --> H + f(p5)
c                         |
c                         --> b(p3)+bbar(p4)
c                            
c--all momenta incoming
c
c--- Matrix elements are taken from:
c--- R.~K.~Ellis, I.~Hinchliffe, M.~Soldate and J.~J.~van der Bij,
c--- %``Higgs Decay To Tau+ Tau-: A Possible Signature Of Intermediate
c--- % Mass Higgs Bosons At The SSC,''
c--- Nucl.\ Phys.\ B {\bf 297}, 221 (1988).
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'dm_params.f' 
      include 'sprods_com.f'
      include 'ewcouple.f' 
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),gg,qg,gq,qq,hdecay
      real(dp):: ehsvm3_dm,ehsvm4_dm,s34
!      real(dp):: mb_eff,mt_eff,massfrun,msqhtautau,msqgamgam
      real(dp):: propsq
      complex(dp):: cprop
      logical:: running_width 
      real(dp):: medwidth_orig
      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

C   Deal with Higgs decay
!      if (hdecaymode == 'tlta') then
!          call htautaudecay(p,3,4,hdecay)
!      elseif (hdecaymode == 'bqba') then
!          call hbbdecay(p,3,4,hdecay)
!      elseif (hdecaymode == 'gaga') then
!          hdecay=msqgamgam(hmass)
!      else
!      write(6,*) 'Unimplemented process in qqb_higgs'
!      stop
!      endif
!      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
!      origmbsq=mbsq
       mbsq=mt**2

!----- overall factor = v^2/dm_lam^6 in effective theory 
!======               = v/dm_lam * (1/prop) in full 
      propsq=1d0
!      gdm=sqrt(g_dmx*g_dmq)
!------ only coupling to DM is defined in full theory
      gdm=g_dmx

!------ store old width 
      running_width=.false.
      medwidth_orig=medwidth       
      if(running_width) then 
!--------- newwidth 
         medwidth=s34/medmass**2*medwidth_orig 
      endif
      call dmsdecay(p,3,4,hdecay) 
      if(effective_th) then 
         hdecay=hdecay*4d0*wmass**2/gwsq*one/(dm_lam**6) 
!         write(6,*) sqrt(4d0*wmass**2/gwsq),mbsq
      else
!         write(6,*) medmass,medwidth,gdm,g_dmx,g_dmq
         cprop=cone/cplx2((s34-medmass**2),medmass*medwidth)
         propsq=abs(cprop)**2 
!========= Note we remove the *v because of the full theory. 
         hdecay=hdecay*propsq*gdm**2
      endif

      
      gg=+avegg*ehsvm3_dm(s(1,2),s(1,5),s(2,5),s34)*hdecay
      qq=+aveqq*ehsvm4_dm(s(1,2),s(1,5),s(2,5),s34)*hdecay
      qg=-aveqg*ehsvm4_dm(s(1,5),s(1,2),s(2,5),s34)*hdecay
      gq=-aveqg*ehsvm4_dm(s(2,5),s(1,5),s(1,2),s34)*hdecay
  
      
      do j=-nf,nf    
      do k=-nf,nf
      msq(j,k)=0d0

      if ((j== 0) .or. (k==0)) then
           if ((j== 0) .and. (k==0)) then
                msq(j,k)=gg
           elseif ((j==0).and.(k.ne.0)) then
                msq(j,k)=gq
           elseif ((j.ne.0).and.(k==0)) then
                msq(j,k)=qg
           endif
      elseif ((j==-k).and. (j.ne.0)) then
           msq(j,k)=qq
      endif

      enddo
      enddo

!---- reset medwidth 
      medwidth=medwidth_orig
      return
      end
      
      
      function ehsvm3_dm(s,t,u,s34)
      implicit none
      include 'types.f'
      real(dp):: ehsvm3_dm
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
c---Matrix element squared Eqn 2.2 of EHSV
      complex(dp):: ehsva2_dm,ehsva4_dm
      real(dp):: s,t,u,s34
      logical:: approx
      parameter(approx=.false.)
c--- approx TRUE uses the heavy fermion approximation to Msq
      


      if (approx) then
      ehsvm3_dm=gwsq/pi*as**3*xn*V/9d0*(
     &        s34**4+s**4+t**4+u**4)/s/t/u/wmass**2
      else
      ehsvm3_dm=
     & abs(ehsva2_dm(s,t,u,s34))**2+abs(ehsva2_dm(u,s,t,s34))**2
     &        +abs(ehsva2_dm(t,u,s,s34))**2
     & +abs(ehsva4_dm(s,t,u,s34))**2 
      ehsvm3_dm=gwsq/pi*as**3*xn*V*s34**4/(s*t*u*wmass**2)*ehsvm3_dm
      endif
      
      return
      end
      
      
      function ehsvm4_dm(s,t,u,s34)
      implicit none
      include 'types.f'
      real(dp):: ehsvm4_dm
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
c---Matrix element squared Eqn 2.6 of EHSV
      complex(dp):: ehsva5_dm
      real(dp):: s,t,u,s34
      logical:: approx
      parameter(approx=.false.)
c--- approx TRUE uses the heavy fermion approximation to Msq


      if (approx) then
      ehsvm4_dm=gwsq/pi*as**3*V/18d0*(u**2+t**2)/s/wmass**2
      else
      ehsvm4_dm=abs(ehsva5_dm(s,t,u,s34))**2
      ehsvm4_dm=gwsq/(4d0*pi)*as**3*V/2d0*(u**2+t**2)/(s*wmass**2)
     & *s34**2/(u+t)**2*ehsvm4_dm
      endif
      
      return
      end
      
      
 
