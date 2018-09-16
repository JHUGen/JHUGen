      subroutine gen_realps(r,p,pswt,*)      
      implicit none
      include 'types.f'
c--- calls appropriate subroutine to choose a 
c--- real-emission phase space point for the process
c--- specified by the variable 'case'   
c---
c---    input: r of random numbers, r
c---    output: momenta p and phase-space weight pswt
c---    note: common block 'npart' is also filled here
c---    
c---    alternate return if generated point should be discarded
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'breit.f'
      include 'dm_params.f'
      include 'frag.f'
      include 'ipsgen.f'
      include 'limits.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'npart.f'
      include 'phasemin.f'
      include 'kprocess.f'
      include 'energy.f'
      include 'taucut.f'
      logical:: vetow_2gam
      integer:: ii
      real(dp):: r(mxdim),p(mxpart,4),pswt,
     & ptmp,m3,m4,m5,wt34,wt345,wt346,wt3456,wtprop,
     & s34,s345,s346,s3456,dot,wtips(4)
      
c--- statement function
      wtprop(s34,wmass,wwidth)=(s34-wmass**2)**2+(wmass*wwidth)**2
        
c--- processes that use "gen3"
      if     ( (kcase==kW_only)
     &    .or. (kcase==kZ_only)
     &    .or. (kcase==kggfus0)
     &    .or. (kcase==kHigaga)
     &   .or.  (kcase==kWcsbar)
     &   .or.  (kcase==kWcs_ms) ) then
        npart=3
c        if (new_pspace) then
c          call gen3a(r,p,pswt,*999)
c        else
          call gen3(r,p,pswt,*999)
c        endif

c--- processes that use "gen3jet"     
      elseif ((kcase==kgamgam) .or. (kcase==kgg2gam)) then
        npart=3
c        call gen3(r,p,pswt,*999)
        call gen3jetgaga(r,p,pswt,*999)

c--- processes that use "gen3m"     
      elseif ( (kcase==ktt_tot)
     &    .or. (kcase==kbb_tot)
     &    .or. (kcase==kcc_tot) ) then
        m3=mass2
        m4=mass2
        m5=0._dp
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)

      elseif ((kcase==ktwojet) .or. (kcase==kdirgam)
     &   .or. (kcase==khflgam)) then
        npart=3
        if(frag) then
c---       this phase space does a better job for the photons
           call gen_photons_jets(r,1,2,p,pswt,*999)
        else
c---       this phase space does a better job for the jets 
           call gen3jet(r,p,pswt,*999)
        endif


c--- processes that use "gen4"     
      elseif ( (kcase==kW_cjet)
     &   .or.  (kcase==kWbfrmc)
     &   .or.  (kcase==kW_tndk)
     &   .or.  (kcase==kepem3j)   ) then
        npart=4
        call gen4(r,p,pswt,*999)
                  
c--- processes that use "gen4"     
      elseif ((kcase==kqg_tbq) .or. (kcase==kqq_tbg)) then
        npart=4
        call gen4(r,p,pswt,*999)
                  
c--- processes that use "gen4mdk"     
      elseif (kcase==k4ftwdk) then
        npart=6
        call gen4mdk(r,p,pswt,*999)

      elseif (kcase==kHi_Zga) then
        npart=4
        call gen_HZgamj(r,p,pswt,*999)
                  
c--- processes that use "gen4mdkrad"     
      elseif (kcase==kdk_4ft) then
        npart=6
        call gen4mdkrad(r,p,pswt,*999)
                  
c--- processes that use "gen4mdk"     
      elseif ( (kcase==kZ_tdkj)
     &     .or.(kcase==kH_tdkj)) then
        npart=7
        call gen5mdk(r,p,pswt,*999)
      
c--- processes that use "gen4_3M"
      elseif ( (kcase==ktottth) ) then
        m3=mt
        m4=mt
        m5=hmass
        npart=4
        call gen4_3m(r,p,m3,m4,m5,pswt,*999)        
          
c--- processes that use "gen5"     
      elseif ( (kcase==kWbbmas) 
     & .or. (kcase==kWttmas)
     & .or. (kcase==kWWqqdk)
     & .or. (kcase==kZHbbdk)
     & .or. (kcase==kWHbbdk)) then
        npart=5
        call gen5(r,p,pswt,*999)
      elseif ((kcase==kZ_tjet) .or. (kcase==kH_tjet)) then
        npart=5
        taumin=(mt/sqrts)**2
        call gen5(r,p,pswt,*999)

      elseif (kcase==kHWWdkW) then
        npart=5
        call gen5h(r,p,pswt,*999)
                  
c--- processes that use "gen6"
      elseif ( 
     &      (kcase==kW_twdk) 
     & .or. (kcase==kWtdkay)
     & .or. (kcase==kHWWjet)
     & .or. (kcase==kHZZjet)
     & .or. (kcase==kWH1jet)
     & .or. (kcase==kZH1jet)
     & ) then
        npart=6
        if ((usescet) .and. (abovecut)) then
          call genVHjjtaucut(r,p,pswt,*999)
c          call gen6(r,p,pswt,*999)
        else
          call gen6(r,p,pswt,*999)
        endif
                  
c--- processes that use "gen7"
      elseif ( 
     &      (kcase==kqq_HWW)
     & .or. (kcase==kqq_HZZ)
     & .or. (kcase==kWH__WW)
     & .or. (kcase==kWH__ZZ)
     & .or. (kcase==kZH__WW)
     & .or. (kcase==kZH__ZZ)
     & .or. (kcase==kHWW2jt)
     & .or. (kcase==kHZZ2jt)
     & .or. (kcase==ktt_ldk)
     & .or. (kcase==ktt_hdk)
     & .or. (kcase==ktt_udk)
     & .or. (kcase==ktthWdk)
     & ) then
        npart=7
        call gen7(r,p,pswt,*999)
                  
c--- processes that use "gen7m"     
      elseif ( (kcase==ktt_bbl) 
     &    .or. (kcase==ktt_bbh)
     &    .or. (kcase==ktt_bbu)
     & ) then
        m3=mt
        m4=mt
        m5=0._dp
        npart=7
      call gen7m(r,p,m3,m4,m5,pswt,*999)

c--- processes that use "gen9"     
      elseif ( (kcase==kqq_ttw)) then 
        npart=9
      call gen9_rap(r,p,pswt,*999)

      elseif ( (kcase==kttwldk)) then 
        npart=9
      call gen9dk_rap(r,p,pswt,*999)

c--- processes that use "gen_njets" with an argument of "2"
      elseif ( (kcase==kW_1jet)
     &    .or. (kcase==kWcjet0)
     &    .or. (kcase==kZ_1jet)
     &    .or. (kcase==kH_1jet)
     &    .or. (kcase==kggfus1)
     &    .or. (kcase==kHgagaj)
     &    .or. (kcase==kgQ__ZQ) ) then
        npart=4
c        if (new_pspace) then
c          call gen4a(r,p,pswt,*999)
c        else
          if ((usescet) .and. (abovecut)) then
            call gen4taucut(r,p,pswt,*999)
c            call gen_njets(r,2,p,pswt,*999)
          else
            call gen_njets(r,2,p,pswt,*999)
          endif
c        endif 
        
c--- processes that use "gen_njets" with an argument of "3"
      elseif ( (kcase==kWbbbar)
     &    .or. (kcase==kW_2jet)
     &    .or. (kcase==kZ_2jet)
     &    .or. (kcase==kZbbbar)
     &    .or. (kcase==kW_bjet)
     &    .or. (kcase==kZ_bjet)
     &    .or. (kcase==kqq_Hqq)
     &    .or. (kcase==kqq_Hgg)
     &    .or. (kcase==kggfus2)
     &    .or. (kcase==kgagajj)) then
        npart=5
        call gen_njets(r,3,p,pswt,*999) 
        
c--- processes that use "gen_Vphotons_jets"     
      elseif ( (kcase==kWgamma)
     &   .or.  (kcase==kZgamma)   ) then
        npart=4
c      if (new_pspace) then
c          call gen_vgamj(r,p,pswt,*999)
c        else
          call gen_Vphotons_jets(r,1,1,p,pswt,*999)
c        endif
                
c--- processes that use "gen_photons_jets"
      elseif (kcase==kgmgmjt) then
        npart=4
        if ((usescet) .and. (abovecut)) then
          call gen4taucut(r,p,pswt,*999)
        else
          call gen_photons_jets(r,2,2,p,pswt,*999)
        endif
        
      elseif (kcase==ktrigam) then
        npart=4
        call gen_photons_jets(r,3,1,p,pswt,*999)

c--- special treatment for Z+gamma+gamma 
      elseif ((kcase==kW_2gam) .or. (kcase==kZ_2gam)) then
        npart=5
        if  (ipsgen == 1) then
            call gen_Vphotons_jets(r,2,1,p,pswt,*999) !AA+AB
        elseif  (ipsgen == 2) then
            call gen_Vphotons_jets_dkrad2(r,2,1,p,pswt,*999) !BB+BC
        elseif  (ipsgen == 3) then
           call gen_Vphotons_jets_dkrad(r,2,1,p,pswt,*999) !CC+AC+CD
        elseif  (ipsgen == 4) then
           call gen_Vphotons_jets_dkrad(r,2,1,p,pswt,*999) !D.e+_dpA.e+_dpBD
           do ii=1,4
              ptmp=p(5,ii)
              p(5,ii)=p(6,ii)
              p(6,ii)=ptmp
           enddo
        else
           write(6,*) 'Parameter ipsgen should be 1 or 2 or 3 or 4'
           write(6,*) 'ipsgen = ',ipsgen
           stop
        endif
        if (kcase==kW_2gam) then
c          if (vetow_2gam(p)) goto 999 ! partition PS according to ipsgen
          s34=2._dp*dot(p,3,4)
          s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
          s346=s34+2._dp*dot(p,3,6)+2._dp*dot(p,4,6)
          s3456=s345+s346-s34+2._dp*dot(p,5,6)
          wt34=wtprop(s34,wmass,wwidth)
          wt345=wtprop(s345,wmass,wwidth)
          wt346=wtprop(s346,wmass,wwidth)
          wt3456=wtprop(s3456,wmass,wwidth)
          wtips(1)=wt345*wt346*wt3456
          wtips(2)=wt34*wt345*wt346
          wtips(3)=wt34*wt346*wt3456
          wtips(4)=wt34*wt345*wt3456
          pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(2)+wtips(3)+wtips(4))
        endif
      elseif((kcase==kdm_jet).or.(kcase==kdm_gam)) then 
         m3=xmass
         m4=xmass
         npart=4
         call gen4m(r,p,m3,m4,0._dp,0._dp,pswt,*999)
c--- special treatment for Z+gamma+jet
      elseif (kcase==kZgajet) then
        npart=5
        if     (ipsgen == 1) then
          call gen_Vphotons_jets(r,1,2,p,pswt,*999)
        elseif (ipsgen == 2) then
          call gen_Vphotons_jets_dkrad(r,1,2,p,pswt,*999)
      else
        write(6,*) 'Parameter ipsgen should be 1 or 2'
        write(6,*) 'ipsgen = ',ipsgen
        stop
        endif
                  
c--- processes that use "gen_stop" with an argument of "1"
      elseif ( (kcase==kttdkay)
     &    .or. (kcase==ktdecay) ) then
        npart=5
        call gen_stop(r,1,p,pswt,*999)

c--- processes that use "gen_stop" with an argument of "2"
      elseif (kcase==kbq_tpq) then
        npart=5
        call gen_stop(r,2,p,pswt,*999)
        
c--- processes that use "gen_stop" with an argument of "2"
      elseif (kcase==kt_bbar) then
        npart=5
        call gen_stop(r,2,p,pswt,*999)
                
c--- DEFAULT: processes that use "gen5"
      else
        npart=5
c        if (new_pspace) then
c          call gen5a(r,p,pswt,*999)
c        else
          call gen5(r,p,pswt,*999)    
c        endif
      endif
            

      return

c--- alternate return      
  999 continue   
      return 1
      
      end
      
