      subroutine gen_realps(r,p,pswt,*)      
c--- calls appropriate subroutine to choose a 
c--- real-emission phase space point for the process
c--- specified by the variable 'case'   
c---
c---    input: r of random numbers, r
c---    output: momenta p and phase-space weight pswt
c---    note: common block 'npart' is also filled here
c---    
c---    alternate return if generated point should be discarded
      implicit none
      include 'constants.f'
      include 'breit.f'
      include 'dm_params.f'
      include 'frag.f'
      include 'ipsgen.f'
      include 'limits.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'npart.f'
      include 'phasemin.f'
      include 'process.f'
      include 'energy.f'
      logical vetow_2gam
      integer ii
      double precision r(mxdim),p(mxpart,4),pswt,
     & ptmp,m3,m4,m5,wt34,wt345,wt346,wt3456,wtprop,
     & s34,s345,s346,s3456,dot,wtips(4)
      
c--- statement function
      wtprop(s34,wmass,wwidth)=(s34-wmass**2)**2+(wmass*wwidth)**2
        
c--- processes that use "gen3"
      if     ( (case .eq. 'W_only')
     .    .or. (case .eq. 'Z_only')
     .    .or. (case .eq. 'ggfus0')
     .    .or. (case .eq. 'Higaga')
     .   .or.  (case .eq. 'Wcsbar')
     .   .or.  (case .eq. 'Wcs_ms') ) then
        npart=3
c        if (new_pspace) then
c          call gen3a(r,p,pswt,*999)
c        else
          call gen3(r,p,pswt,*999)
c        endif

c--- processes that use "gen3jet"     
      elseif (case .eq. 'gamgam') then
        npart=3
c        call gen3(r,p,pswt,*999)
        call gen3jetgaga(r,p,pswt,*999)

c--- processes that use "gen3m"     
      elseif ( (case .eq. 'tt_tot')
     .    .or. (case .eq. 'bb_tot')
     .    .or. (case .eq. 'cc_tot') ) then
        m3=mass2
        m4=mass2
        m5=0d0
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)

      elseif ((case .eq. 'twojet') .or. (case .eq. 'dirgam')
     &   .or. (case .eq. 'hflgam')) then
        npart=3
        if(frag) then
c---       this phase space does a better job for the photons
           call gen_photons_jets(r,1,2,p,pswt,*999)
        else
c---       this phase space does a better job for the jets 
           call gen3jet(r,p,pswt,*999)
        endif


c--- processes that use "gen4"     
      elseif ( (case .eq. 'W_cjet')
     .   .or.  (case .eq. 'Wbfrmc')
     .   .or.  (case .eq. 'W_tndk')
     .   .or.  (case .eq. 'epem3j')   ) then
        npart=4
        call gen4(r,p,pswt,*999)
                  
c--- processes that use "gen4"     
      elseif ((case .eq. 'qg_tbq') .or. (case .eq. 'qq_tbg')) then
        npart=4
        call gen4(r,p,pswt,*999)
                  
c--- processes that use "gen4mdk"     
      elseif (case .eq. '4ftwdk') then
        npart=6
        call gen4mdk(r,p,pswt,*999)

      elseif (case .eq. 'Hi_Zga') then
        npart=4
        call gen_HZgamj(r,p,pswt,*999)
                  
c--- processes that use "gen4mdkrad"     
      elseif (case .eq. 'dk_4ft') then
        npart=6
        call gen4mdkrad(r,p,pswt,*999)
                  
c--- processes that use "gen4mdk"     
      elseif ( (case .eq. 'Z_tdkj')
     &     .or.(case .eq. 'H_tdkj')) then
        npart=7
        call gen5mdk(r,p,pswt,*999)
      
c--- processes that use "gen4_3M"
      elseif ( (case .eq. 'tottth') ) then
        m3=mt
        m4=mt
        m5=hmass
        npart=4
        call gen4_3m(r,p,m3,m4,m5,pswt,*999)        
          
c--- processes that use "gen5"     
      elseif ( (case .eq. 'Wbbmas') 
     . .or. (case .eq. 'Wttmas')
     . .or. (case .eq. 'WWqqdk')) then
        npart=5
        call gen5(r,p,pswt,*999)
      elseif ((case .eq. 'Z_tjet') .or. (case .eq. 'H_tjet')) then
        npart=5
        taumin=(mt/sqrts)**2
        call gen5(r,p,pswt,*999)

      elseif (case .eq. 'HWWdkW') then
        npart=5
        call gen5h(r,p,pswt,*999)
                  
c--- processes that use "gen6"     
      elseif ( 
     .      (case .eq. 'W_twdk') 
     . .or. (case .eq. 'Wtdkay')
     . .or. (case .eq. 'HWWjet')
     . .or. (case .eq. 'HZZjet')
     . ) then
        npart=6
        call gen6(r,p,pswt,*999)
                  
c--- processes that use "gen7"     
      elseif ( 
     .      (case .eq. 'qq_HWW')
     . .or. (case .eq. 'qq_HZZ')
     . .or. (case .eq. 'WH__WW')
     . .or. (case .eq. 'WH__ZZ')
     . .or. (case .eq. 'ZH__WW')
     . .or. (case .eq. 'ZH__ZZ')
     . .or. (case .eq. 'HWW2jt')
     . .or. (case .eq. 'HZZ2jt')
     . .or. (case .eq. 'tt_ldk')
     . .or. (case .eq. 'tt_hdk')
     . .or. (case .eq. 'tt_udk')
     . .or. (case .eq. 'tthWdk')
     . ) then
        npart=7
        call gen7(r,p,pswt,*999)
                  
c--- processes that use "gen7m"     
      elseif ( (case .eq. 'tt_bbl') 
     .    .or. (case .eq. 'tt_bbh')
     .    .or. (case .eq. 'tt_bbu')
     . ) then
        m3=mt
        m4=mt
        m5=0d0
        npart=7
      call gen7m(r,p,m3,m4,m5,pswt,*999)

c--- processes that use "gen9"     
      elseif ( (case .eq. 'qq_ttw')) then 
        npart=9
      call gen9_rap(r,p,pswt,*999)

      elseif ( (case .eq. 'ttwldk')) then 
        npart=9
      call gen9dk_rap(r,p,pswt,*999)

c--- processes that use "gen_njets" with an argument of "2"     
      elseif ( (case .eq. 'W_1jet')
     .    .or. (case .eq. 'Wcjet0')
     .    .or. (case .eq. 'Z_1jet')
     .    .or. (case .eq. 'H_1jet')
     .    .or. (case .eq. 'ggfus1')
     .    .or. (case .eq. 'Hgagaj')
     .    .or. (case .eq. 'gQ__ZQ') ) then
        npart=4
c        if (new_pspace) then
c          call gen4a(r,p,pswt,*999)      
c        else
          call gen_njets(r,2,p,pswt,*999)
c        endif 
        
c--- processes that use "gen_njets" with an argument of "3"
      elseif ( (case .eq. 'Wbbbar')
     .    .or. (case .eq. 'W_2jet')
     .    .or. (case .eq. 'Z_2jet')
     .    .or. (case .eq. 'Zbbbar')
     .    .or. (case .eq. 'W_bjet')
     .    .or. (case .eq. 'Z_bjet')
     .    .or. (case .eq. 'qq_Hqq')
     .    .or. (case .eq. 'qq_Hgg')
     .    .or. (case .eq. 'ggfus2')
     .    .or. (case .eq. 'gagajj')) then
        npart=5
        call gen_njets(r,3,p,pswt,*999) 
        
c--- processes that use "gen_Vphotons_jets"     
      elseif ( (case .eq. 'Wgamma')
     .   .or.  (case .eq. 'Zgamma')   ) then
        npart=4
c      if (new_pspace) then
c          call gen_vgamj(r,p,pswt,*999)
c        else
          call gen_Vphotons_jets(r,1,1,p,pswt,*999)
c        endif
                
c--- processes that use "gen_photons_jets"     
      elseif (case .eq. 'gmgmjt') then
        npart=4
        call gen_photons_jets(r,2,2,p,pswt,*999)

      elseif (case .eq. 'trigam') then
        npart=4
        call gen_photons_jets(r,3,1,p,pswt,*999)

c--- special treatment for Z+gamma+gamma 
      elseif ((case .eq. 'W_2gam') .or. (case .eq. 'Z_2gam')) then
        npart=5
        if  (ipsgen .eq. 1) then
            call gen_Vphotons_jets(r,2,1,p,pswt,*999) !AA+AB
        elseif  (ipsgen .eq. 2) then
            call gen_Vphotons_jets_dkrad2(r,2,1,p,pswt,*999) !BB+BC
        elseif  (ipsgen .eq. 3) then
           call gen_Vphotons_jets_dkrad(r,2,1,p,pswt,*999) !CC+AC+CD
        elseif  (ipsgen .eq. 4) then
           call gen_Vphotons_jets_dkrad(r,2,1,p,pswt,*999) !DD+AD+BD
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
        if (case .eq. 'W_2gam') then
c          if (vetow_2gam(p)) goto 999 ! partition PS according to ipsgen
          s34=2d0*dot(p,3,4)
          s345=s34+2d0*dot(p,3,5)+2d0*dot(p,4,5)
          s346=s34+2d0*dot(p,3,6)+2d0*dot(p,4,6)
          s3456=s345+s346-s34+2d0*dot(p,5,6)
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
      elseif((case.eq.'dm_jet').or.(case.eq.'dm_gam')) then 
         m3=xmass
         m4=xmass
         npart=4
         call gen4m(r,p,m3,m4,0d0,0d0,pswt,*999)
c--- special treatment for Z+gamma+jet
      elseif (case .eq. 'Zgajet') then
        npart=5
        if     (ipsgen .eq. 1) then
          call gen_Vphotons_jets(r,1,2,p,pswt,*999)
        elseif (ipsgen .eq. 2) then
          call gen_Vphotons_jets_dkrad(r,1,2,p,pswt,*999)
      else
        write(6,*) 'Parameter ipsgen should be 1 or 2'
        write(6,*) 'ipsgen = ',ipsgen
        stop
        endif
                  
c--- processes that use "gen_stop" with an argument of "1"
      elseif ( (case .eq. 'ttdkay')
     .    .or. (case .eq. 'tdecay') ) then
        npart=5
        call gen_stop(r,1,p,pswt,*999)

c--- processes that use "gen_stop" with an argument of "2"
      elseif (case .eq. 'bq_tpq') then
        npart=5
        call gen_stop(r,2,p,pswt,*999)
        
c--- processes that use "gen_stop" with an argument of "2"
      elseif (case .eq. 't_bbar') then
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
      
