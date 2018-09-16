      subroutine gen_lops(r,p,pswt,*)      
c--- calls appropriate subroutine to choose a 
c--- leading-order phase space point for the process
c--- specified by the variable 'case'   
c---
c---    input: vector of random numbers, r
c---    output: momenta p and phase-space weight pswt
c---    note: common block 'npart' is also filled here
c---    
c---    alternate return if generated point should be discarded
      implicit none
      include 'constants.f'
      include 'breit.f'
      include 'dm_params.f'
      include 'ipsgen.f'
      include 'limits.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'npart.f'
      include 'phasemin.f'
      include 'process.f'
      include 'zerowidth.f'
      include 'energy.f'
      logical vetow_2gam
      integer ii
      double precision r(mxdim),p(mxpart,4),pswt,
     & ptmp,m3,m4,m5,wt34,wt345,wt346,wt3456,wtprop,
     & s34,s345,s346,s3456,dot,wtips(4)
      
c--- statement function
      wtprop(s34,wmass,wwidth)=(s34-wmass**2)**2+(wmass*wwidth)**2

c--- processes that use "gen2"
      if     ( (case .eq. 'W_only')
     .    .or. (case .eq. 'Z_only')
     .    .or. (case .eq. 'Higaga')
     .    .or. (case .eq. 'Wcsbar')
     .    .or. (case .eq. 'Wcs_ms')
     .    .or. (case .eq. 'vlchk2') ) then
        if (case .eq. 'vlchk2') then
          wsqmin=0d0
          wsqmax=sqrts**2
        endif
        npart=2
c        if (new_pspace) then
c          call gen2a(r,p,pswt,*999)
c        else
          call gen2(r,p,pswt,*999)
c        endif

c--- processes that use "gen2jet"     
      elseif ((case .eq. 'twojet') 
     .   .or. (case .eq. 'dirgam')
     .   .or. (case .eq. 'hflgam')
     .   .or. (case .eq. 'gamgam')) then
        npart=2
        call gen2jet(r,p,pswt,*999)

c--- processes that use "gen2m"     
      elseif ( (case .eq. 'ggfus0')
     .    .or. (case .eq. 'tt_tot')
     .    .or. (case .eq. 'bb_tot')
     .    .or. (case .eq. 'cc_tot') ) then
        npart=2
        call gen2m(r,p,pswt,*999)
          
c--- processes that use "gen3"     
      elseif ( (case .eq. 'W_cjet') 
     .   .or.  (case .eq. 'Wbfrmc')
     .   .or.  (case .eq. 'W_tndk')
     .   .or.  (case .eq. 'vlchwn')
     .   .or.  (case .eq. 'epem3j') ) then
        npart=3
        call gen3(r,p,pswt,*999)

      elseif((case.eq.'dm_jet').or.(case.eq.'dm_gam')) then 
         m3=xmass
         m4=m3
         m5=0d0
         npart=3     
         call gen3m(r,p,m3,m4,m5,pswt,*999)

c--- processes that use "gen3h"     
      elseif (case .eq. 'Hi_Zga') then
        npart=3
        call gen3h(r,p,pswt,*999)

c---  processes that use "gen3jet"     
      elseif ( (case .eq. 'Wgamma') 
     .   .or.  (case .eq. 'Zgamma')
     .   .or.  (case .eq. 'W_frag') 
     .   .or.  (case .eq. 'Z_frag')
     .   .or.  (case .eq. 'vlchk3') ) then
        if (case .eq. 'vlchk3') then
          wsqmin=0d0
          wsqmax=sqrts**2
        endif
        npart=3
        call gen_Vphotons_jets(r,1,0,p,pswt,*999)
        
c--- processes that use "gen3jetgaga"     
      elseif (case .eq. 'gmgmjt') then
        npart=3
c        call gen3jetgaga(r,p,pswt,*999)
c        call gen_photons_jets(r,2,1,p,pswt,*999)
        call gen3(r,p,pswt,*999)

      elseif (case .eq. 'gmgmjj') then
        npart=4
c        call gen3jetgaga(r,p,pswt,*999)
        call gen_photons_jets(r,2,2,p,pswt,*999)

c--- processes that use "gen3jetgaga"     
      elseif (case .eq. 'trigam') then
        npart=3
c        call gen3jetgaga(r,p,pswt,*999)
c        call gen_photons_jets(r,3,0,p,pswt,*999)
        call gen3(r,p,pswt,*999)

      elseif(case .eq.'fourga') then 
        npart=4
        call gen_photons_jets(r,4,0,p,pswt,*999)

c--- processes that use "gen_Vphotons_jets"     
c--- special treatment for W+gamma+jet and Z+gamma+jet 
      elseif ( (case .eq. 'Wgajet') .or. (case .eq. 'Zgajet') ) then
        npart=4
c      if (new_pspace) then
c          call gen_vgamj(r,p,pswt,*999) ! New PS routine
c      else
        if     (ipsgen .eq. 1) then
          call gen_Vphotons_jets(r,1,1,p,pswt,*999)
        elseif (ipsgen .eq. 2) then
          call gen_Vphotons_jets_dkrad(r,1,1,p,pswt,*999)
        else
        write(6,*) 'Parameter ipsgen should be 1 or 2'
        write(6,*) 'ipsgen = ',ipsgen
        stop
        endif
c      endif
       
c--- special treatment for Z+gamma+gamma 
      elseif ( (case .eq. 'Z_2gam') .or. (case .eq. 'W_2gam')) then
        npart=4
        if  (ipsgen .eq. 1) then
            call gen_Vphotons_jets(r,2,0,p,pswt,*999) !AA+AB
        elseif  (ipsgen .eq. 2) then
            call gen_Vphotons_jets_dkrad2(r,2,0,p,pswt,*999) !BB+BC
        elseif  (ipsgen .eq. 3) then
           call gen_Vphotons_jets_dkrad(r,2,0,p,pswt,*999) !CC+AC+CD
        elseif  (ipsgen .eq. 4) then
           call gen_Vphotons_jets_dkrad(r,2,0,p,pswt,*999) !DD+AD+BD
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

c--- special treatment for Z+gamma+gamma+jet 
      elseif ( (case .eq. 'Z2gajt') ) then
        npart=5
        if (ipsgen .eq. 1) then
           call gen_Vphotons_jets(r,2,1,p,pswt,*999)  !AA+AB
        elseif (ipsgen .eq. 2) then
           call gen_Vphotons_jets_dkrad2(r,2,1,p,pswt,*999)  !BB+BC
        elseif (ipsgen .eq. 3) then
           call gen_Vphotons_jets_dkrad(r,2,1,p,pswt,*999)  !CC+AC+CD 
        elseif (ipsgen .eq. 4) then
           call gen_Vphotons_jets_dkrad(r,2,1,p,pswt,*999)  !DD+AD+BD
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

c--- special treatment for Z+gamma+jet+jet 
      elseif ( (case .eq. 'Zga2jt') ) then
        npart=5
        if (ipsgen .eq. 1) then
          call gen_Vphotons_jets(r,1,2,p,pswt,*999)
        elseif (ipsgen .eq. 2) then
          call gen_Vphotons_jets_dkrad(r,1,2,p,pswt,*999)
        else
          write(6,*) 'Parameter ipsgen should be 1 or 2'
          write(6,*) 'ipsgen = ',ipsgen
          stop
        endif

      elseif ((case .eq. 'thrjet') .or. (case .eq. 'gamjet')) then
!        call gen3jet(r,p,pswt,*999)
        call gen_photons_jets(r,1,2,p,pswt,*999)
        npart=3

c--- processes that use "gen3m"
      elseif ( (case .eq. 'tottth') ) then
        m3=mt
        m4=mt
        m5=hmass
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)
          
      elseif ( (case .eq. 'totttz') ) then
        m3=mt
        m4=mt
        m5=zmass
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)
        
c--- processes that use "gen3m"
      elseif (case .eq. 'tt_glu') then
        m3=mass2
        m4=mass2
        m5=0d0
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)

c--- processes that use "gen3m"
      elseif ((case .eq. 'qg_tbq') .or. (case .eq. 'qq_tbg')) then
        m3=mt
        m4=mb
        m5=0d0
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)
      
c--- processes that use gen3mdk
      elseif ( (case .eq. '4ftwdk') .or. (case .eq. 'dk_4ft') ) then
        m3=mt
        m4=mb
        m5=0d0
        npart=5
        call gen3mdk(r,p,m3,m4,m5,pswt,*999)
        
c--- processes that use "gen3m_rap"     
      elseif ( (case .eq. 'vlchm3') ) then
        taumin=(2d0*mt/sqrts)**2
        m3=mt
        m4=mt
        npart=3
        call gen3m_rap(r,p,m3,m4,pswt,*999)
          
c--- processes that use "gen4"     
      elseif (case .eq. 'qgtbqq') then
        npart=4
        call gen4(r,p,pswt,*999)
                  
      elseif (case .eq. 'H_tjet') then
        npart=4
        taumin=((mt+hmass)/sqrts)**2
        call gen4(r,p,pswt,*999)

      elseif (case .eq. 'Z_tjet') then
        npart=4
        if (zerowidth) then
          taumin=((mt+zmass)/sqrts)**2
        else
          taumin=((mt+dsqrt(wsqmin))/sqrts)**2
        endif
        call gen4(r,p,pswt,*999)
        
      elseif((case.eq.'dm_gaj').or.(case.eq.'dm2jet')) then 
         m3=xmass
         m4=xmass
         npart=4
         call gen4m(r,p,m3,m4,0d0,0d0,pswt,*999)

c--- processes that use "gen4handc"     
      elseif ( (case .eq. 'HZZint')
     .    .or. (case .eq. 'HZZH+i')
     .    .or. (case .eq. 'ggZZ4l')
     .    .or. (case .eq. 'HWWint')
     .    .or. (case .eq. 'HWWH+i')
     .    .or. (case .eq. 'ggWW4l') ) then
        npart=4
c        call gen4_intf(r,p,pswt,*999)
        call gen4handc(r,p,pswt,*999)
                  
c--- processes that use "gen4h"     
      elseif ( (case .eq. 'HWW_4l')
     .    .or. (case .eq. 'HWWdkW')
     .    .or. (case .eq. 'HWW2lq')
     .    .or. (case .eq. 'HWW_tb')
     .    .or. (case .eq. 'HZZ_4l')
     .    .or. (case .eq. 'HmZZ4l')
     .    .or. (case .eq. 'HZZ_tb') ) then
        npart=4
        call gen4h(r,p,pswt,*999)
         
c--- processes that use "gen4hvv"     
      elseif ( (case .eq. 'HVV_tb') ) then
        npart=4
        call gen4hvv(r,p,pswt,*999)
         
c--- processes that use "gen4handcvv"     
      elseif ( (case .eq. 'ggVV4l') ) then
        npart=4
        call gen4handcvv(r,p,pswt,*999)
         
c--- processes that use "gen4vv"     
      elseif ( (case .eq. 'ggVVbx') ) then
        npart=4
        call gen4vv(r,p,pswt,*999)
         
c--- processes that use "gen4mdk"     
      elseif ((case .eq. '4ftjet') 
     .   .or. (case .eq. 'Z_tdkj')
     .   .or. (case .eq. 'H_tdkj')) then
        npart=6
        call gen4mdk(r,p,pswt,*999)
      
c--- processes that use "gen5mdk"     
      elseif  (case .eq. 'Ztdk2j') then
        npart=7
        call gen5mdk(r,p,pswt,*999)
      
c--- processes that use "gen5" 
      elseif ( (case .eq. 'W_twdk') 
     . .or.    (case .eq. 'HWWjet') 
     . .or.    (case .eq. 'HZZjet') 
     . .or.    (case .eq. 'Wtdkay') 
     . .or.    (case .eq. 'WW_jet') 
     . .or.    (case .eq. 'ZZ_jet') 
     . .or.    (case .eq. 'Zt2jet') 
     . .or.    (case .eq. 'Wbbjem') 
     . .or.    (case .eq. 'vlchwt') 
     . ) then 
        npart=5
        call gen5(r,p,pswt,*999)

c--- processes that use "gen5" 
      elseif (case .eq. 'HZZqgI')  then
        npart=5
        call gen5(r,p,pswt,*999)
      
c--- processes that use "gen6"     
      elseif ( (case .eq. 'tt_bbl')
     .    .or. (case .eq. 'tt_bbh')
     .    .or. (case .eq. 'tt_bbu')
     .    .or. (case .eq. 'tt_ldk')
     .    .or. (case .eq. 'tt_hdk')
     .    .or. (case .eq. 'tt_udk')
     .    .or. (case .eq. 'tthWdk')
     .    .or. (case .eq. 'Wtbwdk')
     .    .or. (case .eq. 'ttZbbl')
     .    .or. (case .eq. 'tautau')
     .    .or. (case .eq. 'vlchk6')
     .    .or. (case .eq. 'vlchwg')
     .    .or. (case .eq. 'vlchwh')
     .    .or. (case .eq. 'qq_HWW')
     .    .or. (case .eq. 'qq_HZZ')
     .    .or. (case .eq. 'qqWWqq')
     .    .or. (case .eq. 'qqWWss')
     .    .or. (case .eq. 'qqWZqq')
     .    .or. (case .eq. 'qqZZqq')
     .    .or. (case .eq. 'qqVVqq')
     .    .or. (case .eq. 'HWW2jt')
     .    .or. (case .eq. 'HZZ2jt')
     .    .or. (case .eq. 'WH__ZZ')
     .    .or. (case .eq. 'WH__WW')  
     .    .or. (case .eq. 'ZH__WW')
     .    .or. (case .eq. 'ZH__ZZ')
     .    .or. (case .eq. 'WpWp2j')
     .    .or. (case .eq. 'WpmZjj')
     .    .or. (case .eq. 'WpmZbj')
     .    .or. (case .eq. 'WpmZbb')
     .    .or. (case .eq. 'WW2jet')) then
        npart=6
        call gen6(r,p,pswt,*999)

c--- processes that use "gen7"     
      elseif ( (case .eq. 'qq_ttg') ) then
        m3=mt
        m4=mt
        m5=0d0
        npart=7
        call gen7m(r,p,m3,m4,m5,pswt,*999)

      elseif ( (case.eq.'WpWp3j')
     &     .or.(case.eq.'HWW3jt')
     &     .or.(case.eq.'HZZ3jt') )  then
        npart=7
      call gen7(r,p,pswt,*999)

c--- processes that use "gen8"     
      elseif ( (case .eq. 'qq_tth') 
     .    .or. (case .eq. 'vlchk8') ) then
        npart=8
        call gen8(r,p,pswt,*999)
  
c--- processes that use "gen8"     
      elseif ( (case .eq. 'qq_ttz') 
     .    .or. (case .eq. 'qqtthz') ) then
        npart=8
c        call gen8(r,p,pswt,*999)
        call genttvdk(r,p,pswt,*999)
  
      elseif ( (case .eq. 'qq_ttw')  .or. (case .eq. 'ttwldk') ) then
        npart=8
        taumin=(2d0*mt/sqrts)**2
        call gen8_rap(r,p,pswt,*999)
          
      elseif (case .eq. 'tth_ww') then 
        npart=10
        call gen10(r,p,pswt,*999)
c--- processes that use "gen_njets" with an argument of "1"     
      elseif ( (case .eq. 'W_1jet')
     .    .or. (case .eq. 'Wcjet0')
     .    .or. (case .eq. 'Z_1jet')
     .    .or. (case .eq. 'H_1jet')
     .    .or. (case .eq. 'httjet')
     .    .or. (case .eq. 'attjet')
     .    .or. (case .eq. 'ggfus1')
     .    .or. (case .eq. 'Hgagaj')
     .    .or. (case .eq. 'gQ__ZQ') ) then
        npart=3
        call gen_njets(r,1,p,pswt,*999)
         
c--- processes that use "gen_njets" with an argument of "2"
      elseif ( (case .eq. 'Wbbbar')
     .    .or. (case .eq. 'W_2jet')
     .    .or. (case .eq. 'Z_2jet')
     .    .or. (case .eq. 'Zbbbar')
     .    .or. (case .eq. 'qq_Hqq')
     .    .or. (case .eq. 'qq_Hgg')
     .    .or. (case .eq. 'ggfus2')
     .    .or. (case .eq. 'gagajj')
     .    .or. (case .eq. 'W_bjet')
     .    .or. (case .eq. 'Wcjetg')
     .    .or. (case .eq. 'Z_bjet') ) then
        npart=4
        call gen_njets(r,2,p,pswt,*999)         
        
c--- processes that use "gen_njets" with an argument of "3"
      elseif ( (case .eq. 'W_3jet') 
     .    .or. (case .eq. 'Wbbjet') 
     .    .or. (case .eq. 'Z_3jet') 
     .    .or. (case .eq. 'Zbbjet') 
     .    .or. (case .eq. 'Wb2jet') 
     .    .or. (case .eq. 'qqHqqg')
     .    .or. (case .eq. 'ggfus3')
     .    .or. (case .eq. 'Zbjetg')
     .    .or. (case .eq. 'vlchk5') ) then
        npart=5
        call gen_njets(r,3,p,pswt,*999)      
c--- processes that use "gen_stop" with an argument of "1" (number of extra jets)
      elseif ( (case .eq. 'bq_tpq')
     .    .or. (case .eq. 'ttdkay')
     .    .or. (case .eq. 't_bbar')
     .    .or. (case .eq. 'tdecay') ) then
        npart=4
        call gen_stop(r,1,p,pswt,*999)
      
c--- DEFAULT: processes that use "gen4"
      else
        if ((case .eq. 'vlchk4') .or. (case .eq. 'vlchkm')) then
          wsqmin=0d0
          wsqmax=sqrts**2
        endif
        npart=4
        call gen4(r,p,pswt,*999)
      endif

c--- ensure no entry beyond npart+2
      p(npart+3,:)=0d0
      
      return

c--- alternate return
  999 continue   
      return 1
      
      end
      
