      subroutine init_is_functions()
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'is_functions_com.f'
      include 'plabel.f'
      integer:: j

c--- hadrons      
      do j=1,mxpart      
      if ( (plabel(j) == 'pp') .or. (plabel(j) == 'pj')
     & .or.(plabel(j) == 'bq') .or. (plabel(j) == 'ba')
     & .or.(plabel(j) == 'qj') ) then
         ishadarray(j) = .true. 
      else
         ishadarray(j) = .false. 
      endif
      enddo



c---  hadrons (b-quark sepical cases
      do j=1,mxpart      
      if ( (plabel(j) == 'pp') .or. (plabel(j) == 'pj')
     & .or.(plabel(j) == 'bq') .or. (plabel(j) == 'ba')
     & .or.(plabel(j) == 'qb') .or. (plabel(j) == 'ab')
     & .or.(plabel(j) == 'qj') ) then
         ishadspecarray(j) = .true. 
      else
         ishadspecarray(j) = .false. 
      endif
      enddo

c--- photons
      do j=1,mxpart      
      if (plabel(j) == 'ga') then
         isphotarray(j) = .true. 
      else
         isphotarray(j) = .false. 
      endif
      enddo

c--- charged leptons
      do j=1,mxpart      
      if (     (plabel(j) == 'el') .or. (plabel(j) == 'ea')
     &    .or. (plabel(j) == 'ml') .or. (plabel(j) == 'ma')
     &    .or. (plabel(j) == 'tl') .or. (plabel(j) == 'ta')) then
         isleptarray(j) = .true. 
      else
         isleptarray(j) = .false. 
      endif
      enddo

c--- electrons
      do j=1,mxpart      
      if (     (plabel(j) == 'el') .or. (plabel(j) == 'ea')) then
         iselectronarray(j) = .true. 
      else
         iselectronarray(j) = .false. 
      endif
      enddo

c--- muons
      do j=1,mxpart      
      if (     (plabel(j) == 'ml') .or. (plabel(j) == 'ma')) then
         ismuonarray(j) = .true. 
      else
         ismuonarray(j) = .false. 
      endif
      enddo

c--- neutrinos
      do j=1,mxpart      
      if (     (plabel(j) == 'nl') .or. (plabel(j) == 'na')
     &    .or. (plabel(j) == 'nm') .or. (plabel(j) == 'bm')
     &    .or. (plabel(j) == 'nt') .or. (plabel(j) == 'bt')) then
         isneutarray(j) = .true. 
      else
         isneutarray(j) = .false. 
      endif
      enddo

c--- dark matter
      do j=1,mxpart      
      if ((plabel(j) == 'xm') .or. (plabel(j) == 'xa')) then
         isdmarray(j) = .true. 
      else
         isdmarray(j) = .false. 
      endif
      enddo
      
      return
      end
      


      function is_hadronic(i)
       implicit none
      include 'types.f'
      logical:: is_hadronic
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'is_functions_com.f'
      integer:: i

c--- return value from array      
      is_hadronic=ishadarray(i)
      
      return 
      end


      function is_spechadronic(i)
       implicit none
      include 'types.f'
      logical:: is_spechadronic
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'is_functions_com.f'
      integer:: i

c--- return value from array      
      is_spechadronic=ishadspecarray(i)
      
      return 
      end


      function is_photon(i)
       implicit none
      include 'types.f'
      logical:: is_photon
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'is_functions_com.f'
      integer:: i

c--- return value from array      
      is_photon=isphotarray(i)
      
      return 
      end


      function is_lepton(i)
       implicit none
      include 'types.f'
      logical:: is_lepton
c--- note: checks for charged leptons only      
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'is_functions_com.f'
      integer:: i

c--- return value from array      
      is_lepton=isleptarray(i)
      
      return 
      end


      function is_electron(i)
       implicit none
      include 'types.f'
      logical:: is_electron
c--- note: checks for electrons only      
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'is_functions_com.f'
      integer:: i

c--- return value from array      
      is_electron=iselectronarray(i)
      
      return 
      end


      function is_muon(i)
       implicit none
      include 'types.f'
      logical:: is_muon
c--- note: checks for electrons only      
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'is_functions_com.f'
      integer:: i

c--- return value from array      
      is_muon=ismuonarray(i)
      
      return 
      end


      function is_neutrino(i)
       implicit none
      include 'types.f'
      logical:: is_neutrino
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'is_functions_com.f'
      integer:: i

c--- return value from array      
      is_neutrino=isneutarray(i)
      
      return 
      end


      function is_darkmatter(i)
       implicit none
      include 'types.f'
      logical:: is_darkmatter
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'is_functions_com.f'
      integer:: i

c--- return value from array      
      is_darkmatter=isdmarray(i)
      
      return 
      end
