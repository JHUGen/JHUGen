      subroutine init_is_functions()
      implicit none
      include 'constants.f'
      include 'is_functions_com.f'
      include 'plabel.f'
      integer j

c--- hadrons      
      do j=1,mxpart      
      if ( (plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'pj')
     & .or.(plabel(j) .eq. 'bq') .or. (plabel(j) .eq. 'ba')
     & .or.(plabel(j) .eq. 'qj') ) then
         ishadarray(j) = .true. 
      else
         ishadarray(j) = .false. 
      endif
      enddo

c--- photons
      do j=1,mxpart      
      if (plabel(j) .eq. 'ga') then
         isphotarray(j) = .true. 
      else
         isphotarray(j) = .false. 
      endif
      enddo

c--- charged leptons
      do j=1,mxpart      
      if (     (plabel(j) .eq. 'el') .or. (plabel(j) .eq. 'ea')
     &    .or. (plabel(j) .eq. 'ml') .or. (plabel(j) .eq. 'ma')
     &    .or. (plabel(j) .eq. 'tl') .or. (plabel(j) .eq. 'ta')) then
         isleptarray(j) = .true. 
      else
         isleptarray(j) = .false. 
      endif
      enddo

c--- electrons
      do j=1,mxpart      
      if (     (plabel(j) .eq. 'el') .or. (plabel(j) .eq. 'ea')) then
         iselectronarray(j) = .true. 
      else
         iselectronarray(j) = .false. 
      endif
      enddo

c--- muons
      do j=1,mxpart      
      if (     (plabel(j) .eq. 'ml') .or. (plabel(j) .eq. 'ma')) then
         ismuonarray(j) = .true. 
      else
         ismuonarray(j) = .false. 
      endif
      enddo

c--- neutrinos
      do j=1,mxpart      
      if (     (plabel(j) .eq. 'nl') .or. (plabel(j) .eq. 'na')
     &    .or. (plabel(j) .eq. 'nm') .or. (plabel(j) .eq. 'bm')
     &    .or. (plabel(j) .eq. 'nt') .or. (plabel(j) .eq. 'bt')) then
         isneutarray(j) = .true. 
      else
         isneutarray(j) = .false. 
      endif
      enddo

c--- dark matter
      do j=1,mxpart      
      if ((plabel(j) .eq. 'xm') .or. (plabel(j) .eq. 'xa')) then
         isdmarray(j) = .true. 
      else
         isdmarray(j) = .false. 
      endif
      enddo
      
      return
      end
      


      logical function is_hadronic(i)
      implicit none
      include 'constants.f'
      include 'is_functions_com.f'
      integer i

c--- return value from array      
      is_hadronic=ishadarray(i)
      
      return 
      end


      logical function is_photon(i)
      implicit none
      include 'constants.f'
      include 'is_functions_com.f'
      integer i

c--- return value from array      
      is_photon=isphotarray(i)
      
      return 
      end


      logical function is_lepton(i)
c--- note: checks for charged leptons only      
      implicit none
      include 'constants.f'
      include 'is_functions_com.f'
      integer i

c--- return value from array      
      is_lepton=isleptarray(i)
      
      return 
      end


      logical function is_electron(i)
c--- note: checks for electrons only      
      implicit none
      include 'constants.f'
      include 'is_functions_com.f'
      integer i

c--- return value from array      
      is_electron=iselectronarray(i)
      
      return 
      end


      logical function is_muon(i)
c--- note: checks for electrons only      
      implicit none
      include 'constants.f'
      include 'is_functions_com.f'
      integer i

c--- return value from array      
      is_muon=ismuonarray(i)
      
      return 
      end


      logical function is_neutrino(i)
      implicit none
      include 'constants.f'
      include 'is_functions_com.f'
      integer i

c--- return value from array      
      is_neutrino=isneutarray(i)
      
      return 
      end


      logical function is_darkmatter(i)
      implicit none
      include 'constants.f'
      include 'is_functions_com.f'
      integer i

c--- return value from array      
      is_darkmatter=isdmarray(i)
      
      return 
      end
