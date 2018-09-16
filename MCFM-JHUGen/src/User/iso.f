!--- This is a driver function for photon isolation, the relevant photon isolation 
!--- critera are selected and applied 

      function iso(p,phot_id,isub,nd)
       implicit none
      include 'types.f'
      logical:: iso
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'frag.f'
      include 'useet.f'
      include 'phot_dip.f'
      include 'z_dip.f'
      include 'lastphot.f'
      real(dp):: p(mxpart,4)
      integer:: isub,nd,imode
      integer:: phot_id ! refers to which photon we are isolating
      logical:: first 
      logical:: photo_iso_phys,photo_iso_z!,photo_iso,iso_old
      data first/.true./
      save first
!$omp threadprivate(first)

c--- imode: set the behaviour of the isolation requirement (default: imode=0)
c---   imode=0 : scaling cut for epsilon_h < 1, fixed cut for epsilon_h > 1
c---   imode=1 : scaling cut for all epsilon_h
c---   imode=2 : fixed cut for all epsilon_h
      imode=0     

      iso=.true.

! Check if no isolation required
      if ((abs(cone_ang) < 1.e-4_dp).or.(abs(epsilon_h) < 1.e-4_dp)) then
        if (first) then
        write(6,*)'****************************************************'
        write(6,*)'*                                                  *'
        write(6,*)'*         No photon isolation cuts applied         *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        first=.false.
      endif
        return
      endif

c--- for imode=0, decide which cut to use depending on epsilon_h
!--   epsilon_h < 1 : epsilon_h corresponds to a pt fraction i.e. pt(had) < epsilon_h pt_gamma 
!--   epsilon_h > 1 : treat it as an E_t max in cone   i.e. pt(had) < epsilon_h 
       if (imode == 0) then
         if (epsilon_h < 0.9999_dp) then
         imode=1
         else
         imode=2
         endif
       endif

!===== NEW ISOLATION : WORK OUT WHICH ROUTINE TO CALL 
 
!========== are we doing fragmentation and is this the
!========== photon produced by fragmentation? 
       if((fragint_mode) .and. (phot_id == lastphot)) then 
          iso=photo_iso_z(p,phot_id,z_frag,imode,isub) 
!========== are we doing a photon dipole and is this the
!========== photon produced by fragmentation? 
       elseif(phot_dip(nd) .and. (phot_id == lastphot)) then 
          iso=photo_iso_z(p,phot_id,z_dip(nd),imode,isub) 
!========== we are in default mode, BORN, virt, real, non-photon dipole 
       else
          iso=photo_iso_phys(p,phot_id,imode,isub) 
       endif

!===== OLD ISOLATION 
c       iso_old=photo_iso(p,isub,phot_dip(nd),phot_id,nd,imode)       
c       if (iso .neqv. iso_old) then
c         write(6,*)'WARNING:',iso,'(new) vs',iso_old,' (old) for nd=',nd
c       else
c         write(6,*) 'OKAY'
c       endif
       
c--- write out isolation parameters
       if    (first .and. (imode == 1)) then 
        write(6,*)'************** Photons Isolated     ****************'
        write(6,*)'*                                                  *'
        write(6,99)'*    E_t(had) in cone',cone_ang,' < ',epsilon_h,
     &   ' E_t(phot)     *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        first=.false.
      elseif (first .and. (imode == 2)) then
        write(6,*)'************** Photons Isolated     ****************'
        write(6,*)'*                                                  *'
        write(6,96)'* E_t (had) in cone',cone_ang,' < ',epsilon_h,
     &   'GeV    *'
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        first=.false.
      endif
      
!      endif

      return 

 99   format(1x,a21,f6.2,a3,f6.2,a16)
 96   format(1x,a19,f6.2,a4,f6.2,a17)
      end
