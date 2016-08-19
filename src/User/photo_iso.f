      logical function photo_iso(p,isub,phot_dip,phot_id,nd,imode) 
      implicit none
!----------------------------------------------------------------------------
!-                           Photon Isolation                               -
!- C. Williams Jan 2011                                                     - 
!-                                                                          -
!- This function implements isolation cuts on photons                       -
!- The general requirement is that pt_had < epsilon_h pt_photon             -  
!- Within a cone of radius cone_ang                                         - 
!- z_c = 1/1+epsilon_h corresponds to the lower cut off in z_frag           -
!----------------------------------------------------------------------------
      include 'constants.f'
      include 'frag.f'
      include 'npart.f'
      include 'z_dip.f'
      double precision p(mxpart,4),p_incone(4),R,pt_incone,pt 
      double precision z_c,opepsilon_h,z_kin,Rjga,tiny
      integer isub,j,nu,nd,imode
      logical phot_dip
      integer phot_id ! refers to which photon we are isolating
      logical is_hadronic 
      parameter(tiny=1d-8)

      photo_iso = .true. 

      do nu=1,4
         p_incone(nu)=0d0
      enddo
      z_kin = 0d0
      pt_incone = 0d0 

      opepsilon_h = one+epsilon_h

      if     (imode .eq. 1) then
!--- isolation using scaling cut 
         z_c = one/opepsilon_h
      elseif (imode .eq. 2) then
!--- isolation using fixed cut
         z_c=pt(phot_id,p)/(epsilon_h+pt(phot_id,p))
c--- This code removed pending improvements
c        if(phot_dip.eqv..true.) then 
c!---- Photon dipole need to rescale pt 
c           z_c=(z_dip(nd)*pt(phot_id,p))
c     &        /(epsilon_h+(z_dip(nd)*pt(phot_id,p)))
c        elseif(z_frag.gt.tiny) then 
c           z_c=(z_frag*pt(phot_id,p))
c     &        /(epsilon_h+(z_frag*pt(phot_id,p)))
c        else 
c           z_c=pt(phot_id,p)/(epsilon_h+pt(phot_id,p))
c        endif
      else
        write(6,*) 'Unknown isolation parameter: imode=',imode
        stop
      endif
         

!---- Define hadronic four-momentum (and pt) in cone 
      do j=3,npart+2-isub
         if(is_hadronic(j)) then 
            Rjga=R(p,j,phot_id) 
            if(Rjga .lt. cone_ang) then 
               do nu=1,4
                  p_incone(nu)=p_incone(nu)+p(j,nu) 
               enddo
               pt_incone=pt_incone+pt(j,p) 
            endif
         endif
      enddo
c--- CHECK: definition of pt_incone above??????? CHECK

!---- isub = 0 Can have (currently in MCFM - tree level Fragmentation or NLO Direct) 
!----- for Frag z_frag > 0.001d0 use this to separate pieces
      if((isub .eq. 0) .and. (z_frag .lt. tiny)) then 
!---- LO/NLO Direct
         z_kin = pt(phot_id,p)/(pt_incone+pt(phot_id,p))
         
         if(z_kin .lt. z_c) then 
            photo_iso = .false. 
            
            return 
         endif
         
      elseif((isub.eq.0) .and. (z_frag .gt. tiny)) then 
!---- Frag, case 1 no radiation in_cone only check z  
         if(pt_incone .lt. tiny) then 
            if(z_frag .lt. z_c) then 
               photo_iso = .false. 
               return 
            endif
         else
!---- Radiation in cone ! currently never used, need to check when ness
            z_kin=z_frag*pt(phot_id,p)/(pt(phot_id,p)+z_frag*pt_incone)
            if(z_kin .lt. z_c) then 
               photo_iso =.false.
               return 
            endif
         endif
         
      elseif((isub .eq. 1) .and. (phot_dip .eqv. .false.)) then 
!---- isub = 1 phot_dip = .false. z is calculated from parton kinematics         
         
         z_kin = pt(phot_id,p)/(pt_incone+pt(phot_id,p))
      
         if(z_kin .lt. z_c) then 
            photo_iso = .false. 
            return 
         endif
         
      elseif((isub .eq. 1) .and. (phot_dip .eqv. .true.)) then 
!---- isub = 1 phot_dip = .true. z_frag = z
!     case 1 no radiation in_cone only check z  
         if(pt_incone .lt. tiny) then 
            if(z_dip(nd) .lt. z_c) then 
               photo_iso = .false. 
               return 
            endif
         else
!---- Radiation in cone !  
            z_kin=pt(phot_id,p)/(pt(phot_id,p)+pt_incone)
            if(z_kin .lt. z_c) then 
               photo_iso =.false.             
               return 
            endif
         endif
         
      
      endif
      
      return 
      end

      
