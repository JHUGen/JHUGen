      function photo_iso_phys(p,phot_id,imode,isub)
       implicit none
      include 'types.f'
      logical:: photo_iso_phys
       
!----------------------------------------------------------------------------
!-                NEW!       Photon Isolation                               -
!- C. Williams April 2013                                                    - 
!-                                                                          -
!- This function implements isolation cuts on photons                       -
!- The general requirement is that pt_had < epsilon_h pt_photon             -  
!- Within a cone of radius cone_ang                                         - 

!---- this new version ONLY WORKS with physical momenta, i.e p should correspond 
!---- to either a Born, Virt or real phase space point. 
!---- in adition, dipoles which do not arise from a collinear photon can also use this routine 
!---- as such the isolation is solely defined in terms of the partons in a cone
!---- around the photon there is no need for knowledge of z_frag or z_dip
!----------------------------------------------------------------------------
 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'frag.f'
      include 'npart.f'
      real(dp):: p(mxpart,4) 
      integer:: imode,phot_id 
      real(dp):: z_c,opeps,Rjga,pt,p_incone(4),pt_inc
      real(dp):: z_kin,R
      logical:: is_hadronic
      integer:: j,nu,isub
!------ initialize 
      photo_iso_phys = .true. 
      p_incone(:)=0._dp 
      

!====== Parameter from input file
      opeps=one+epsilon_h 
      
      if(imode==1) then 
!======= this is scaling isolation E_T_max < epsilon_h * pt_gamma 
         z_c=one/opeps 
      elseif(imode==2) then 
!======= this is fixed cut  i.e. E_T_max < 10 GeV etc. 
         z_c=pt(phot_id,p)/(epsilon_h+pt(phot_id,p))
      else 
         write(6,*) 'Unknown isolation parameter : imode = ',imode
      endif
      
!========== Calculate the hadronic four-momenta in the cone 
      do j=3,npart+2-isub
         if(is_hadronic(j)) then 
            Rjga=R(p,j,phot_id) 
            if(Rjga < cone_ang) then 
               do nu=1,4
                  p_incone(nu)=p_incone(nu)+p(j,nu) 
               enddo 
            endif
         endif
      enddo
!====== calculate the PT of this quantity 
      pt_inc=sqrt(p_incone(1)**2+p_incone(2)**2) 

!====== now calculate the hadronic energy fraction of the photon pt 
      z_kin = pt(phot_id,p)/(pt_inc+pt(phot_id,p))
      
c      write(6,*) 'photo_iso_phys: isub=',isub
c      write(6,*) 'photo_iso_phys: phot_id=',phot_id
c      write(6,*) 'photo_iso_phys: pt(phot_id,p)',pt(phot_id,p)
c      write(6,*) 'photo_iso_phys: pt_inc',pt_inc
c      write(6,*) 'photo_iso_phys: z_kin=',z_kin
c      write(6,*) 'photo_iso_phys: z_c  =',z_c

!====== Finally compare to z_c    
      if(z_kin < z_c) then 
            photo_iso_phys = .false. 
!            write(6,*) z_kin,z_c,pt_inc,pt(phot_id,p),Rjga
!            call writeout(p) 
!            pause
         endif
      return 
      end 

      

      
