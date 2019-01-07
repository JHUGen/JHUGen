      logical function photo_iso_z(p,phot_id,z,imode,isub)
      implicit none 
!----------------------------------------------------------------------------
!-                NEW!       Photon Isolation                               -
!- C. Williams April 2013                                                    - 
!-                                                                          -
!- This function implements isolation cuts on photons                       -
!- The general requirement is that pt_had < epsilon_h pt_photon             -  
!- Within a cone of radius cone_ang                                         - 

!---- this new version ONLY WORKS with some kind of photon dipole or fragmentaiton 
!, i.e p should correspond to a photon dipole or fragmentation phase space point. 
!---- to either a Born, Virt or real phase space point. 
!---- as such the isolation is solely defined in terms of the partons in a cone
!---- around the photon there is no need for knowledge of z_frag or z_dip
!----------------------------------------------------------------------------
      include 'constants.f'
      include 'frag.f'
      include 'npart.f'
      double precision p(mxpart,4) 
      integer imode,phot_id 
      double precision z_c,opeps,Rjga,pt,p_incone(4),pt_inc
      double precision x_had,R,z
      logical is_hadronic
      integer j,nu,isub
!------ initialize 
      photo_iso_z = .true. 
      p_incone(:)=0d0 

!====== Parameter from input file
      opeps=one+epsilon_h 
      
      if(imode.eq.1) then 
!======= this is scaling isolation: E_T_max < epsilon_h * pt_gamma 
         z_c=one/opeps 
      elseif(imode.eq.2) then 
!======= this is fixed cut  i.e. E_T_max < 10 GeV etc. 
         z_c=pt(phot_id,p)/(epsilon_h+pt(phot_id,p))
      else 
         write(6,*) 'Unknown isolation parameter : imode = ',imode
      endif

!========== Calculate the hadronic four-momenta in the cone 
      
      do j=3,npart+2-isub
         if(is_hadronic(j)) then 
            Rjga=R(p,j,phot_id) 
            if(Rjga .lt. cone_ang) then 
               do nu=1,4
                  p_incone(nu)=p_incone(nu)+p(j,nu) 
               enddo 
            endif
         endif
      enddo
!====== calculate the PT of this quantity 
      pt_inc=dsqrt(p_incone(1)**2+p_incone(2)**2) 

!===== NOW DEFINE x_had (this is different from the photo_iso routine, 
!====== see eq. (5.9) of hep-ph/0204023)
      x_had=pt(phot_id,p)/(pt(phot_id,p)+z*pt_inc)
!===== nb for the vast majority of the time there will be no such radiation in cone, in this
!===== case pt_inc=0d0 and x_had =1 

c       write(6,*) 'photo_iso_z: phot_id=',phot_id
c       write(6,*) 'photo_iso_z: pt(phot_id,p)',pt(phot_id,p)
c       write(6,*) 'photo_iso_z: remnant: ',pt(phot_id,p)*(1d0/z-1d0)
c       write(6,*) 'photo_iso_z: z,pt_inc',z,pt_inc
c       write(6,*) 'photo_iso_z: x_had',x_had
c       write(6,*) 'photo_iso_z: z_c',z_c

!==== constraint is now on the product of x and z 
      if(x_had*z.lt.z_c) then 
         photo_iso_z=.false. 
!         write(6,*) z_c,pt_inc,pt(phot_id,p),Rjga
!         write(6,*) phot_id,x_had,z
!            call writeout(p) 
!            pause
   
      endif
      
      return 
      end 
