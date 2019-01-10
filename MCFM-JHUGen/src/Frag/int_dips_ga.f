!---Integrated Dipoles for photon-fermion splitting --------------------------
!-                                                                           -
!---- Author C Williams Dec 2010 --------------------------------------------
!-                                                                           -   
!---- Only final initial and final-final are currently implemented since the -
!---- Photon is required to be a final state particle.                       -
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!                      Final - Initial  
!-----------------------------------------------------------------------------
!
! Corresponds to 5.175 of Catani Seymour with delta_ab=0 since b=photon a=fermion


      double precision function fi_gaq(x,p,L,a,b,vorz)
      implicit none
      include 'constants.f'
      include 'betacut.f'
      double precision x,p(mxpart,4),omx,n(4),L
      integer nu,vorz,j,a,b
      double precision n_sq,adb,adn,lmx
      double precision bit,invbit
      double precision dlg,L_ab,K_ab,V_ab,P_gaq

     
      fi_gaq=0d0
      do nu=1,4
      n(nu)=0d0
      enddo
      adn=0d0
      adb=0d0
      n_sq=0d0

      do j=1,2
         do nu=1,4
            n(nu)=n(nu)-p(j,nu)
         enddo
      enddo
            
     

!---- Determine momenta and invariants
!---- DEBUG CHANGED N-definition for diphoton
!      do j=3,mxpart
!         if(plabel(j) .eq. 'ga') then 
            Do nu=1,4              
               n(nu)=n(nu)-p(a,nu)
            enddo
!         endif
!      enddo
         
      do nu=1,4
         if (nu.ne.4) then 
            adb=adb-p(a,nu)*p(b,nu)
            adn=adn-p(a,nu)*n(nu)
            n_sq=n_sq-n(nu)*n(nu)
         else
            adb=adb+p(a,nu)*p(b,nu)
            adn=adn+p(a,nu)*n(nu)
            n_sq=n_sq+n(nu)*n(nu)
         endif
      enddo

c      write(6,*) adn,n_sq

      omx=one-x
      P_gaq=(one+omx**2)/x

C      write(6,*) 'read vorz',vorz
C      write(6,*) 'split ',P_gaq

      if((vorz .eq. 1).or.(vorz .eq. 3)) then        
         return
      elseif(vorz .eq. 2) then 
!---- split regular pieces into 3, V_ab = regular pieces of 5.82
!----- K_ab = regular pieces of 5.155 
!------L_ab = regular pieces of 5.176
         lmx=dlog(omx)
         dlg=dlog(-n_sq*adb/(two*adn**2))


         invbit=n_sq*x/(two*omx*adn)
         bit=1d0/invbit
              

!---- DROP 1/EPS FOR NOW (-epinv+L)=>L      
         V_ab=P_gaq*(lmx)+(L)*P_gaq+x

         K_ab=P_gaq*lmx 

         L_ab=-P_gaq*dlg

         fi_gaq=V_ab+K_ab+L_ab
!-- beta cut terms
       
        
         if(bit .gt. bfi) then
            fi_gaq=fi_gaq+P_gaq*dlog(invbit*bfi)
         endif

         return 
      endif
      return 
      end
         
!-----------------------------------------------------------------------------
!                      Final - Final 
!-----------------------------------------------------------------------------
!
! Corresponds to 5.132 of Catani Seymour 
      
      double precision function ff_gaq(z,L,vorz) 
      implicit none
      include 'constants.f'
      include 'betacut.f' 
      double precision z,L,omz,lz,lmz
      double precision P_gaq
      integer vorz

      omz=one-z
      lz=dlog(z)
      lmz=dlog(omz)
      ff_gaq=0d0

      P_gaq=(one+omz**2)/z

      if((vorz .eq. 1) .or. (vorz .eq.3)) then
         return
      elseif(vorz .eq. 2) then 
         ff_gaq=P_gaq*((L)+(lz+lmz))+z 
!----- beta cut terms 
         ff_gaq=ff_gaq+dlog(bff)*P_gaq 
         return
      endif


      
      return 
      end
