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


      function fi_gaq(x,p,L,a,b,vorz)
      implicit none
      include 'types.f'
      real(dp):: fi_gaq
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'betacut.f'
      real(dp):: x,p(mxpart,4),omx,n(4),L
      integer:: nu,vorz,j,a,b
      real(dp):: n_sq,adb,adn,lmx
      real(dp):: bit,invbit
      real(dp):: dlg,L_ab,K_ab,V_ab,P_gaq

     
      fi_gaq=0._dp
      do nu=1,4
      n(nu)=0._dp
      enddo
      adn=0._dp
      adb=0._dp
      n_sq=0._dp

      do j=1,2
         do nu=1,4
            n(nu)=n(nu)-p(j,nu)
         enddo
      enddo
            
     

!---- Determine momenta and invariants
!---- DEBUG CHANGED N-definition for diphoton
!      do j=3,mxpart
!         if(plabel(j) == 'ga') then 
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

      if((vorz == 1).or.(vorz == 3)) then        
         return
      elseif(vorz == 2) then 
!---- split regular pieces into 3, V_ab = regular pieces of 5.82
!----- K_ab = regular pieces of 5.155 
!------L_ab = regular pieces of 5.176
         lmx=log(omx)
         dlg=log(-n_sq*adb/(two*adn**2))


         invbit=n_sq*x/(two*omx*adn)
         bit=1._dp/invbit
              

!---- DROP 1/EPS FOR NOW (-epinv+L)=>L      
         V_ab=P_gaq*(lmx)+(L)*P_gaq+x

         K_ab=P_gaq*lmx 

         L_ab=-P_gaq*dlg

         fi_gaq=V_ab+K_ab+L_ab
!-- beta cut terms
       
        
         if(bit > bfi) then
            fi_gaq=fi_gaq+P_gaq*log(invbit*bfi)
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
      
      function ff_gaq(z,L,vorz) 
      implicit none
      include 'types.f'
      real(dp):: ff_gaq
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'betacut.f' 
      real(dp):: z,L,omz,lz,lmz
      real(dp):: P_gaq
      integer:: vorz

      omz=one-z
      lz=log(z)
      lmz=log(omz)
      ff_gaq=0._dp

      P_gaq=(one+omz**2)/z

      if((vorz == 1) .or. (vorz ==3)) then
         return
      elseif(vorz == 2) then 
         ff_gaq=P_gaq*((L)+(lz+lmz))+z 
!----- beta cut terms 
         ff_gaq=ff_gaq+log(bff)*P_gaq 
         return
      endif


      
      return 
      end
