************************************************************
*   Transforms p1+p2 ----> p3+...pn into                   *
*              q1+q2 -----> q3 + qn-1                      *  
*  with an identified photon in the final state            * 
*   C Williams    Dec 2010                                 *
************************************************************


      subroutine transformfrag(p,q,z,ip,jp,kp) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'npart.f' 
      real(dp):: p(mxpart,4),q(mxpart,4),z,omz
      integer:: ip,jp,kp,j,nu,ipart
      real(dp):: n(4),idn,jdn,kdn
      real(dp):: pdotp,pdot_ar
      real(dp):: kt(4),ks(4),k(4),kdk,ksdks
      real(dp):: kdp(mxpart),ksdp(mxpart)
    
      

      do j=1,npart+2
         do nu=1,4
            q(j,nu)=0._dp
         enddo
      enddo

     
c---- Intial-Initial not coded
      if((ip <= 2) .and. (kp <= 2)) return 
      
c---- Initial-Final not coded
      if((ip <= 2) .and. (kp > 2)) return 
     
c---- Final-Initial 
     
      if((ip > 2) .and. (kp <= 2)) then 
         
         idn=0._dp
         jdn=0._dp
         kdn=0._dp
         do nu=1,4
            n(nu)=-p(1,nu)-p(2,nu)-p(ip,nu)
            if(nu .ne. 4) then 
               idn=idn-p(ip,nu)*n(nu)
               jdn=jdn-p(jp,nu)*n(nu)
               kdn=kdn-p(kp,nu)*n(nu) 
            else
               idn=idn+p(ip,nu)*n(nu)
               jdn=jdn+p(jp,nu)*n(nu)
               kdn=kdn+p(kp,nu)*n(nu)
            endif
         enddo

         

         z=idn/(jdn+idn)
        
          
         omz=one-z
            

         if (npart == 3) then 
           
c---- for 3 particles in the final state have i==photon j=final state emitted and one other FS mom 
c---- Need to keep kp fixed so FS mom gets shifted to conserve momentum 
            do nu=1,4
               q(1,nu)=p(1,nu)
               q(2,nu)=p(2,nu)
            enddo
            
            ipart=3
            do j=3,npart+2
               
               if (j == ip) then 
                  do nu=1,4
                     q(ipart,nu) = (one/z)*p(ip,nu)
                  enddo           
               elseif (j == jp) then 
                  goto 131 
               else
                  do nu=1,4
                     q(ipart,nu)=p(jp,nu)+p(j,nu)-(omz/z)*p(ip,nu)
                  enddo            
              endif
               ipart=ipart+1
 131           continue          
            enddo 
c            
         else 
           
                                
            k=0._dp
            kt=0._dp
            ks=0._dp
            
            
            kdk=0._dp
     
            do nu=1,4
               q(1,nu)=p(1,nu)
               q(2,nu)=p(2,nu)
               k(nu)=n(nu)-p(jp,nu)
               kt(nu)=n(nu)+(one-one/z)*p(ip,nu)
               ks(nu)=k(nu)+kt(nu)
               if(nu .ne.4) then 
                  kdk=kdk-kt(nu)*kt(nu)
               else
                  kdk=kdk+kt(nu)*kt(nu)
               endif
            enddo
           
            

            kdk=pdotp(k,k)
            ksdks=pdotp(ks,ks)
            
            ipart=3
            do j=3,npart+2
              
               if(j == jp) then 
                  goto 19
               elseif(j == ip) then 
                  do nu=1,4
                     q(ipart,nu)=one/z*p(ip,nu)
                  enddo
               else
                  kDp(j)=pdot_ar(k,p,j)
                  ksDp(j)=pdot_ar(ks,p,j)
                  do nu=1,4
                     q(ipart,nu)=p(j,nu)-two*ksDp(j)*ks(nu)/ksDks
     &                 +two*kdp(j)*kt(nu)/kdk                     
                  enddo              
               endif
               ipart=ipart+1
 19            continue
            enddo
         endif
      
         return 
         
c---  Final-Final 

      elseif((ip > 2) .and. (kp > 2)) then 

         if(npart==3) then 

!----- SPECIAL FORM FOR 2=>3 DIPOLE 
         ipart=1
      
      
         do j=1,npart+2
            if (j == ip) then 
               do nu=1,4
                  q(ipart,nu) = (one/z)*p(ip,nu)
               enddo           
            elseif (j == jp) then 
               goto 13
            elseif (j == kp) then 
               do nu=1,4
                  q(ipart,nu)=p(jp,nu)+p(kp,nu)-(omz/z)*p(ip,nu)
               enddo
            else
               do nu=1,4
                  q(ipart,nu)=p(j,nu)
               enddo
            endif
            ipart=ipart+1
 13         continue          
         enddo 
      

      else

         write(6,*) 'LORENTZ TRANSFORM NEEDS TO BE IMPLEMENTED' 
         stop 
      endif

      endif
     
      

      return 
      end

      function pdotp(p1,p2)
      implicit none
      include 'types.f'
      real(dp):: pdotp
      
      real(dp):: p1(4),p2(4)
      integer:: j
      pdotp=0._dp
      do j=1,3
         pdotp=pdotp-p1(j)*p2(j)
      enddo
      pdotp=pdotp+p1(4)*p2(4)
      return 
      end 

      function pdot_ar(p1,p2,k)
      implicit none
      include 'types.f'
      real(dp):: pdot_ar
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: k,j
      real(dp):: p1(4),p2(mxpart,4)
      pdot_ar=0._dp
      do j=1,3
         pdot_ar=pdot_ar-p1(j)*p2(k,j)
      enddo
      pdot_ar=pdot_ar+p1(4)*p2(k,4)
      return 
      end


      function checkv(p,ip,jp,kp)
       implicit none
      include 'types.f'
      logical:: checkv
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'betacut.f'
      real(dp):: p(mxpart,4),idj,idn
      real(dp):: vt,n(4),z,n_sq,omz,jdn,az
      integer:: ip,jp,nu,kp
      
      checkv=.true.
      
    
      
      idn=0._dp
      jdn=0._dp
      idj=0._dp
     
      do nu=1,4
         n(nu)=-p(1,nu)-p(2,nu)-p(ip,nu)
         if(nu .ne. 4) then 
            idn=idn-p(ip,nu)*n(nu)
            jdn=jdn-p(jp,nu)*n(nu)
            idj=idj-p(ip,nu)*p(jp,nu)
           
         else
            idn=idn+p(ip,nu)*n(nu)
            jdn=jdn+p(jp,nu)*n(nu)
            idj=idj+p(ip,nu)*p(jp,nu)
            
         endif
      enddo

      n_sq=0._dp
      do nu=1,4
         if(nu.ne.4) then 
            n_sq=n_sq-n(nu)*n(nu)
         else
            n_sq=n_sq+n(nu)*n(nu)
         endif
      enddo
      
      
      
      z=idn/(jdn+idn)          
      omz=one-z

      az=two*omz*idn/(z*n_sq)

    
      vt=idj/idn

      
     
      
      if(vt < az*bfi) then 
         checkv=.true. 
      else
         checkv=.false. 
      endif
      
      return 
      end


      
      function check_nv(p,ip,jp,kp)
       implicit none
      include 'types.f'
      logical:: check_nv
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'betacut.f'
      real(dp):: p(mxpart,4),idj,idn
      real(dp):: vt,n(4),n_sq,jdn 
      integer:: ip,jp,nu,kp
      
      check_nv=.true.
      
    
      
      idn=0._dp
      jdn=0._dp
      idj=0._dp
     
      do nu=1,4
         n(nu)=-p(1,nu)-p(2,nu)-p(ip,nu)
         if(nu .ne. 4) then 
            idn=idn-p(ip,nu)*n(nu)
            jdn=jdn-p(jp,nu)*n(nu)
            idj=idj-p(ip,nu)*p(jp,nu)
           
         else
            idn=idn+p(ip,nu)*n(nu)
            jdn=jdn+p(jp,nu)*n(nu)
            idj=idj+p(ip,nu)*p(jp,nu)
            
         endif
      enddo

      n_sq=0._dp
      do nu=1,4
         if(nu.ne.4) then 
            n_sq=n_sq-n(nu)*n(nu)
         else
            n_sq=n_sq+n(nu)*n(nu)
         endif
      enddo

      vt=idj/idn 
      
      

      if(vt>bfi) then 
         check_nv =.false.
         return 
      else
         check_nv=.true. 
      endif
      return 
      end
      
      
