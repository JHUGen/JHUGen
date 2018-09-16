c---- Perturbative Fragmentation functions 
     
      subroutine NP_fragsetI(z,Msq,parton_id,NP_part) 
      
      implicit none

      double precision z,Msq
      integer parton_id 
      double precision NP_part

      include 'zMsq_grid.f'
      integer n_spec
      parameter(n_spec=6)
      integer i,j
      double precision XQDUM1(num_z,num_M2_4Flav,n_spec)
      double precision XQDUM2(num_z,num_M2_5Flav,n_spec) 

c---- Grid of BFG

      include 'gridNonPert__setI.f'

     
      Integer IQMax,IQ,it_p
      parameter(IQMax=22,it_p=4)
      double precision X(num_z)
      double precision Y(IQMax),Y_4flav(IQMax),Y_5flav(IQMax) 
      double precision F(num_z,IQMax) 
c      double precision M2(IQMax) 
      double precision Z2A(it_p),Z1A(it_p),OA(it_p,it_p)
      double precision a,b,out_1,out_2
      integer jz,jy

      double precision zmin,zmax,Mmin
 
      logical first,no_b
      first=.true.
    

      if(first) then 
         do i=1,num_z
            X(i)=log10(z_grid(i))
         enddo 
         do i=1,num_M2_4Flav
            Y_4flav(i)=log10(M2_4Flav(i))
         enddo 
         do i=1,num_M2_5Flav
            Y_5flav(i)=log10(M2_5Flav(i))
         enddo
      first = .false.
      endif

      zmin=z_grid(1) 
      zmax=z_grid(num_z)
      
c--- Check that M**2 isnt too small 
      if(Msq .lt. 20.26d0) then 
         Mmin = M2_4Flav(1)
         no_b = .true.
      else
         Mmin = M2_5Flav(1) 
         no_b = .false. 
      endif 
      if (Msq .lt. Mmin) then 
         write(6,*) 'WARNING ',Msq,' is too small ' 
         stop
      endif

      if(no_b) then 
         IQ = num_M2_4Flav
         Do j=1,IQ 
            do i=1,num_z
               F(i,j)=XQDUM1(i,j,parton_id)
            enddo
            Y(j)=Y_4Flav(j) 
         enddo
      else
         IQ = num_M2_5Flav
         Do j=1,IQ 
            do i=1,num_z
               F(i,j)=XQDUM2(i,j,parton_id)
            enddo
            Y(j)=Y_5Flav(j) 
         enddo
      endif
      

      NP_part = 0.0
     
        
  
      
      if((z .gt. zmin) .and. (z .lt. zmax)) then 
         a=log10(z)
         b=log10(Msq)
         call locate(X,num_z,a,jz)
         call locate(Y,IQ,b,jy) 
         do i=1,it_p
c-----      Special Case where jz = 1
            if (jz .eq. 1) then 
               Z1A(i)=X(i) 
               do j=1,it_p
                  if (jy .eq. 1) then 
                     Z2A(j)=Y(j)
                     OA(i,j)=F(i,j) 
                  elseif (jy .eq. (IQ-1)) then 
                     Z2A(j)=Y(jy-3+j)
                     OA(i,j)=F(i,jy-3+j) 
                  elseif (jy .eq. IQ) then 
                     Z2A(j)=Y(jy-4+j) 
                     OA(i,j)=F(i,jy-4+j) 
                  else
                     Z2A(j) =Y(jy-2+j)
                     OA(i,j)=F(i,jy-2+j) 
                  endif
               enddo
c-----      Special case where jz = num_z -1
            elseif (jz .eq. (num_z-1)) then
               Z1A(i)=X(jz-3+i) 
               do j=1,it_p
                  if (jy .eq. 1) then  
                     Z2A(j)=Y(j)
                     OA(i,j)=F(jz-3+i,j)
                  elseif (jy .eq. (IQ-1)) then 
                     Z2A(j)=Y(jy-3+j)
                     OA(i,j)=F(jz-3+i,jy-3+j)
                  elseif (jy .eq. IQ) then 
                     Z2A(j) = Y(jy-4+j) 
                     OA(i,j) = F(jz-3+i,jy-3+j)
                  else 
                     Z2A(j)= Y(jy-2+j)
                     OA(i,j)=F(jz-3+i,jy-2+j)
                  endif
               enddo
c-----         Special case where jz = num_z 
            elseif (jz .eq. num_z) then 
               Z1A(i)=X(jz-4+i)
               do j=1,it_p
                  if (jy .eq. 1) then 
                     Z2A(j)=Y(j)
                     OA(i,j) = F(jz-4+i,j)
                  elseif (jy .eq. (IQ-1)) then 
                     Z2A(j) = Y(jy-3+j) 
                     OA(i,j) = F(jz-4+i,jy-3+j) 
                  elseif (jy .eq. IQ ) then 
                     Z2A(j) = Y(jy-4+j) 
                     OA(i,j) = F(jz-4+i,jy-3+j) 
                  else
                     Z2A(j)=Y(jy-2+j) 
                     OA(i,j) = F(jz-4+i,jy-3+j)
                  endif
               enddo
c----          General case
            else
               Z1A(i)=X(jz-2+i) 
               do j=1,it_p
                  if (jy .eq. 1) then 
                     Z2A(j)=Y(j)
                     OA(i,j)=F(jz-2+i,j)
                  elseif (jy .eq. (IQ-1)) then 
                     Z2A(j)=Y(jy-3+j)
                     OA(i,j)=F(jz-2+i,jy-3+j)
                  elseif (jy .eq. IQ) then 
                     Z2A(j)=Y(jy-4+j)
                     OA(i,j) = F(jz-2+i,jy-4+j)
                  else
                     Z2A(j) = Y(jy-2+j)
                     OA(i,j) = F(jz-2+i,jy-2+j)
                  endif
               enddo
            endif
         enddo


            call DPOLIN2(Z1A,Z2A,OA,4,4,A,B,out_1,out_2) 
            NP_part=out_1 
           
      endif

      
      end subroutine 

      
