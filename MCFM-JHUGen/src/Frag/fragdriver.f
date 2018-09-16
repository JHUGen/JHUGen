c---- This is a wrapper routine for the fragmentation functions 
c---- Given M^2 and z as inputs the routine returns D(z,M^2) 
c---- In addition one can set choose either Set I or Set II of BFG (hep-ph/9704447)

      subroutine get_frag(z,Msq2,set_flag,parton_flag,frag_func) 
      implicit none
      include 'types.f'
       
      real(dp):: z,Msq2 
      integer:: set_flag, parton_id,parton_flag,j
      real(dp):: frag_func
      real(dp):: pert_part,NP_part
      integer,parameter::parton_mcfm(6)=(/0,1,2,3,4,5/) 
      integer,parameter::parton_BFG(6)=(/1,3,2,4,5,6/)

      parton_id=100

c--- Convert MCFM notation for partons into BFG 
      do j=1,6 
         if(parton_flag == parton_mcfm(j)) then 
            parton_id = parton_BFG(j) 
         endif
      enddo

      if(parton_id > 6) then 
         write(6,*) 'WARNING PARTON NOT RECOGNISED AS INPUT' 
      endif

c---- Perturbative part 
      
      Call pert_frag(z,Msq2,parton_id,pert_part) 
      
c---- Non Perturbative part 
      if(set_flag == 1) then 
         Call NP_fragsetI(z,Msq2,parton_id,NP_part)
      elseif(set_flag == 2) then 
         Call NP_fragsetII(z,Msq2,parton_id,NP_part)
      else 
         write(6,*) 'WARNING no NP set dectected' 
         NP_part=0._dp
      endif

c---- reset parton id 
      parton_id=100
c---- Output      
      frag_func=pert_part+NP_part 
  
      return
      end subroutine 
      
