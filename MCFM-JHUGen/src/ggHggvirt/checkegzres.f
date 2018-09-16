      subroutine checkEGZres
      implicit none
      include 'constants.f'
      include 'epinv.f'
      include 'masses.f'
      include 'scale.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'nflav.f'
      include 'scheme.f'
      include 'deltar.f'
      double precision p(mxpart,4),Hqarbvsqanal,Hqaqavsqanal,
     . HAQggvsqanal,Hggggvsqanal,check,check_a,check_b,check_i
      
c--- Set up momenta, Higgs mass and scale for special point

c--- These are per the paper
c      p(1,4)=+0.30674037867d0 
c      p(1,1)=-0.17738694693d0
c      p(1,2)=-0.01664472021d0
c      p(1,3)=-0.24969277974d0
c  
c      p(2,4)=+0.34445032281d0
c      p(2,1)=+0.14635282800d0
c      p(2,2)=-0.10707762397d0
c      p(2,3)=+0.29285022975d0
c    
c      p(3,4)=0.5d0
c      p(3,1)=0d0
c      p(3,2)=0d0
c      p(3,3)=0d0
c      
c      p(4,4)=0.5d0
c      p(4,1)=0d0
c      p(4,2)=0d0
c      p(4,3)=0d0
c      
c      p(5,4)=+0.22091667641d0
c      p(5,1)=+0.08911915938d0
c      p(5,2)=+0.19733901856d0
c      p(5,3)=+0.04380941793d0
c    
c      p(6,4)=+0.12789262211d0
c      p(6,1)=-0.05808504045d0
c      p(6,2)=-0.07361667438d0
c      p(6,3)=-0.08696686795d0

c--- These are direct from Giulia's code (more s.f.)
      p(1,4)=+0.306740378669096d0    
      p(1,1)=-0.1773869469327211d0
      p(1,2)=-0.01664472021061108d0
      p(1,3)=-0.2496927797375384d0
  
      p(2,4)=+0.3444503228094516d0
      p(2,1)=+0.1463528280021222d0
      p(2,2)=-0.1070776239705892d0
      p(2,3)=+0.2928502297491496d0
    
      p(3,4)=0.5d0
      p(3,1)=0d0
      p(3,2)=0d0
      p(3,3)=0d0
      
      p(4,4)=0.5d0
      p(4,1)=0d0
      p(4,2)=0d0
      p(4,3)=0d0
      
      p(5,4)=+0.2209166764078266d0
      p(5,1)=+0.08911915938215637d0
      p(5,2)=+0.1973390185605316d0
      p(5,3)=+0.04380941793341205d0
    
      p(6,4)=+0.1278926221136257d0
      p(6,1)=-0.05808504045155738d0
      p(6,2)=-0.0736166743793313d0
      p(6,3)=-0.0869668679450233d0
      
      hmass=1d0
      scale=1d0
      musq=1d0
      nflav=5
      scheme='tH-V'
      deltar=1d0
      
      call dotem(6,p,s)
      call spinoru(6,p,za,zb)

c--- only comparing finite piece      
      epinv=0d0
       
      call HqarbLO(6,2,5,1,check)
      write(6,99) 'Proc A: EGZ check B',check
      check=Hqarbvsqanal(5,6,1,2)
      write(6,99) 'Proc A: EGZ check V',check

      call HqaqaLO(6,2,5,1,check,check_a,check_b,check_i)
      write(6,99) 'Proc B: EGZ check B',check
      check=Hqarbvsqanal(5,6,1,2)+Hqarbvsqanal(1,6,5,2)
     .     +Hqaqavsqanal(5,6,1,2)
      write(6,99) 'Proc B: EGZ check V',check
      
      call HQAggLO(1,2,5,6,check,check_a,check_b,check_i)
      write(6,99) 'Proc C: EGZ check B',check
      check=HAQggvsqanal(1,2,5,6)
      write(6,99) 'Proc C: EGZ check V',check
      
      call HggggLO(1,2,5,6,check,check_a,check_b,check_i)
      write(6,99) 'Proc D: EGZ check B',check
      check=Hggggvsqanal(1,2,5,6)
      write(6,99) 'Proc D: EGZ check V',check
      
c      pause
      
      return
      
   99 format(a20,e24.15)
   
      end
      
