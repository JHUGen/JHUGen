c--- This subroutine is a basic implementation of hep-ph/9801442 valid for 
c---- f(p1)+f(p2) --> gamma(p3)+gamma(p4)+f(p5) C Williams Dec 2010

      subroutine basic_di_frix(pin,passed,isub) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'kprocess.f'
      include 'frag.f'
      include 'kpart.f'
      real(dp):: pin(mxpart,4)
      integer:: isub
      logical:: passed,in_cone
      real(dp):: ret_ET,R
      real(dp):: R45,R35,ET3,ET4,ET5,alpha_3,alpha_4

      passed=.true.
      
c--- This routine also called for 'gmgmjt' -> always need second branch
      if ( ((kpart==klord) .or. (kpart==kvirt)
     &.or. (isub == 1) ) .and.
     &     ((kcase==kgamgam) .or. (kcase==kgg2gam)
     & .or. (kcase==kHigaga)) ) then
c---- LO, virtual and dipole pieces for gamgam and Higaga have no jets --> passed cut
         return
 
c---  Real radiation (isub=0) and 'gmgmjt' should be checked
      elseif (isub == 0) then 
         
c---- Calculate angles and energies
         ET3=ret_ET(pin,3)
         ET4=ret_ET(pin,4)
         ET5=ret_ET(pin,5)
         
         R35=R(pin,3,5)
         R45=R(pin,4,5)
         
         alpha_3=ET3/(1._dp-cos(cone_ang)) 
         alpha_4=ET4/(1._dp-cos(cone_ang))
         
         if((R35 < cone_ang) .and. (R45 < cone_ang)) then  
            passed = in_cone(R35,ET5,alpha_3) 
            if (passed) then 
               passed = in_cone(R45,ET5,alpha_4)
            endif
            return 
         elseif (R35 < cone_ang) then 
            passed=in_cone(R35,ET5,alpha_3) 
            return 
         elseif (R45 < cone_ang) then 
            passed=in_cone(R45,ET5,alpha_4) 
            return 
         else
            passed = .true. 
            return 
         endif
                        
      endif

      return 
      end
      
