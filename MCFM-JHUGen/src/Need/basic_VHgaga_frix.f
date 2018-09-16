c--- This subroutine is a basic implementation of hep-ph/9801442 valid for 
c---- f(p1)+f(p2) -->V->(l(p3)+la(p4))+gamma(p5)+gamma(p6)+f(p7) C Williams July 11

      subroutine basic_VHgaga_frix(pin,passed,isub) 
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
      real(dp):: R57,R67,ET5,ET6,ET7,alpha_5,alpha_6

      passed=.true.
      
c--- This routine also called for 'gmgmjt' -> always need second branch
      if ( ((kpart==klord) .or. (kpart==kvirt)
     &.or. (isub == 1) ) .and.
     &  ((kcase==kWHgaga) .or. (kcase==kZHgaga)) ) then
c---- LO, virtual and dipole pieces for gamgam and Higaga have no jets --> passed cut
         return
 
c---  Real radiation (isub=0) and 'gmgmjt' should be checked
      elseif (isub == 0) then 
         
c---- Calculate angles and energies
         ET5=ret_ET(pin,5)
         ET6=ret_ET(pin,6)
         ET7=ret_ET(pin,7)
         
         R57=R(pin,5,7)
         R67=R(pin,6,7)
         
         alpha_5=ET5/(1._dp-cos(cone_ang)) 
         alpha_6=ET6/(1._dp-cos(cone_ang))
         
         if((R57 < cone_ang) .and. (R67 < cone_ang)) then  
            passed = in_cone(R57,ET7,alpha_5) 
            if (passed) then 
               passed = in_cone(R67,ET7,alpha_6)
            endif
            return 
         elseif (R57 < cone_ang) then 
            passed=in_cone(R57,ET7,alpha_5) 
            return 
         elseif (R67 < cone_ang) then 
            passed=in_cone(R67,ET7,alpha_6) 
            return 
         else
            passed = .true. 
            return 
         endif
                        
      endif

      return 
      end
      
