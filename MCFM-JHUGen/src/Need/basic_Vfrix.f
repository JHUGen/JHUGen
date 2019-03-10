c--- This subroutine is a basic implementaion of hep-ph/9801442 valid for 
c---- f(p1)+f(p2) --> W/Z(-->l(p3)+l(p4)) + gamma(p5)+f(p6) 
c---- routine frixione.f C Williams 9th Nov 2010

      subroutine basic_Vfrix(pin,passed,isub) 
      implicit none 
      include 'constants.f' 
      include 'frag.f'
      include 'process.f'
      include 'part.f'
      double precision pin(mxpart,4)
      integer isub
      logical passed,in_cone
      double precision ret_ET,R
      double precision R56,ET6,ET5,alpha
      
      passed=.true.

c--- This routine also called for 'Zgajet' -> always need to check
      if ((case .eq. 'Wgamma') .or. (case .eq. 'Zgamma')) then
c---- LO, virtual and dipole pieces have no jets --> passed cut
c---  Real radiation (isub=0) should be checked
        if ( (part .eq. 'lord') .or. (part .eq. 'virt')
     &  .or. (isub .eq. 1) ) return
      endif
     
c---- Calculate angles and energies
         ET5=ret_ET(pin,5)
         ET6=ret_ET(pin,6)
         R56=R(pin,5,6)
         
         alpha=ET5/(1d0-dcos(cone_ang))

c---- Calcualte prefac for check 
                 
         
         if (R56 .lt. cone_ang) then 
            passed=in_cone(R56,ET6,alpha) 
         else
            passed = .true. 
         endif
     
      return 
      end
      

