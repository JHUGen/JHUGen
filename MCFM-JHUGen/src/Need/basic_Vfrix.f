c--- This subroutine is a basic implementaion of hep-ph/9801442 valid for 
c---- f(p1)+f(p2) --> W/Z(-->l(p3)+l(p4)) + gamma(p5)+f(p6) 
c---- routine frixione.f C Williams 9th Nov 2010

      subroutine basic_Vfrix(pin,passed,isub) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'frag.f'
      include 'kprocess.f'
      include 'kpart.f'
      real(dp):: pin(mxpart,4)
      integer:: isub
      logical:: passed,in_cone
      real(dp):: ret_ET,R
      real(dp):: R56,ET6,ET5,alpha
      
      passed=.true.

c--- This routine also called for 'Zgajet' -> always need to check
      if ((kcase==kWgamma) .or. (kcase==kZgamma)) then
c---- LO, virtual and dipole pieces have no jets --> passed cut
c---  Real radiation (isub=0) should be checked
        if ( (kpart==klord) .or. (kpart==kvirt)
     &  .or. (isub == 1) ) return
      endif
     
c---- Calculate angles and energies
         ET5=ret_ET(pin,5)
         ET6=ret_ET(pin,6)
         R56=R(pin,5,6)
         
         alpha=ET5/(1._dp-cos(cone_ang))

c---- Calcualte prefac for check 
                 
         
         if (R56 < cone_ang) then 
            passed=in_cone(R56,ET6,alpha) 
         else
            passed = .true. 
         endif
     
      return 
      end
      

