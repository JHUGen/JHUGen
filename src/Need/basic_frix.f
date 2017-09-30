c--- This subroutine is a basic implementaion of hep-ph/9801442 valid for 
c---- f(p1)+f(p2) --> gamma(p3)+f(p4)+f(p5) it should only be used for checking the more general 
c---- routine frixione.f C Williams 9th Nov 2010

      subroutine basic_frix(pin,passed,isub) 
      implicit none 
      include 'constants.f'
      include 'frag.f'
      include 'process.f'
      include 'part.f'
      double precision pin(mxpart,4)
      integer isub
      logical passed,in_cone
      double precision ret_ET,R
      double precision R34,R35,ET3,ET4,ET5,alpha

      passed=.true.
      
     

c--- This routine also called for 'gamjet' -> always need second branch
      if ( ((part .eq. 'lord') .or. (part .eq. 'virt')
     &.or. (isub .eq. 1)) .and. (case .eq. 'dirgam') ) then
c---  LO, virtual and dipole pieces for 'dirgam'       
         ET3=ret_ET(pin,3)
         ET4=ret_ET(pin,4)
                 
         R34=R(pin,3,4)
         
         alpha=ET3/(1d0-dcos(cone_ang))
         
c---- Calculate prefac for check          
         if (R34 .lt. cone_ang) then 
           passed=in_cone(R34,ET4,alpha) 
         else
           passed = .true. 
         endif
         return 


c---   isub=0 real piece for 'dirgam', always for 'gamjet'
      elseif (isub .eq. 0) then 
         
c---- Calculate angles and energies
         ET3=ret_ET(pin,3)
         ET4=ret_ET(pin,4)
         ET5=ret_ET(pin,5)
         
         R34=R(pin,3,4)
         R35=R(pin,3,5)
         
         alpha=ET3/(1d0-dcos(cone_ang)) 
         
         
         if((R34 .lt. cone_ang) .and. (R35 .lt. cone_ang)) then 
c----     2 partons cannot be unresolved at NLO 
            passed = .false. 
c---- This code is not needed
c           if (R34 .lt. R35) then
c             minR=R34
c             maxR=R35
c             minET=ET4
c           else
c             minR=R35
c             maxR=R34
c             minET=ET5
c           endif
c           passed=in_cone(minR,minET,alpha)
c           if (passed) passed=in_cone(maxR,ET4+ET5,alpha)
         elseif (R34 .lt. cone_ang) then 
            passed=in_cone(R34,ET4,alpha) 
         elseif (R35 .lt. cone_ang) then 
            passed=in_cone(R35,ET5,alpha) 
         else
            passed = .true. 
         endif
                        
      endif

      return 
      end
      
      

  
      
      logical function in_cone(Rij,Ejet,pref)
      implicit none
      double precision Rij,Ejet,pref

      if ( Ejet .lt. pref*(1d0-dcos(Rij))) then 
         in_cone = .true. 
      else
         in_cone = .false. 
      endif

      return 
      end
