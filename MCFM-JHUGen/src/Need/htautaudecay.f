      subroutine htautaudecay(p,jm,jp,msq)
************************************************************************
*     Author: J.M. Campbell, August 2012                               *
*                                                                      *
*     matrix element squared for the process of                        *
*     Higgs decay  H --> tau^-(jm)+tau^+(jp)                                 *
*     with bottom mass included                                        *
************************************************************************
      implicit none
      include 'constants.f'
      include 'masses.f'
      integer jm,jp
      double precision p(mxpart,4),s56,msq,msqhtautau
      
      s56=2d0*(p(jm,4)*p(jp,4)-p(jm,1)*p(jp,1)
     &        -p(jm,2)*p(jp,2)-p(jm,3)*p(jp,3))+2d0*mtau**2
      
      msq=msqhtautau(s56)
      
      return
      end
      
      
      
      double precision function msqhtautau(s)
      implicit none
      double precision s
      include 'masses.f'
      include 'ewcouple.f'
 
      msqhtautau=gwsq*mtausq/(4d0*wmass**2)*2d0*(s-4d0*mtau**2)

      return
      end
      
      
