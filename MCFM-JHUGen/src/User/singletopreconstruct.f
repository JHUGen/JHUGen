      subroutine singletopreconstruct(p,failed,costheta,ylight)
      implicit none
      include 'types.f'
c--- Given the usual momentum array, try to reconstruct which
c--- set of momenta reconstruct the top antitop quark;
c--- given those assignments, the routine returns the following quantities:
c---
c---    costheta      the angle between the lepton and the light jet,
c---                    computed in the top rest frame
c---    ylight            the rapidity of the light jet
c---
c--- For real radiation events, the jet algorithm may result in 
c--- events for which the invariant mass is not exactly equal to mt,
c--- despite the phase space producing tops exactly on-shell. In that
c--- case, allow |s-mt| <= 'toler' [GeV]  
c---
c--- If the event does not contain exactly one b-jet and one light jet
c--- then failed is equal to .true.
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'plabel.f'
      logical:: failed
      integer:: i
      real(dp):: tiny,p(mxpart,4),
     & ptop(4),ptoprest(4),plep(4),plight(4),pleprest(4),plightrest(4),
     & mtcand,costheta,ylight,yrap,ptopalt(4),mtcandalt
      parameter (tiny=1.e-8_dp)
c      real(dp):: small,toler
c      parameter (small=1.e-4_dp,toler=5d9)

c--- default: reconstruction is okay
      failed=.false.      
       
c--- note: this routine was written for the case 'notag=1' in chooser.f
c--- (the default is notag=0) and inclusive=F, with the result that
c--- only 2 jets are ever accepted in the event      
      if (p(5,4) < tiny) then
c        write(6,*) 'No b-quark present'
      failed=.true.
      return
      endif

      if (p(7,4) < tiny) then
c        write(6,*) 'No light jet present'
      failed=.true.
      return
      endif

c--- this requires that two jets are b and light-jet
      do i=1,4
        ptop(i)=p(3,i)+p(4,i)+p(5,i)
        ptopalt(i)=p(3,i)+p(4,i)+p(5,i)+p(7,i)
        plight(i)=p(7,i)
      if (plabel(3) == 'el') then
        plep(i)=p(3,i)
      else
        plep(i)=p(4,i)
      endif
      enddo
      mtcand=sqrt(ptop(4)**2-ptop(1)**2-ptop(2)**2-ptop(3)**2)
      mtcandalt=sqrt(ptopalt(4)**2-ptopalt(1)**2
     &               -ptopalt(2)**2-ptopalt(3)**2)
c      write(6,*) 'Event accepted, top candidate mass=',mtcand,mtcandalt

c--- check to see if 345+7 is a better top candidate than 345:
c--- if so, assume light get radiated in decay, so no angle to construct
      if (abs(mtcandalt-mt) < abs(mtcand-mt)) then
c        write(6,*) 'Radiation in decay: ',mtcandalt,' vs ',mtcand
      failed=.true.
      return
      endif
      
      ptoprest(4)=mtcand
      do i=1,3
      ptoprest(i)=0._dp
      enddo
      
      call boostx(plep,ptop,ptoprest,pleprest)
      call boostx(plight,ptop,ptoprest,plightrest)

c--- cos(theta*)
      costheta=(pleprest(1)*plightrest(1)
     &         +pleprest(2)*plightrest(2)
     &         +pleprest(3)*plightrest(3))
     &        /sqrt(pleprest(1)**2+pleprest(2)**2+pleprest(3)**2)
     &        /sqrt(plightrest(1)**2+plightrest(2)**2+plightrest(3)**2)
c      write(6,*) 'cos(theta*)=',costheta
      
c--- rapidity of light jet
      ylight=yrap(7,p)
      
      return
      end
      
