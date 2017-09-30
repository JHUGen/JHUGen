      subroutine wampd(mq,qwidth,p1,pn,pe,pb,t1,amp)
c     t(t1/p1) --> Pn(pn)+Pe^+(pe)+b(pb)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer p1,pe,pn,pb,t1
      double complex amp(2)
      integer i,j
      double precision propd,propt,tsq,mq,qwidth

      propd=dsqrt((s(pe,pn)-wmass**2)**2+(wmass*wwidth)**2)
      tsq  =s(pe,pn)+s(pe,pb)+s(pn,pb)
      propt=dsqrt((tsq-mq**2)**2+(mq*qwidth)**2) 

c      write(6,*) propd,tsq,propt

      do i=1,2
      amp(i)=0d0
      enddo

c--- label on amplitudes represents heavy quark helicity
c---  amp(1)= negative helicity, amp(2)= positive helicity    
      amp(1) =
     &  - za(pn,pb)*zb(pe,t1)

      amp(2) =
     &  - 1/(zb(p1,t1))*za(pn,pb)*zb(p1,pe)*mq

      do j=1,2
      amp(j)=amp(j)/propd/propt
      enddo

      return
      end 
