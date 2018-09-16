      subroutine wamp(mq,qwidth,p1,p2,ie,in,jn,je,jb,amp)
      implicit none
      include 'types.f'
c--- Routine modified 1/6/07 to handle different quark masses & widths
c     g(-p1)+s(-p2) --> e^-(ie)+n(in)+Pn(jn)+Pe^+(je)+b(jb)
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: p1,p2,ie,in,je,jn,jb
      complex(dp):: amp(2)
      integer:: j
      real(dp):: propp,propd,propt,taugt,tsq,mq,qwidth
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      common/zabprods/zab,zba
!$omp threadprivate(/zabprods/)

      propp=sqrt((s(ie,in)-wmass**2)**2+(wmass*wwidth)**2)
      propd=sqrt((s(je,jn)-wmass**2)**2+(wmass*wwidth)**2)
      taugt=s(p1,je)+s(p1,jn)+s(p1,jb)

      tsq  =s(je,jn)+s(je,jb)+s(jn,jb)
      propt=sqrt((tsq-mq**2)**2+(mq*qwidth)**2) 

c--- label on amplitudes represents gluon helicity
c---  amp(1)= negative helicity, amp(2)= positive helicity    
      amp(1) =
     &  + 1/(zb(p1,p2))*za(p1,ie)*za(jn,jb)*zb(p2,in)*zb(p2,je)*mq**2*
     & taugt**(-1)
     &  - 1/(zb(p1,p2))*za(ie,in)*za(jn,jb)*zb(p2,in)**2*zab(p1,je)*
     & taugt**(-1)

      amp(2) =
     &  - 1/(za(p1,p2))*za(p2,ie)*za(jn,jb)*zb(p1,je)*zb(p2,in)*mq**2*
     & taugt**(-1)
     &  + 1/(za(p1,p2))*za(jn,jb)*zb(p1,in)*zab(ie,je)
     &  + 1/(za(p1,p2))*za(jn,jb)*zb(p2,in)*zab(p2,je)*zab(ie,p1)*
     & taugt**(-1)

      do j=1,2
      amp(j)=amp(j)/propp/propd/propt
      enddo

      return
      end 
