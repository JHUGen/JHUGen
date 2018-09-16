!=== C. Williams July 2013
!=== routine which takes in a momentum, P and defines two orthogonal 
!=== spinors p and eta which are can be used in double cuts
      subroutine defpeta_bubs(P,p1,p2) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: P(4),p1(4),p2(4)
      real(dp):: Psq,Pdoteta
      real(dp):: ran2
      real(dp):: phi1,eta1,Et1
      real(dp):: r1,r2,r3,tau
      integer:: i

      Psq=P(4)**2-P(1)**2-P(2)**2-P(3)**2
      r1=ran2()
      r2=ran2()
      r3=ran2()

!==== first define a random light-like vector 
      Et1=r1*Psq/2._dp
      eta1=(-5._dp*(one-r2)+5._dp*r2)
      phi1=twopi*r3

      p1(1)=Et1*cos(phi1)
      p1(2)=Et1*sin(phi1)
      p1(3)=Et1*sinh(eta1)
      p1(4)=Et1*cosh(eta1)

      Pdoteta=0._dp
      do i=1,3
         Pdoteta=Pdoteta-2._dp*(P(i)*p1(i))
      enddo
      Pdoteta=Pdoteta+2._dp*(P(4)*p1(4))

!====== now rescale that massless vector by Psq/(2 P.p1)
      tau=Psq/Pdoteta
      
      do i=1,4
         p1(i)=tau*p1(i)
         p2(i)=P(i)-p1(i)
      enddo
      
!==== check conditions 
!======
      
!==== write 
c      write(6,*) 'P in ',P(1),P(2),P(3),P(4),'mass',Psq
c      write(6,*) 'p1+p2 ',p1(1)+p2(1),p1(2)+p2(2),p1(3)+p2(3)
c     &,p1(4)+p2(4)
c      write(6,*) 'p1 out ',p1(1),p1(2),p1(3),p1(4),'mass'
c     &,p1(4)**2-p1(3)**2-p1(2)**2-p1(1)**2
c      write(6,*) 'p2 out ',p2(1),p2(2),p2(3),p2(4),'mass' 
c     &,p2(4)**2-p2(3)**2-p2(2)**2-p2(1)**2
     
      return 
      end

      subroutine gen_peta_sprods(p,n,inp,za,zb) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'zprods_decl.f'
      include 'utoolscc.f' 
      real(dp):: p(mxpart,4),p1(4),p2(4),pmo(mxpart,4) 
      real(dp):: inP(4)
      integer:: n,i,nu
!==== n is the number of final state particles 

      write(6,*) 'this routine no longer in use'
      stop
      call defpeta_bubs(inp,p1,p2)
      do i=1,n
         do nu=1,4
            pmo(i,nu)=p(i,nu) 
         enddo
      enddo
      
      do nu=1,4
         pmo(n+1,nu)=p1(nu)
         pmo(n+2,nu)=p2(nu)
      enddo

c         call writeout(pmo)
c         write(6,*) inp 
         
      if (utoolscc) then
        call spinoruorz(n+2,pmo,zb,za) 
      else      
        call spinoruorz(n+2,pmo,za,zb) 
      endif
      
      
      return
      end
