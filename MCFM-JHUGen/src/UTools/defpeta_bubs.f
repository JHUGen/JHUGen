!=== C. Williams July 2013
!=== routine which takes in a momentum, P and defines two orthogonal 
!=== spinors p and eta which are can be used in double cuts
      subroutine defpeta_bubs(P,p1,p2) 
      implicit none 
      include 'constants.f'
      double precision P(4),p1(4),p2(4)
      double precision Psq,Pdoteta
      double precision ran2
      double precision phi1,eta1,Et1
      double precision r1,r2,r3,tau
      integer i

      Psq=P(4)**2-P(1)**2-P(2)**2-P(3)**2
      r1=ran2()
      r2=ran2()
      r3=ran2()

!==== first define a random light-like vector 
      Et1=r1*Psq/2d0
      eta1=(-5d0*(one-r2)+5d0*r2)
      phi1=twopi*r3

      p1(1)=Et1*dcos(phi1)
      p1(2)=Et1*dsin(phi1)
      p1(3)=Et1*dsinh(eta1)
      p1(4)=Et1*dcosh(eta1)

      Pdoteta=0d0
      do i=1,3
         Pdoteta=Pdoteta-2d0*(P(i)*p1(i))
      enddo
      Pdoteta=Pdoteta+2d0*(P(4)*p1(4))

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
      include 'constants.f' 
      include 'zprods_decl.f'
      include 'utoolscc.f' 
      double precision p(mxpart,4),p1(4),p2(4),pmo(mxpart,4) 
      double precision inP(4)
      integer n,i,nu
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
