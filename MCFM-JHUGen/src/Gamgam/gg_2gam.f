      subroutine gg_2gam(p,msq)
      implicit none
      include 'types.f'
C-----Matrix element for g(-p1)+g(-p2)->gamma(p3)+gamma(p4)      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  facgg,msqgggaga,Qsum
      real(dp),parameter::statfac=0.5_dp

c--set msq=0 to initalize
      msq(:,:)=0._dp

      call dotem(3,p,s)

      Qsum=Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2 
      facgg=4._dp*esq*gsq/(16._dp*pisq)*Qsum

      msq(0,0)=avegg*V*facgg**2*msqgggaga(s(1,2),s(1,3),s(2,3))*statfac
      
      return
      end


      function msqgggaga(s,t,u)
      implicit none
      include 'types.f'
      real(dp):: msqgggaga
C---amplitudes taken from Eqs. 2,3 of Bern,Dixon,Schmidt
C---hep-ph/0206194
      
      real(dp):: s,t,u
      complex(dp):: m1(2,2,2,2)
      integer:: h1,h2,h3,h4

      call m1fill(s,t,u,m1)
      msqgggaga=0._dp
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      msqgggaga=msqgggaga+abs(m1(h1,h2,h3,h4))**2
      enddo
      enddo
      enddo
      enddo
      return
      end

      
      subroutine m1fill(s,t,u,m1)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: s,t,u
      complex(dp):: m1(2,2,2,2)
      integer:: h1,h2,h3,h4

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      m1(h1,h2,h3,h4)=cone
      enddo
      enddo
      enddo
      enddo

      m1(1,1,2,2)=m1(1,1,2,2)
     & *(-0.5_dp*(t**2+u**2)/s**2*(log(t/u)**2+pisq)
     & -(t-u)/s*log(t/u)-cone)
      m1(1,2,1,2)=m1(1,2,1,2)
     & *(-0.5_dp*(t**2+s**2)/u**2*log(-t/s)**2
     & -(t-s)/u*log(-t/s)-cone
     & -impi*((t**2+s**2)/u**2*log(-t/s)+(t-s)/u))
      m1(2,1,1,2)=m1(2,1,1,2)
     & *(-0.5_dp*(u**2+s**2)/t**2*log(-u/s)**2
     & -(u-s)/t*log(-u/s)-cone
     & -impi*((u**2+s**2)/t**2*log(-u/s)+(u-s)/t))

C implement parity
      m1(2,2,1,1)=m1(1,1,2,2)
      m1(2,1,2,1)=m1(1,2,1,2)
      m1(1,2,2,1)=m1(2,1,1,2)

      return
      end
