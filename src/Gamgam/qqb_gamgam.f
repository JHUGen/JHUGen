      subroutine qqb_gamgam(p,msq)
C=====
C-----Matrix element for f(-p1)+f(-p2)->gamma(p3)+gamma(p4)
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'noglue.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     .  qa,aq,gg,statfac,facgg,msqgggaga,Qsum
      parameter(statfac=0.5d0)

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
      enddo
      enddo

      call dotem(3,p,s)

      fac=8d0*xn*esq**2*statfac

      Qsum=+Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2 
      facgg=4d0*esq*gsq/(16d0*pisq)*Qsum

      if (omitgg) then
        gg=0d0
      else
        gg=avegg*V*facgg**2*msqgggaga(s(1,2),s(1,3),s(2,3))*statfac
      endif
      
      qa=fac*aveqq*(s(1,3)/s(2,3)+s(2,3)/s(1,3))
      aq=qa

      do j=-nf,nf
      k=-j
C--qa      
      if (j .gt. 0) then
        msq(j,k)=Q(j)**4*qa
C--aq      
      elseif (k .gt. 0) then
        msq(j,k)=Q(k)**4*aq
C--gg
      elseif (j .eq. 0) then
        msq(j,k)=gg
      endif
 
      enddo

      return
      end


      double precision function msqgggaga(s,t,u)
C---amplitudes taken from Eqs. 2,3 of Bern,Dixon,Schmidt
C---hep-ph/0206194
      implicit none
      double precision s,t,u
      double complex m1(2,2,2,2)
      integer h1,h2,h3,h4

      call m1fill(s,t,u,m1)
      msqgggaga=0d0
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      msqgggaga=msqgggaga+cdabs(m1(h1,h2,h3,h4))**2
      enddo
      enddo
      enddo
      enddo
      return
      end

      
      subroutine m1fill(s,t,u,m1)
      implicit none
      include 'constants.f'
      double precision s,t,u
      double complex m1(2,2,2,2)
      integer h1,h2,h3,h4

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
     & *(-0.5d0*(t**2+u**2)/s**2*(log(t/u)**2+pisq)
     & -(t-u)/s*log(t/u)-cone)
      m1(1,2,1,2)=m1(1,2,1,2)
     & *(-0.5d0*(t**2+s**2)/u**2*log(-t/s)**2
     & -(t-s)/u*log(-t/s)-cone
     & -impi*((t**2+s**2)/u**2*log(-t/s)+(t-s)/u))
      m1(2,1,1,2)=m1(2,1,1,2)
     & *(-0.5d0*(u**2+s**2)/t**2*log(-u/s)**2
     & -(u-s)/t*log(-u/s)-cone
     & -impi*((u**2+s**2)/t**2*log(-u/s)+(u-s)/t))

C implement parity
      m1(2,2,1,1)=m1(1,1,2,2)
      m1(2,1,2,1)=m1(1,2,1,2)
      m1(1,2,2,1)=m1(2,1,1,2)

      return
      end
