      subroutine qqb_zgam(p,msq)
      implicit none
      include 'types.f'
      
C-----Author Keith Ellis, October 2002
C----- updated: John Campbell, August 2011 (anomalous couplings)
c---- Matrix element for Z/gamma+gamma production
C---- averaged over initial colours and spins
c---
c     q(-p1)+qbar(-p2)-->(e^-(p3)+e^+(p4))+gamma(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'zprods_decl.f'
      integer:: j,k,h12,h34,h5
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,qqb(2),qbq(2)
      complex(dp):: qbqamp(2,2,2,2),qqbamp(2,2,2,2)
      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      
c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call spinoru(5,p,za,zb)

      fac=aveqq*8._dp*esq**3*xn
           
      call zamps(1,2,3,4,5,za,zb,qbqamp)
      call zamps(2,1,3,4,5,za,zb,qqbamp)

      do j=1,2
      qbq(j)=0._dp
      qqb(j)=0._dp
      do h12=1,2
      do h34=1,2
      do h5=1,2
      qbq(j)=qbq(j)+abs(qbqamp(j,h12,h34,h5))**2
      qqb(j)=qqb(j)+abs(qqbamp(j,h12,h34,h5))**2
      enddo
      enddo
      enddo
      qqb(j)=fac*qqb(j)
      qbq(j)=fac*qbq(j)
      enddo

      do j=-nf,nf
      k=-j
          if ((j == 0) .and. (k == 0)) then
            msq(j,k)=0._dp
          elseif ((j > 0) .and. (k < 0)) then
            msq(j,k)=qqb(jj(j))
          elseif ((j < 0) .and. (k > 0)) then
            msq(j,k)=qbq(jj(k))
          endif
      enddo

      return
      end


      subroutine zamps(p1,p2,p3,p4,p5,za,zb,qbqamp)
      implicit none
      include 'types.f'
c--- returns amplitudes including effects of initial state radiation (ai),
c--- final state radiation (af), anomalous ZZ\gamma couplings (anomZZ)
c--- and anomalous Z\gamma\gamma couplings (anomgamZ, anomZgam)
c--- note that af=0 when zerowidth=.true.
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'zerowidth.f'
      include 'anomcoup.f'
      integer:: p1,p2,p3,p4,p5,j
      real(dp):: s12,s34,xfac
      complex(dp)::ai(2,2,2),af(2,2,2),prp34,prp12,qbqamp(2,2,2,2),
     & anomZZ(2,2,2),anomgamZ(2,2,2),anomZgam(2,2,2),cprop
      s12=s(p1,p2)      
      s34=s(p3,p4)      

c--- apply a dipole form factor to anomalous couplings (only if tevscale > 0)
c--- and also rescaling by factors of MZ, c.f. Eqs. (5),(18) of hep-ph/0002138
c--- note that dipole form factors use the third power for h1,h3 (dim. 6)
c--- and the fourth power for h2,h4 (dim. 8)
      if (tevscale > 0._dp) then
        xfac=1._dp/(1._dp+s12/(tevscale*1d3)**2)
      else
        xfac=1._dp
      endif
      
      h1tZ=xfac**3*h1Z/zmass**2
      h2tZ=xfac**4*h2Z/zmass**4
      h3tZ=xfac**3*h3Z/zmass**2
      h4tZ=xfac**4*h4Z/zmass**4
      h1tgam=xfac**3*h1gam/zmass**2
      h2tgam=xfac**4*h2gam/zmass**4
      h3tgam=xfac**3*h3gam/zmass**2
      h4tgam=xfac**4*h4gam/zmass**4

C basic
      ai(1,1,2)=+za(p1,p3)**2*zb(p3,p4)/(za(p1,p5)*za(p2,p5)*s34)
C (3<-->4)
      ai(1,2,2)=+za(p1,p4)**2*zb(p4,p3)/(za(p1,p5)*za(p2,p5)*s34)
C Flip_2: (za<-->zb),(1<-->2),(3<-->4),
      ai(1,1,1)=+zb(p2,p4)**2*za(p4,p3)/(zb(p1,p5)*zb(p2,p5)*s34)
C (za<-->zb),(1<-->2)
      ai(1,2,1)=+zb(p2,p3)**2*za(p3,p4)/(zb(p1,p5)*zb(p2,p5)*s34)

c--- anomalous couplings, hep-ph/0002138 Eq.(19) * (-i)
c--- diagrams qq~ -> Z* -> Z+gam
C basic
      anomZZ(1,1,2)=1._dp/4._dp/s34
     & *((im*h1tZ+h3tZ)*2._dp*za(p1,p3)*zb(p2,p5)*zb(p4,p5)
     &  +(im*h2tZ+h4tZ)*za(p1,p2)*za(p3,p5)*zb(p5,p4)*zb(p2,p5)**2)
C (3<-->4)
      anomZZ(1,2,2)=1._dp/4._dp/s34
     & *((im*h1tZ+h3tZ)*2._dp*za(p1,p4)*zb(p2,p5)*zb(p3,p5)
     &  +(im*h2tZ+h4tZ)*za(p1,p2)*za(p4,p5)*zb(p5,p3)*zb(p2,p5)**2)
C Flip_2: (za<-->zb),(1<-->2),(3<-->4) and h1t -> -h1t, h2t -> -h2t
      anomZZ(1,1,1)=1._dp/4._dp/s34
     & *((-im*h1tZ+h3tZ)*2._dp*zb(p2,p4)*za(p1,p5)*za(p3,p5)
     &  +(-im*h2tZ+h4tZ)*zb(p2,p1)*zb(p4,p5)*za(p5,p3)*za(p1,p5)**2)
C (za<-->zb),(1<-->2)
      anomZZ(1,2,1)=1._dp/4._dp/s34
     & *((-im*h1tZ+h3tZ)*2._dp*zb(p2,p3)*za(p1,p5)*za(p4,p5)
     &  +(-im*h2tZ+h4tZ)*zb(p2,p1)*zb(p3,p5)*za(p5,p4)*za(p1,p5)**2)

c--- diagrams qq~ -> gam* -> Z+gam
C basic
      anomgamZ(1,1,2)=1._dp/4._dp/s34
     & *((im*h1tgam+h3tgam)*2._dp*za(p1,p3)*zb(p2,p5)*zb(p4,p5)
     &  +(im*h2tgam+h4tgam)*za(p1,p2)*za(p3,p5)*zb(p5,p4)*zb(p2,p5)**2)
C (3<-->4)
      anomgamZ(1,2,2)=1._dp/4._dp/s34
     & *((im*h1tgam+h3tgam)*2._dp*za(p1,p4)*zb(p2,p5)*zb(p3,p5)
     &  +(im*h2tgam+h4tgam)*za(p1,p2)*za(p4,p5)*zb(p5,p3)*zb(p2,p5)**2)
C Flip_2: (za<-->zb),(1<-->2),(3<-->4) and h1t -> -h1t, h2t -> -h2t
      anomgamZ(1,1,1)=1._dp/4._dp/s34
     & *((-im*h1tgam+h3tgam)*2._dp*zb(p2,p4)*za(p1,p5)*za(p3,p5)
     &  +(-im*h2tgam+h4tgam)*zb(p2,p1)*zb(p4,p5)*za(p5,p3)*za(p1,p5)**2)
C (za<-->zb),(1<-->2)
      anomgamZ(1,2,1)=1._dp/4._dp/s34
     & *((-im*h1tgam+h3tgam)*2._dp*zb(p2,p3)*za(p1,p5)*za(p4,p5)
     &  +(-im*h2tgam+h4tgam)*zb(p2,p1)*zb(p3,p5)*za(p5,p4)*za(p1,p5)**2)

c--- diagrams qq~ -> Z -> gam*+gam
C basic
      anomZgam(1,1,2)=1._dp/4._dp/s12
     & *((im*h1tZ+h3tZ)*2._dp*za(p3,p1)*zb(p4,p5)*zb(p2,p5)
     &  +(im*h2tZ+h4tZ)*za(p3,p4)*za(p1,p5)*zb(p5,p2)*zb(p4,p5)**2)
C (3<-->4)
      anomZgam(1,2,2)=1._dp/4._dp/s12
     & *((im*h1tZ+h3tZ)*2._dp*za(p4,p1)*zb(p3,p5)*zb(p2,p5)
     &  +(im*h2tZ+h4tZ)*za(p4,p3)*za(p1,p5)*zb(p5,p2)*zb(p3,p5)**2)
C Flip_2: (za<-->zb),(1<-->2),(3<-->4) and h1t -> -h1t, h2t -> -h2t
      anomZgam(1,1,1)=1._dp/4._dp/s12
     & *((-im*h1tZ+h3tZ)*2._dp*zb(p4,p2)*za(p3,p5)*za(p1,p5)
     &  +(-im*h2tZ+h4tZ)*zb(p4,p3)*zb(p2,p5)*za(p5,p1)*za(p3,p5)**2)
C (za<-->zb),(1<-->2)
      anomZgam(1,2,1)=1._dp/4._dp/s12
     & *((-im*h1tZ+h3tZ)*2._dp*zb(p3,p2)*za(p4,p5)*za(p1,p5)
     &  +(-im*h2tZ+h4tZ)*zb(p3,p4)*zb(p2,p5)*za(p5,p1)*za(p4,p5)**2)
      
      if (zerowidth) then

      af(1,1,1)=czip
      af(2,1,1)=czip
      af(1,1,2)=czip
      af(1,2,1)=czip

      else

C basic
      af(1,1,2)=+za(p1,p3)**2*zb(p1,p2)/(za(p3,p5)*za(p4,p5)*s12)
C (1<-->2)
      af(2,1,2)=+za(p2,p3)**2*zb(p2,p1)/(za(p3,p5)*za(p4,p5)*s12)
C Flip_2 (za<-->zb),(1<-->2),(3<-->4)
      af(1,1,1)=+zb(p2,p4)**2*za(p2,p1)/(zb(p3,p5)*zb(p4,p5)*s12)
C Flip_2 (za<-->zb),(3<-->4)
      af(2,1,1)=+zb(p1,p4)**2*za(p1,p2)/(zb(p3,p5)*zb(p4,p5)*s12)

      endif

C     This is complex conjugation
      ai(2,2,1)=conjg(ai(1,1,2))
      ai(2,1,1)=conjg(ai(1,2,2))
      ai(2,1,2)=conjg(ai(1,2,1))
      ai(2,2,2)=conjg(ai(1,1,1))

      af(2,2,1)=conjg(af(1,1,2))
      af(1,2,1)=conjg(af(2,1,2))
      af(2,2,2)=conjg(af(1,1,1))
      af(1,2,2)=conjg(af(2,1,1))

      anomZZ(2,2,1)=conjg(anomZZ(1,1,2))
      anomZZ(2,1,1)=conjg(anomZZ(1,2,2))
      anomZZ(2,1,2)=conjg(anomZZ(1,2,1))
      anomZZ(2,2,2)=conjg(anomZZ(1,1,1))

      anomgamZ(2,2,1)=conjg(anomgamZ(1,1,2))
      anomgamZ(2,1,1)=conjg(anomgamZ(1,2,2))
      anomgamZ(2,1,2)=conjg(anomgamZ(1,2,1))
      anomgamZ(2,2,2)=conjg(anomgamZ(1,1,1))

      anomZgam(2,2,1)=conjg(anomZgam(1,1,2))
      anomZgam(2,1,1)=conjg(anomZgam(1,2,2))
      anomZgam(2,1,2)=conjg(anomZgam(1,2,1))
      anomZgam(2,2,2)=conjg(anomZgam(1,1,1))


c--   calculate propagator factors
      prp34=s34/cplx2((s34-zmass**2),zmass*zwidth)
      prp12=s12/cplx2((s12-zmass**2),zmass*zwidth)
c--   correction factor due to final-state Z not being exactly on-shell
      cprop=(s12-s34)/cplx2((s12-zmass**2),zmass*zwidth)
            
      do j=1,2
      qbqamp(j,1,1,1)=Q(j)*(Q(j)*q1+L(j)*l1*prp34)*ai(1,1,1)
     &                 +q1*(Q(j)*q1+L(j)*l1*prp12)*af(1,1,1)
     &                 +L(j)*l1*prp34*anomZZ(1,1,1)*cprop
     &                 +Q(j)*l1*prp34*anomgamZ(1,1,1)
     &                 +L(j)*q1*prp12*anomZgam(1,1,1)
      qbqamp(j,1,1,2)=Q(j)*(Q(j)*q1+L(j)*l1*prp34)*ai(1,1,2)
     &                 +q1*(Q(j)*q1+L(j)*l1*prp12)*af(1,1,2)
     &                 +L(j)*l1*prp34*anomZZ(1,1,2)*cprop
     &                 +Q(j)*l1*prp34*anomgamZ(1,1,2)
     &                 +L(j)*q1*prp12*anomZgam(1,1,2)

      qbqamp(j,2,2,1)=Q(j)*(Q(j)*q1+R(j)*r1*prp34)*ai(2,2,1)
     &                 +q1*(Q(j)*q1+R(j)*r1*prp12)*af(2,2,1)
     &                 +R(j)*r1*prp34*anomZZ(2,2,1)*cprop
     &                 +Q(j)*r1*prp34*anomgamZ(2,2,1)
     &                 +R(j)*q1*prp12*anomZgam(2,2,1)
      qbqamp(j,2,2,2)=Q(j)*(Q(j)*q1+R(j)*r1*prp34)*ai(2,2,2)
     &                 +q1*(Q(j)*q1+R(j)*r1*prp12)*af(2,2,2)
     &                 +R(j)*r1*prp34*anomZZ(2,2,2)*cprop
     &                 +Q(j)*r1*prp34*anomgamZ(2,2,2)
     &                 +R(j)*q1*prp12*anomZgam(2,2,2)

      qbqamp(j,1,2,1)=Q(j)*(Q(j)*q1+L(j)*r1*prp34)*ai(1,2,1)
     &                 +q1*(Q(j)*q1+L(j)*r1*prp12)*af(1,2,1)
     &                 +L(j)*r1*prp34*anomZZ(1,2,1)*cprop
     &                 +Q(j)*r1*prp34*anomgamZ(1,2,1)
     &                 +L(j)*q1*prp12*anomZgam(1,2,1)
      qbqamp(j,1,2,2)=Q(j)*(Q(j)*q1+L(j)*r1*prp34)*ai(1,2,2)
     &                 +q1*(Q(j)*q1+L(j)*r1*prp12)*af(1,2,2)
     &                 +L(j)*r1*prp34*anomZZ(1,2,2)*cprop
     &                 +Q(j)*r1*prp34*anomgamZ(1,2,2)
     &                 +L(j)*q1*prp12*anomZgam(1,2,2)

      qbqamp(j,2,1,1)=Q(j)*(Q(j)*q1+R(j)*l1*prp34)*ai(2,1,1)
     &                 +q1*(Q(j)*q1+R(j)*l1*prp12)*af(2,1,1)
     &                 +R(j)*l1*prp34*anomZZ(2,1,1)*cprop
     &                 +Q(j)*l1*prp34*anomgamZ(2,1,1)
     &                 +R(j)*q1*prp12*anomZgam(2,1,1)
      qbqamp(j,2,1,2)=Q(j)*(Q(j)*q1+R(j)*l1*prp34)*ai(2,1,2)
     &                 +q1*(Q(j)*q1+R(j)*l1*prp12)*af(2,1,2)
     &                 +R(j)*l1*prp34*anomZZ(2,1,2)*cprop
     &                 +Q(j)*l1*prp34*anomgamZ(2,1,2)
     &                 +R(j)*q1*prp12*anomZgam(2,1,2)

      enddo
      return
      end

