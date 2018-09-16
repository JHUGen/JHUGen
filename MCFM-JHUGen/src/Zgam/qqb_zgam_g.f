      subroutine qqb_zgam_g(p,msq)
      implicit none
      include 'types.f'
      
C-----Author Keith Ellis, September 2002
C----- updated: John Campbell, August 2011 (anomalous couplings)
c---- Matrix element for Z/gamma+gamma production
C---- averaged over initial colours and spins
c---
c     q(-p1)+qbar(-p2)-->(e^-(p3)+e^+(p4))+gamma(p5)+g(p6)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j,k,h5,h6
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,s345,s34,
     & qqb(2),qbq(2),gq(2),gqb(2),qg(2),qbg(2)
      complex(dp):: pr345,prp34,cprop,
     & qbqi(2,2,2,2),qbqf(2,2,2,2),qqbi(2,2,2,2),qqbf(2,2,2,2),
     & qbgi(2,2,2,2),qbgf(2,2,2,2),qgi(2,2,2,2),qgf(2,2,2,2),
     & gqi(2,2,2,2),gqf(2,2,2,2),gqbi(2,2,2,2),gqbf(2,2,2,2),     
     & qbqZZ(2,2,2,2),qbqgamZ(2,2,2,2),qbqZgam(2,2,2,2),
     & qqbZZ(2,2,2,2),qqbgamZ(2,2,2,2),qqbZgam(2,2,2,2),
     & qbgZZ(2,2,2,2),qbggamZ(2,2,2,2),qbgZgam(2,2,2,2),
     & qgZZ(2,2,2,2),qggamZ(2,2,2,2),qgZgam(2,2,2,2),
     & gqZZ(2,2,2,2),gqgamZ(2,2,2,2),gqZgam(2,2,2,2),
     & gqbZZ(2,2,2,2),gqbgamZ(2,2,2,2),gqbZgam(2,2,2,2)   
      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer,parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      
c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call spinoru(6,p,za,zb)

      s34=s(3,4)
      s345=s(3,4)+s(3,5)+s(4,5)
      
      fac=16._dp*esq**3*xn*cf*gsq

c--   calculate propagator factors
      prp34=s34/cplx2((s34-zmass**2),zmass*zwidth)
      pr345=s345/cplx2((s345-zmass**2),zmass*zwidth)
c--   correction factor due to final-state Z not being exactly on-shell
      cprop=(s345-s34)/cplx2((s345-zmass**2),zmass*zwidth)

      call zamps_g(1,2,3,4,5,6,za,zb,qbqi,qbqf,qbqZZ,qbqgamZ,qbqZgam)
      call zamps_g(2,1,3,4,5,6,za,zb,qqbi,qqbf,qqbZZ,qqbgamZ,qqbZgam)
      call zamps_g(6,2,3,4,5,1,za,zb,gqi,gqf,gqZZ,gqgamZ,gqZgam)
      call zamps_g(2,6,3,4,5,1,za,zb,gqbi,gqbf,gqbZZ,gqbgamZ,gqbZgam)
      call zamps_g(6,1,3,4,5,2,za,zb,qgi,qgf,qgZZ,qggamZ,qgZgam)
      call zamps_g(1,6,3,4,5,2,za,zb,qbgi,qbgf,qbgZZ,qbggamZ,qbgZgam)


      do j=1,2
C---initialize to zero
      qbq(j)=0._dp
      qqb(j)=0._dp
      gq(j)=0._dp
      gqb(j)=0._dp
      qbg(j)=0._dp
      qg(j)=0._dp

c----sum over helicities
      do h5=1,2
      do h6=1,2
      qbq(j)=qbq(j)
     &       +abs(Q(j)*(Q(j)*q1+L(j)*l1*prp34)*qbqi(1,1,h5,h6)
     &               +q1*(Q(j)*q1+L(j)*l1*pr345)*qbqf(1,1,h5,h6)
     &               +L(j)*l1*prp34*qbqZZ(1,1,h5,h6)*cprop
     &               +Q(j)*l1*prp34*qbqgamZ(1,1,h5,h6)
     &               +L(j)*q1*pr345*qbqZgam(1,1,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+R(j)*r1*prp34)*qbqi(2,2,h5,h6)
     &               +q1*(Q(j)*q1+R(j)*r1*pr345)*qbqf(2,2,h5,h6)
     &               +R(j)*r1*prp34*qbqZZ(2,2,h5,h6)*cprop
     &               +Q(j)*r1*prp34*qbqgamZ(2,2,h5,h6)
     &               +R(j)*q1*pr345*qbqZgam(2,2,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+L(j)*r1*prp34)*qbqi(1,2,h5,h6)
     &               +q1*(Q(j)*q1+L(j)*r1*pr345)*qbqf(1,2,h5,h6)
     &               +L(j)*r1*prp34*qbqZZ(1,2,h5,h6)*cprop
     &               +Q(j)*r1*prp34*qbqgamZ(1,2,h5,h6)
     &               +L(j)*q1*pr345*qbqZgam(1,2,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+R(j)*l1*prp34)*qbqi(2,1,h5,h6)
     &               +q1*(Q(j)*q1+R(j)*l1*pr345)*qbqf(2,1,h5,h6)
     &               +R(j)*l1*prp34*qbqZZ(2,1,h5,h6)*cprop
     &               +Q(j)*l1*prp34*qbqgamZ(2,1,h5,h6)
     &               +R(j)*q1*pr345*qbqZgam(2,1,h5,h6))**2

      qqb(j)=qqb(j)
     &       +abs(Q(j)*(Q(j)*q1+L(j)*l1*prp34)*qqbi(1,1,h5,h6)
     &               +q1*(Q(j)*q1+L(j)*l1*pr345)*qqbf(1,1,h5,h6)
     &               +L(j)*l1*prp34*qqbZZ(1,1,h5,h6)*cprop
     &               +Q(j)*l1*prp34*qqbgamZ(1,1,h5,h6)
     &               +L(j)*q1*pr345*qqbZgam(1,1,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+R(j)*r1*prp34)*qqbi(2,2,h5,h6)
     &               +q1*(Q(j)*q1+R(j)*r1*pr345)*qqbf(2,2,h5,h6)
     &               +R(j)*r1*prp34*qqbZZ(2,2,h5,h6)*cprop
     &               +Q(j)*r1*prp34*qqbgamZ(2,2,h5,h6)
     &               +R(j)*q1*pr345*qqbZgam(2,2,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+L(j)*r1*prp34)*qqbi(1,2,h5,h6)
     &               +q1*(Q(j)*q1+L(j)*r1*pr345)*qqbf(1,2,h5,h6)
     &               +L(j)*r1*prp34*qqbZZ(1,2,h5,h6)*cprop
     &               +Q(j)*r1*prp34*qqbgamZ(1,2,h5,h6)
     &               +L(j)*q1*pr345*qqbZgam(1,2,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+R(j)*l1*prp34)*qqbi(2,1,h5,h6)
     &               +q1*(Q(j)*q1+R(j)*l1*pr345)*qqbf(2,1,h5,h6)
     &               +R(j)*l1*prp34*qqbZZ(2,1,h5,h6)*cprop
     &               +Q(j)*l1*prp34*qqbgamZ(2,1,h5,h6)
     &               +R(j)*q1*pr345*qqbZgam(2,1,h5,h6))**2

      gqb(j)=gqb(j)
     &       +abs(Q(j)*(Q(j)*q1+L(j)*l1*prp34)*gqbi(1,1,h5,h6)
     &               +q1*(Q(j)*q1+L(j)*l1*pr345)*gqbf(1,1,h5,h6)
     &               +L(j)*l1*prp34*gqbZZ(1,1,h5,h6)*cprop
     &               +Q(j)*l1*prp34*gqbgamZ(1,1,h5,h6)
     &               +L(j)*q1*pr345*gqbZgam(1,1,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+R(j)*r1*prp34)*gqbi(2,2,h5,h6)
     &               +q1*(Q(j)*q1+R(j)*r1*pr345)*gqbf(2,2,h5,h6)
     &               +R(j)*r1*prp34*gqbZZ(2,2,h5,h6)*cprop
     &               +Q(j)*r1*prp34*gqbgamZ(2,2,h5,h6)
     &               +R(j)*q1*pr345*gqbZgam(2,2,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+L(j)*r1*prp34)*gqbi(1,2,h5,h6)
     &               +q1*(Q(j)*q1+L(j)*r1*pr345)*gqbf(1,2,h5,h6)
     &               +L(j)*r1*prp34*gqbZZ(1,2,h5,h6)*cprop
     &               +Q(j)*r1*prp34*gqbgamZ(1,2,h5,h6)
     &               +L(j)*q1*pr345*gqbZgam(1,2,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+R(j)*l1*prp34)*gqbi(2,1,h5,h6)
     &               +q1*(Q(j)*q1+R(j)*l1*pr345)*gqbf(2,1,h5,h6)
     &               +R(j)*l1*prp34*gqbZZ(2,1,h5,h6)*cprop
     &               +Q(j)*l1*prp34*gqbgamZ(2,1,h5,h6)
     &               +R(j)*q1*pr345*gqbZgam(2,1,h5,h6))**2

      gq(j)=gq(j)
     &       +abs(Q(j)*(Q(j)*q1+L(j)*l1*prp34)*gqi(1,1,h5,h6)
     &               +q1*(Q(j)*q1+L(j)*l1*pr345)*gqf(1,1,h5,h6)
     &               +L(j)*l1*prp34*gqZZ(1,1,h5,h6)*cprop
     &               +Q(j)*l1*prp34*gqgamZ(1,1,h5,h6)
     &               +L(j)*q1*pr345*gqZgam(1,1,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+R(j)*r1*prp34)*gqi(2,2,h5,h6)
     &               +q1*(Q(j)*q1+R(j)*r1*pr345)*gqf(2,2,h5,h6)
     &               +R(j)*r1*prp34*gqZZ(2,2,h5,h6)*cprop
     &               +Q(j)*r1*prp34*gqgamZ(2,2,h5,h6)
     &               +R(j)*q1*pr345*gqZgam(2,2,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+L(j)*r1*prp34)*gqi(1,2,h5,h6)
     &               +q1*(Q(j)*q1+L(j)*r1*pr345)*gqf(1,2,h5,h6)
     &               +L(j)*r1*prp34*gqZZ(1,2,h5,h6)*cprop
     &               +Q(j)*r1*prp34*gqgamZ(1,2,h5,h6)
     &               +L(j)*q1*pr345*gqZgam(1,2,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+R(j)*l1*prp34)*gqi(2,1,h5,h6)
     &               +q1*(Q(j)*q1+R(j)*l1*pr345)*gqf(2,1,h5,h6)
     &               +R(j)*l1*prp34*gqZZ(2,1,h5,h6)*cprop
     &               +Q(j)*l1*prp34*gqgamZ(2,1,h5,h6)
     &               +R(j)*q1*pr345*gqZgam(2,1,h5,h6))**2

      qbg(j)=qbg(j)
     &       +abs(Q(j)*(Q(j)*q1+L(j)*l1*prp34)*qbgi(1,1,h5,h6)
     &               +q1*(Q(j)*q1+L(j)*l1*pr345)*qbgf(1,1,h5,h6)
     &               +L(j)*l1*prp34*qbgZZ(1,1,h5,h6)*cprop
     &               +Q(j)*l1*prp34*qbggamZ(1,1,h5,h6)
     &               +L(j)*q1*pr345*qbgZgam(1,1,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+R(j)*r1*prp34)*qbgi(2,2,h5,h6)
     &               +q1*(Q(j)*q1+R(j)*r1*pr345)*qbgf(2,2,h5,h6)
     &               +R(j)*r1*prp34*qbgZZ(2,2,h5,h6)*cprop
     &               +Q(j)*r1*prp34*qbggamZ(2,2,h5,h6)
     &               +R(j)*q1*pr345*qbgZgam(2,2,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+L(j)*r1*prp34)*qbgi(1,2,h5,h6)
     &               +q1*(Q(j)*q1+L(j)*r1*pr345)*qbgf(1,2,h5,h6)
     &               +L(j)*r1*prp34*qbgZZ(1,2,h5,h6)*cprop
     &               +Q(j)*r1*prp34*qbggamZ(1,2,h5,h6)
     &               +L(j)*q1*pr345*qbgZgam(1,2,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+R(j)*l1*prp34)*qbgi(2,1,h5,h6)
     &               +q1*(Q(j)*q1+R(j)*l1*pr345)*qbgf(2,1,h5,h6)
     &               +R(j)*l1*prp34*qbgZZ(2,1,h5,h6)*cprop
     &               +Q(j)*l1*prp34*qbggamZ(2,1,h5,h6)
     &               +R(j)*q1*pr345*qbgZgam(2,1,h5,h6))**2

      qg(j)=qg(j)
     &       +abs(Q(j)*(Q(j)*q1+L(j)*l1*prp34)*qgi(1,1,h5,h6)
     &               +q1*(Q(j)*q1+L(j)*l1*pr345)*qgf(1,1,h5,h6)
     &               +L(j)*l1*prp34*qgZZ(1,1,h5,h6)*cprop
     &               +Q(j)*l1*prp34*qggamZ(1,1,h5,h6)
     &               +L(j)*q1*pr345*qgZgam(1,1,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+R(j)*r1*prp34)*qgi(2,2,h5,h6)
     &               +q1*(Q(j)*q1+R(j)*r1*pr345)*qgf(2,2,h5,h6)
     &               +R(j)*r1*prp34*qgZZ(2,2,h5,h6)*cprop
     &               +Q(j)*r1*prp34*qggamZ(2,2,h5,h6)
     &               +R(j)*q1*pr345*qgZgam(2,2,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+L(j)*r1*prp34)*qgi(1,2,h5,h6)
     &               +q1*(Q(j)*q1+L(j)*r1*pr345)*qgf(1,2,h5,h6)
     &               +L(j)*r1*prp34*qgZZ(1,2,h5,h6)*cprop
     &               +Q(j)*r1*prp34*qggamZ(1,2,h5,h6)
     &               +L(j)*q1*pr345*qgZgam(1,2,h5,h6))**2
     &       +abs(Q(j)*(Q(j)*q1+R(j)*l1*prp34)*qgi(2,1,h5,h6)
     &               +q1*(Q(j)*q1+R(j)*l1*pr345)*qgf(2,1,h5,h6)
     &               +R(j)*l1*prp34*qgZZ(2,1,h5,h6)*cprop
     &               +Q(j)*l1*prp34*qggamZ(2,1,h5,h6)
     &               +R(j)*q1*pr345*qgZgam(2,1,h5,h6))**2

      enddo
      enddo
      enddo

      do j=-nf,nf
      do k=-nf,nf
          if ((j == 0) .and. (k == 0)) then
            msq(j,k)=0._dp
          elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=aveqg*fac*gqb(-kk(k))
          elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=aveqg*fac*gq(kk(k))
          elseif ((j > 0) .and. (k == -j)) then
            msq(j,k)=aveqq*fac*qqb(jj(j))
          elseif ((j < 0) .and. (k == -j)) then
            msq(j,k)=aveqq*fac*qbq(kk(k))
          elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=aveqg*fac*qg(jj(j))
          elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=aveqg*fac*qbg(-jj(j))
          else 
            msq(j,k)=0._dp
          endif
      enddo
      enddo

      return
      end


      subroutine zamps_g(p1,p2,p3,p4,p5,p6,za,zb,ai,af,
     & anomZZ,anomgamZ,anomZgam)
      implicit none
      include 'types.f'
c--- returns amplitudes corresponding to initial state radiation (ai),
c--- final state radiation (af), anomalous ZZ\gamma couplings (anomZZ)
c--- and anomalous Z\gamma\gamma couplings (anomgamZ, anomZgam)
c--- note that af=0 when zerowidth=.true.
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zerowidth.f'
      include 'masses.f'
      include 'anomcoup.f'
      integer:: p1,p2,p3,p4,p5,p6,h12,h34,h5,h6
      real(dp):: s126,s34,s345,s156,s256,xfac
      complex(dp):: ai(2,2,2,2),af(2,2,2,2),
     & anomZZ(2,2,2,2),anomgamZ(2,2,2,2),anomZgam(2,2,2,2),zazb,zbza
      complex(dp):: z6354,z2354,z6453,z3162,z1264
      complex(dp):: z3456,z1263,z4162,z1254,z6254,z3165,z1253
      complex(dp):: z6253,z4165,z3152,z3156,z5264,z4152,z4156,z5263

      zazb(p1,p2,p3,p4)=+za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zbza(p1,p2,p3,p4)=+zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)

      s126=s(p1,p2)+s(p1,p6)+s(p2,p6)      
      s34=s(p3,p4)      
      s345=s(p3,p4)+s(p3,p5)+s(p4,p5)      
      s156=s(p1,p5)+s(p1,p6)+s(p5,p6)      
      s256=s(p2,p5)+s(p2,p6)+s(p5,p6)      

      z1263=zazb(p1,p2,p6,p3)
      z1253=zazb(p1,p2,p5,p3)

      z1264=zazb(p1,p2,p6,p4)
      z1254=zazb(p1,p2,p5,p4)

      z2354=zazb(p2,p3,p5,p4)

      z3152=zazb(p3,p1,p5,p2)
      z3162=zazb(p3,p1,p6,p2)

      z3156=zazb(p3,p1,p5,p6)
      z3456=zazb(p3,p4,p5,p6)

      z3165=zazb(p3,p1,p6,p5)

      z4162=zazb(p4,p1,p6,p2)
      z4152=zazb(p4,p1,p5,p2)

      z4165=zazb(p4,p1,p6,p5)
      z4156=zazb(p4,p1,p5,p6)

      z5264=zazb(p5,p2,p6,p4)
      z5263=zazb(p5,p2,p6,p3)
 
      z6354=zazb(p6,p3,p5,p4)
      z6254=zazb(p6,p2,p5,p4)

      z6453=zazb(p6,p4,p5,p3)
      z6253=zazb(p6,p2,p5,p3)

c--- apply a dipole form factor to anomalous couplings (only if tevscale > 0)
c--- and also rescaling by factors of MZ, c.f. Eqs. (5),(18) of hep-ph/0002138
c--- note that dipole form factors use the third power for h1,h3 (dim. 6)
c--- and the fourth power for h2,h4 (dim. 8)
      if (tevscale > 0._dp) then
        xfac=1._dp/(1._dp+s126/(tevscale*1d3)**2)
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

C---  c.f. hep/9803250 Eqs (4.9) - (4.12) multiplied by -i
      ai(1,1,1,1)=zb(p1,p2)*zb(p2,p4)**2
     & /(zb(p4,p3)*zb(p1,p5)*zb(p2,p5)*zb(p1,p6)*zb(p2,p6))

      ai(1,2,1,1)=zb(p1,p2)*zb(p2,p3)**2
     & /(zb(p3,p4)*zb(p1,p5)*zb(p2,p5)*zb(p1,p6)*zb(p2,p6))

      ai(1,1,2,2)=za(p1,p2)*za(p1,p3)**2
     & /(za(p4,p3)*za(p1,p5)*za(p2,p5)*za(p1,p6)*za(p2,p6))

      ai(2,1,1,1)=zb(p1,p2)*zb(p1,p4)**2
     & /(zb(p3,p4)*zb(p1,p5)*zb(p2,p5)*zb(p1,p6)*zb(p2,p6))

      ai(1,1,2,1)=-z1254*z3162
     & /(s34*za(p1,p5)*za(p2,p5)*zb(p1,p6)*zb(p2,p6))
     & +za(p1,p3)*zb(p2,p5)*z6254/(s34*s256*za(p2,p5)*zb(p2,p6))
     & +za(p1,p6)*zb(p2,p4)*z3165/(s34*s156*za(p1,p5)*zb(p1,p6))

      ai(1,2,2,1)=-z1253*z4162
     & /(s34*za(p1,p5)*za(p2,p5)*zb(p1,p6)*zb(p2,p6))
     & +za(p1,p4)*zb(p2,p5)*z6253/(s34*s256*za(p2,p5)*zb(p2,p6))
     & +za(p1,p6)*zb(p2,p3)*z4165/(s34*s156*za(p1,p5)*zb(p1,p6))

      ai(1,1,1,2)=-z1264*z3152
     & /(s34*za(p1,p6)*za(p2,p6)*zb(p1,p5)*zb(p2,p5))
     & +za(p1,p5)*zb(p2,p4)*z3156/(s34*s156*za(p1,p6)*zb(p1,p5))
     & +za(p1,p3)*zb(p2,p6)*z5264/(s34*s256*za(p2,p6)*zb(p2,p5))

      ai(1,2,1,2)=-z1263*z4152
     & /(s34*za(p1,p6)*za(p2,p6)*zb(p1,p5)*zb(p2,p5))
     & +za(p1,p5)*zb(p2,p3)*z4156/(s34*s156*za(p1,p6)*zb(p1,p5))
     & +za(p1,p4)*zb(p2,p6)*z5263/(s34*s256*za(p2,p6)*zb(p2,p5))

      ai(1,2,2,2)=conjg(ai(2,1,1,1))
      ai(2,1,2,2)=conjg(ai(1,2,1,1))
      ai(2,2,1,2)=conjg(ai(1,1,2,1))
      ai(2,2,2,1)=conjg(ai(1,1,1,2))

      ai(2,2,2,2)=conjg(ai(1,1,1,1))
      ai(2,2,1,1)=conjg(ai(1,1,2,2))
      ai(2,1,2,1)=conjg(ai(1,2,1,2))
      ai(2,1,1,2)=conjg(ai(1,2,2,1))

      if (zerowidth) then

      do h12=1,2
      do h34=1,2
      do h5=1,2
      do h6=1,2
      af(h12,h34,h5,h6)=czip
      enddo
      enddo
      enddo
      enddo

      else

      af(1,1,1,1)=+(z6354*zb(p2,p6)+z1264*zb(p1,p2))*zb(p2,p4)
     & /(s345*zb(p1,p6)*zb(p2,p6)*zb(p3,p5)*zb(p4,p5))
      af(1,2,1,1)=-(z6453*zb(p2,p6)+z1263*zb(p1,p2))*zb(p2,p3)
     & /(s345*zb(p1,p6)*zb(p2,p6)*zb(p3,p5)*zb(p4,p5))
      af(2,1,1,1)=-(z2354*zb(p1,p2)+z6354*zb(p1,p6))*zb(p1,p4)
     & /(s345*zb(p1,p6)*zb(p2,p6)*zb(p3,p5)*zb(p4,p5))
      af(1,1,2,2)=+(z3456*za(p1,p6)-z3162*za(p1,p2))*za(p1,p3)
     & /(s345*za(p1,p6)*za(p2,p6)*za(p3,p5)*za(p4,p5))

      af(1,1,1,2)=-z1264**2 
     & /(s345*za(p1,p6)*za(p2,p6)*zb(p3,p5)*zb(p4,p5))
      af(1,2,1,2)=+z1263**2 
     & /(s345*za(p1,p6)*za(p2,p6)*zb(p3,p5)*zb(p4,p5))
      af(1,1,2,1)=-z3162**2
     & /(s345*za(p3,p5)*za(p4,p5)*zb(p1,p6)*zb(p2,p6))
      af(1,2,2,1)=+z4162**2  
     & /(s345*za(p3,p5)*za(p4,p5)*zb(p1,p6)*zb(p2,p6))

      af(1,2,2,2)=conjg(af(2,1,1,1))
      af(2,1,2,2)=conjg(af(1,2,1,1))
      af(2,2,1,2)=conjg(af(1,1,2,1))
      af(2,2,2,1)=conjg(af(1,1,1,2))

      af(2,2,2,2)=conjg(af(1,1,1,1))
      af(2,2,1,1)=conjg(af(1,1,2,2))
      af(2,1,2,1)=conjg(af(1,2,1,2))
      af(2,1,1,2)=conjg(af(1,2,2,1))

      endif

c--- anomalous couplings, hep-ph/0002138 Eq.(21) * (-i)
c--- diagrams qq~ -> Z* -> Z+gam
c--- basic
      anomZZ(1,1,2,2)=1._dp/(4._dp*s34*za(p1,p6)*za(p2,p6))
     & *((im*h1tZ+h3tZ)*2._dp*za(p1,p3)*zb(p4,p5)*zazb(p1,p2,p6,p5)
     &  +(im*h2tZ+h4tZ)*za(p3,p5)*zb(p5,p4)*zazb(p1,p2,p6,p5)**2)
c--- anomalous couplings, hep-ph/0002138 Eq.(22) * (-i)
      anomZZ(1,1,2,1)=1._dp/(4._dp*s34*zb(p1,p6)*zb(p2,p6))
     & *((im*h1tZ+h3tZ)*2._dp*zb(p4,p5)*zb(p2,p5)*zazb(p3,p1,p6,p2)
     &  +(im*h2tZ+h4tZ)*s345*za(p3,p5)*zb(p5,p4)*zb(p2,p5)**2)
C (3<-->4) ! but manually perform symmetry of (1,2,2,2) to (2,1,1,1)
           ! in order to use the same relations as for the SM amplitudes
c      anomZZ(1,2,2,2)=1._dp/(4._dp*s34*za(p1,p6)*za(p2,p6))
c     & *((im*h1tZ+h3tZ)*2._dp*za(p1,p4)*zb(p3,p5)*zazb(p1,p2,p6,p5)
c     &  +(im*h2tZ+h4tZ)*za(p4,p5)*zb(p5,p3)*zazb(p1,p2,p6,p5)**2)
      anomZZ(2,1,1,1)=1._dp/(4._dp*s34*zb(p1,p6)*zb(p2,p6))
     & *((-im*h1tZ+h3tZ)*2._dp*zb(p1,p4)*za(p3,p5)*zbza(p1,p2,p6,p5)
     &  +(-im*h2tZ+h4tZ)*zb(p4,p5)*za(p5,p3)*zbza(p1,p2,p6,p5)**2)
      anomZZ(1,2,2,1)=1._dp/(4._dp*s34*zb(p1,p6)*zb(p2,p6))
     & *((im*h1tZ+h3tZ)*2._dp*zb(p3,p5)*zb(p2,p5)*zazb(p4,p1,p6,p2)
     &  +(im*h2tZ+h4tZ)*s345*za(p4,p5)*zb(p5,p3)*zb(p2,p5)**2)
C Flip_2: (za<-->zb),(1<-->2),(3<-->4) and h1t -> -h1t, h2t -> -h2t
      anomZZ(1,1,1,1)=1._dp/(4._dp*s34*zb(p2,p6)*zb(p1,p6))
     & *((-im*h1tZ+h3tZ)*2._dp*zb(p2,p4)*za(p3,p5)*zbza(p2,p1,p6,p5)
     &  +(-im*h2tZ+h4tZ)*zb(p4,p5)*za(p5,p3)*zbza(p2,p1,p6,p5)**2)
      anomZZ(1,1,1,2)=1._dp/(4._dp*s34*za(p2,p6)*za(p1,p6))
     & *((-im*h1tZ+h3tZ)*2._dp*za(p3,p5)*za(p1,p5)*zbza(p4,p2,p6,p1)
     &  +(-im*h2tZ+h4tZ)*s345*zb(p4,p5)*za(p5,p3)*za(p1,p5)**2)
C (za<-->zb),(1<-->2) and h1t -> -h1t, h2t -> -h2t
      anomZZ(1,2,1,1)=1._dp/(4._dp*s34*zb(p2,p6)*zb(p1,p6))
     & *((-im*h1tZ+h3tZ)*2._dp*zb(p2,p3)*za(p4,p5)*zbza(p2,p1,p6,p5)
     &  +(-im*h2tZ+h4tZ)*zb(p3,p5)*za(p5,p4)*zbza(p2,p1,p6,p5)**2)
      anomZZ(1,2,1,2)=1._dp/(4._dp*s34*za(p2,p6)*za(p1,p6))
     & *((-im*h1tZ+h3tZ)*2._dp*za(p4,p5)*za(p1,p5)*zbza(p3,p2,p6,p1)
     &  +(-im*h2tZ+h4tZ)*s345*zb(p3,p5)*za(p5,p4)*za(p1,p5)**2)
      
      anomZZ(1,2,2,2)=conjg(anomZZ(2,1,1,1))
      anomZZ(2,1,2,2)=conjg(anomZZ(1,2,1,1))
      anomZZ(2,2,1,2)=conjg(anomZZ(1,1,2,1))
      anomZZ(2,2,2,1)=conjg(anomZZ(1,1,1,2))

      anomZZ(2,2,2,2)=conjg(anomZZ(1,1,1,1))
      anomZZ(2,2,1,1)=conjg(anomZZ(1,1,2,2))
      anomZZ(2,1,2,1)=conjg(anomZZ(1,2,1,2))
      anomZZ(2,1,1,2)=conjg(anomZZ(1,2,2,1))

c--- diagrams qq~ -> gam* -> Z+gam
c--- basic
      anomgamZ(1,1,2,2)=1._dp/(4._dp*s34*za(p1,p6)*za(p2,p6))
     & *((im*h1tgam+h3tgam)*2._dp*za(p1,p3)*zb(p4,p5)*zazb(p1,p2,p6,p5)
     &  +(im*h2tgam+h4tgam)*za(p3,p5)*zb(p5,p4)*zazb(p1,p2,p6,p5)**2)
c--- anomalous couplings, hep-ph/0002138 Eq.(22) * (-i)
      anomgamZ(1,1,2,1)=1._dp/(4._dp*s34*zb(p1,p6)*zb(p2,p6))
     & *((im*h1tgam+h3tgam)*2._dp*zb(p4,p5)*zb(p2,p5)*zazb(p3,p1,p6,p2)
     &  +(im*h2tgam+h4tgam)*s345*za(p3,p5)*zb(p5,p4)*zb(p2,p5)**2)
C (3<-->4) ! but manually perform symmetry of (1,2,2,2) to (2,1,1,1)
           ! in order to use the same relations as for the SM amplitudes
c      anomgamZ(1,2,2,2)=1._dp/(4._dp*s34*za(p1,p6)*za(p2,p6))
c     & *((im*h1tgam+h3tgam)*2._dp*za(p1,p4)*zb(p3,p5)*zazb(p1,p2,p6,p5)
c     &  +(im*h2tgam+h4tgam)*za(p4,p5)*zb(p5,p3)*zazb(p1,p2,p6,p5)**2)
      anomgamZ(2,1,1,1)=1._dp/(4._dp*s34*zb(p1,p6)*zb(p2,p6))
     & *((-im*h1tgam+h3tgam)*2._dp*zb(p1,p4)*za(p3,p5)*zbza(p1,p2,p6,p5)
     &  +(-im*h2tgam+h4tgam)*zb(p4,p5)*za(p5,p3)*zbza(p1,p2,p6,p5)**2)
      anomgamZ(1,2,2,1)=1._dp/(4._dp*s34*zb(p1,p6)*zb(p2,p6))
     & *((im*h1tgam+h3tgam)*2._dp*zb(p3,p5)*zb(p2,p5)*zazb(p4,p1,p6,p2)
     &  +(im*h2tgam+h4tgam)*s345*za(p4,p5)*zb(p5,p3)*zb(p2,p5)**2)
C Flip_2: (za<-->zb),(1<-->2),(3<-->4) and h1t -> -h1t, h2t -> -h2t
      anomgamZ(1,1,1,1)=1._dp/(4._dp*s34*zb(p2,p6)*zb(p1,p6))
     & *((-im*h1tgam+h3tgam)*2._dp*zb(p2,p4)*za(p3,p5)*zbza(p2,p1,p6,p5)
     &  +(-im*h2tgam+h4tgam)*zb(p4,p5)*za(p5,p3)*zbza(p2,p1,p6,p5)**2)
      anomgamZ(1,1,1,2)=1._dp/(4._dp*s34*za(p2,p6)*za(p1,p6))
     & *((-im*h1tgam+h3tgam)*2._dp*za(p3,p5)*za(p1,p5)*zbza(p4,p2,p6,p1)
     &  +(-im*h2tgam+h4tgam)*s345*zb(p4,p5)*za(p5,p3)*za(p1,p5)**2)
C (za<-->zb),(1<-->2) and h1t -> -h1t, h2t -> -h2t
      anomgamZ(1,2,1,1)=1._dp/(4._dp*s34*zb(p2,p6)*zb(p1,p6))
     & *((-im*h1tgam+h3tgam)*2._dp*zb(p2,p3)*za(p4,p5)*zbza(p2,p1,p6,p5)
     &  +(-im*h2tgam+h4tgam)*zb(p3,p5)*za(p5,p4)*zbza(p2,p1,p6,p5)**2)
      anomgamZ(1,2,1,2)=1._dp/(4._dp*s34*za(p2,p6)*za(p1,p6))
     & *((-im*h1tgam+h3tgam)*2._dp*za(p4,p5)*za(p1,p5)*zbza(p3,p2,p6,p1)
     &  +(-im*h2tgam+h4tgam)*s345*zb(p3,p5)*za(p5,p4)*za(p1,p5)**2)
      
      anomgamZ(1,2,2,2)=conjg(anomgamZ(2,1,1,1))
      anomgamZ(2,1,2,2)=conjg(anomgamZ(1,2,1,1))
      anomgamZ(2,2,1,2)=conjg(anomgamZ(1,1,2,1))
      anomgamZ(2,2,2,1)=conjg(anomgamZ(1,1,1,2))

      anomgamZ(2,2,2,2)=conjg(anomgamZ(1,1,1,1))
      anomgamZ(2,2,1,1)=conjg(anomgamZ(1,1,2,2))
      anomgamZ(2,1,2,1)=conjg(anomgamZ(1,2,1,2))
      anomgamZ(2,1,1,2)=conjg(anomgamZ(1,2,2,1))

c--- diagrams qq~ -> Z -> gam*+gam
C basic
      anomZgam(1,1,2,2)=-1._dp/(4._dp*s126*za(p1,p6)*za(p2,p6))
     & *zb(p4,p5)*(za(p1,p2)*zb(p2,p5)+za(p1,p6)*zb(p6,p5))
     & *((im*h1tgam+h3tgam)*2._dp*za(p1,p3)
     &  +(im*h2tgam+h4tgam)*za(p3,p4)*zb(p4,p5)*za(p1,p5))
      anomZgam(1,1,2,1)=-1._dp/(4._dp*s126*zb(p1,p6)*zb(p2,p6))
     & *zb(p4,p5)*zb(p2,p5)
     & *((im*h1tgam+h3tgam)
     &   *2._dp*(za(p3,p1)*zb(p1,p2)+za(p3,p6)*zb(p6,p2))
     &  +(im*h2tgam+h4tgam)
     &   *za(p3,p4)*zb(p4,p5)*(za(p5,p1)*zb(p1,p2)+za(p5,p6)*zb(p6,p2)))
C (3<-->4) ! but manually perform symmetry of (1,2,2,2) to (2,1,1,1)
           ! in order to use the same relations as for the SM amplitudes
      anomZgam(2,1,1,1)=-1._dp/(4._dp*s126*zb(p1,p6)*zb(p2,p6))
     & *za(p3,p5)*(zb(p1,p2)*za(p2,p5)+zb(p1,p6)*za(p6,p5))
     & *((-im*h1tgam+h3tgam)*2._dp*zb(p1,p4)
     &  +(-im*h2tgam+h4tgam)*zb(p4,p3)*za(p3,p5)*zb(p1,p5))
      anomZgam(1,2,2,1)=-1._dp/(4._dp*s126*zb(p1,p6)*zb(p2,p6))
     & *zb(p3,p5)*zb(p2,p5)
     & *((im*h1tgam+h3tgam)
     &   *2._dp*(za(p4,p1)*zb(p1,p2)+za(p4,p6)*zb(p6,p2))
     &  +(im*h2tgam+h4tgam)
     &   *za(p4,p3)*zb(p3,p5)*(za(p5,p1)*zb(p1,p2)+za(p5,p6)*zb(p6,p2)))
C Flip_2: (za<-->zb),(1<-->2),(3<-->4) and h1t -> -h1t, h2t -> -h2t
      anomZgam(1,1,1,1)=-1._dp/(4._dp*s126*zb(p1,p6)*zb(p2,p6))
     & *za(p3,p5)*(zb(p2,p1)*za(p1,p5)+zb(p2,p6)*za(p6,p5))
     & *((-im*h1tgam+h3tgam)*2._dp*zb(p2,p4)
     &  +(-im*h2tgam+h4tgam)*zb(p4,p3)*za(p3,p5)*zb(p2,p5))
      anomZgam(1,1,1,2)=-1._dp/(4._dp*s126*za(p1,p6)*za(p2,p6))
     & *za(p3,p5)*za(p1,p5)
     & *((-im*h1tgam+h3tgam)
     &   *2._dp*(zb(p4,p2)*za(p2,p1)+zb(p4,p6)*za(p6,p1))
     &  +(-im*h2tgam+h4tgam)
     &   *zb(p4,p3)*za(p3,p5)*(zb(p5,p2)*za(p2,p1)+zb(p5,p6)*za(p6,p1)))
C (za<-->zb),(1<-->2) and h1t -> -h1t, h2t -> -h2t
      anomZgam(1,2,1,1)=-1._dp/(4._dp*s126*zb(p1,p6)*zb(p2,p6))
     & *za(p4,p5)*(zb(p2,p1)*za(p1,p5)+zb(p2,p6)*za(p6,p5))
     & *((-im*h1tgam+h3tgam)*2._dp*zb(p2,p3)
     &  +(-im*h2tgam+h4tgam)*zb(p3,p4)*za(p4,p5)*zb(p2,p5))
      anomZgam(1,2,1,2)=-1._dp/(4._dp*s126*za(p1,p6)*za(p2,p6))
     & *za(p4,p5)*za(p1,p5)
     & *((-im*h1tgam+h3tgam)
     &   *2._dp*(zb(p3,p2)*za(p2,p1)+zb(p3,p6)*za(p6,p1))
     &  +(-im*h2tgam+h4tgam)
     &   *zb(p3,p4)*za(p4,p5)*(zb(p5,p2)*za(p2,p1)+zb(p5,p6)*za(p6,p1)))
     
      anomZgam(1,2,2,2)=conjg(anomZgam(2,1,1,1))
      anomZgam(2,1,2,2)=conjg(anomZgam(1,2,1,1))
      anomZgam(2,2,1,2)=conjg(anomZgam(1,1,2,1))
      anomZgam(2,2,2,1)=conjg(anomZgam(1,1,1,2))

      anomZgam(2,2,2,2)=conjg(anomZgam(1,1,1,1))
      anomZgam(2,2,1,1)=conjg(anomZgam(1,1,2,2))
      anomZgam(2,1,2,1)=conjg(anomZgam(1,2,1,2))
      anomZgam(2,1,1,2)=conjg(anomZgam(1,2,2,1))
            
      return
      end

