      subroutine qqb_zgam_frag(p,msq)
      implicit none
C-----Author John Campbell
C-----June 2000
!---- Modifed to include Fragmentation by 
!---- C Williams Dec 2010
c----Matrix element for Z production
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->(e^-(p3)+e^+(p4))+g(p5)
c---
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'frag.f'
      integer j,k,i
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      double complex AqqbZg(2,2,2),AqbqZg(2,2,2),AqgZq(2,2,2),
     .               AqbgZqb(2,2,2),AgqbZqb(2,2,2),AgqZq(2,2,2),prop
      double precision fsq,D(0:5)
      common/D/D
!$omp threadprivate(/D/)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo
      
      fsq=frag_scale**2
c---- Generate array D(j) corresponding to MCFM notation 0=gluon 1=up 2=down ....
        do i=0,5
           D(i)=0d0
           if     (fragset .eq. 'BFGset_I') then
              call get_frag(z_frag,fsq,1,i,D(i))   
           elseif (fragset .eq. 'BFGsetII') then  
              call get_frag(z_frag,fsq,2,i,D(i))   
           elseif (fragset .eq. 'GdRG__LO') then 
            call GGdR_frag(z_frag,i,D(i),0)
           else
              write(6,*) 'Unrecognized fragmentation set name: ',fragset
              stop        
           endif
        enddo


      call dotem(5,p,s)
      call spinoru(5,p,za,zb)
      
c---protect from soft and collinear singularities
c      if  ((-s(1,5) .lt. cutoff) .or. (-s(2,5) .lt. cutoff)) return

C-----Protect from photon pole by cutting off at some value about 10 GeV
c      if (s(3,4) .lt. 4d0*mbsq) return

      prop=s(3,4)/Dcmplx((s(3,4)-zmass**2),zmass*zwidth)
      fac=4d0*V*esq**2*gsq

c      qqbZg= +aveqq*s(3,4)**2*fac*z1jet(1,2,3,4,5)
c      gqbZqb=-aveqg*s(3,4)**2*fac*z1jet(5,2,3,4,1)
c      qgZq=  -aveqg*s(3,4)**2*fac*z1jet(1,5,3,4,2)      
c      qbqZg= +aveqq*s(3,4)**2*fac*z1jet(2,1,3,4,5)
c      qbgZqb=-aveqg*s(3,4)**2*fac*z1jet(5,1,3,4,2)
c      gqZq=  -aveqg*s(3,4)**2*fac*z1jet(2,5,3,4,1)

      call zgamps(1,2,3,4,5,za,zb,AqqbZg)
      call zgamps(5,2,3,4,1,za,zb,AgqbZqb)
      call zgamps(1,5,3,4,2,za,zb,AqgZq)
      call zgamps(2,1,3,4,5,za,zb,AqbqZg)
      call zgamps(5,1,3,4,2,za,zb,AqbgZqb)
      call zgamps(2,5,3,4,1,za,zb,AgqZq)
      
      do j=-nf,nf
      do k=-nf,nf
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j .eq. 0) .and. (k .eq. 0)) then
          msq(j,k)=0d0
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          msq(j,k)=cdabs((Q(j)*q1+L(j)*l1*prop)*AqqbZg(1,1,1))**2
     .            +cdabs((Q(j)*q1+L(j)*l1*prop)*AqqbZg(1,1,2))**2
     .            +cdabs((Q(j)*q1+L(j)*r1*prop)*AqqbZg(1,2,1))**2
     .            +cdabs((Q(j)*q1+L(j)*r1*prop)*AqqbZg(1,2,2))**2
     .            +cdabs((Q(j)*q1+R(j)*l1*prop)*AqqbZg(2,1,1))**2
     .            +cdabs((Q(j)*q1+R(j)*l1*prop)*AqqbZg(2,1,2))**2
     .            +cdabs((Q(j)*q1+R(j)*r1*prop)*AqqbZg(2,2,1))**2
     .            +cdabs((Q(j)*q1+R(j)*r1*prop)*AqqbZg(2,2,2))**2
          msq(j,k)=msq(j,k)*aveqq*fac/s(3,4)**2*D(0)
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          msq(j,k)=cdabs((Q(k)*q1+L(k)*l1*prop)*AqbqZg(1,1,1))**2
     .            +cdabs((Q(k)*q1+L(k)*l1*prop)*AqbqZg(1,1,2))**2
     .            +cdabs((Q(k)*q1+L(k)*r1*prop)*AqbqZg(1,2,1))**2
     .            +cdabs((Q(k)*q1+L(k)*r1*prop)*AqbqZg(1,2,2))**2
     .            +cdabs((Q(k)*q1+R(k)*l1*prop)*AqbqZg(2,1,1))**2
     .            +cdabs((Q(k)*q1+R(k)*l1*prop)*AqbqZg(2,1,2))**2
     .            +cdabs((Q(k)*q1+R(k)*r1*prop)*AqbqZg(2,2,1))**2
     .            +cdabs((Q(k)*q1+R(k)*r1*prop)*AqbqZg(2,2,2))**2
          msq(j,k)=msq(j,k)*aveqq*fac/s(3,4)**2*D(0)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=cdabs((Q(j)*q1+L(j)*l1*prop)*AqgZq(1,1,1))**2
     .            +cdabs((Q(j)*q1+L(j)*l1*prop)*AqgZq(1,1,2))**2
     .            +cdabs((Q(j)*q1+L(j)*r1*prop)*AqgZq(1,2,1))**2
     .            +cdabs((Q(j)*q1+L(j)*r1*prop)*AqgZq(1,2,2))**2
     .            +cdabs((Q(j)*q1+R(j)*l1*prop)*AqgZq(2,1,1))**2
     .            +cdabs((Q(j)*q1+R(j)*l1*prop)*AqgZq(2,1,2))**2
     .            +cdabs((Q(j)*q1+R(j)*r1*prop)*AqgZq(2,2,1))**2
     .            +cdabs((Q(j)*q1+R(j)*r1*prop)*AqgZq(2,2,2))**2
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2*D(j)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=cdabs((Q(-j)*q1+L(-j)*l1*prop)*AqbgZqb(1,1,1))**2
     .            +cdabs((Q(-j)*q1+L(-j)*l1*prop)*AqbgZqb(1,1,2))**2
     .            +cdabs((Q(-j)*q1+L(-j)*r1*prop)*AqbgZqb(1,2,1))**2
     .            +cdabs((Q(-j)*q1+L(-j)*r1*prop)*AqbgZqb(1,2,2))**2
     .            +cdabs((Q(-j)*q1+R(-j)*l1*prop)*AqbgZqb(2,1,1))**2
     .            +cdabs((Q(-j)*q1+R(-j)*l1*prop)*AqbgZqb(2,1,2))**2
     .            +cdabs((Q(-j)*q1+R(-j)*r1*prop)*AqbgZqb(2,2,1))**2
     .            +cdabs((Q(-j)*q1+R(-j)*r1*prop)*AqbgZqb(2,2,2))**2
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2*D(abs(j))
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=cdabs((Q(k)*q1+L(k)*l1*prop)*AgqZq(1,1,1))**2
     .            +cdabs((Q(k)*q1+L(k)*l1*prop)*AgqZq(1,1,2))**2
     .            +cdabs((Q(k)*q1+L(k)*r1*prop)*AgqZq(1,2,1))**2
     .            +cdabs((Q(k)*q1+L(k)*r1*prop)*AgqZq(1,2,2))**2
     .            +cdabs((Q(k)*q1+R(k)*l1*prop)*AgqZq(2,1,1))**2
     .            +cdabs((Q(k)*q1+R(k)*l1*prop)*AgqZq(2,1,2))**2
     .            +cdabs((Q(k)*q1+R(k)*r1*prop)*AgqZq(2,2,1))**2
     .            +cdabs((Q(k)*q1+R(k)*r1*prop)*AgqZq(2,2,2))**2
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2*D(k)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=cdabs((Q(-k)*q1+L(-k)*l1*prop)*AgqbZqb(1,1,1))**2
     .            +cdabs((Q(-k)*q1+L(-k)*l1*prop)*AgqbZqb(1,1,2))**2
     .            +cdabs((Q(-k)*q1+L(-k)*r1*prop)*AgqbZqb(1,2,1))**2
     .            +cdabs((Q(-k)*q1+L(-k)*r1*prop)*AgqbZqb(1,2,2))**2
     .            +cdabs((Q(-k)*q1+R(-k)*l1*prop)*AgqbZqb(2,1,1))**2
     .            +cdabs((Q(-k)*q1+R(-k)*l1*prop)*AgqbZqb(2,1,2))**2
     .            +cdabs((Q(-k)*q1+R(-k)*r1*prop)*AgqbZqb(2,2,1))**2
     .            +cdabs((Q(-k)*q1+R(-k)*r1*prop)*AgqbZqb(2,2,2))**2
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2*D(abs(k))
      endif

   19 continue
      enddo
      enddo
      return
      end
 

c      double precision function z1jet(j1,j2,j3,j4,j5)
c      implicit none 
c      include 'constants.f'
c      include 'sprods_com.f'
c      integer j1,j2,j3,j4,j5
c      double precision s12,s15,s25

c      s12=s(j1,j2)
c      s15=s(j1,j5)
c      s25=s(j2,j5)
c---calculate the propagator
c      z1jet=(s(j1,j4)**2+s(j2,j3)**2)/(s25*s15*s(j3,j4))
cc      z1jet=
cc     . (s12*(2d0*s(j1,j4)*s(j2,j3)+s(j1,j4)*s(j3,j5)+s(j2,j3)*s(j4,j5))
cc     . +s15*(s(j1,j4)*s(j2,j3)+s(j1,j4)*s(j3,j5)-s(j2,j3)*s(j2,j4))
cc     . +s25*(s(j1,j4)*s(j2,j3)+s(j2,j3)*s(j4,j5)-s(j1,j3)*s(j1,j4)))
cc     . /(s15*s25*s(j3,j4)**2)
c      return
c      end


      subroutine zgamps(j1,j2,j3,j4,j5,za,zb,amps)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      double complex amps(2,2,2)
      integer h1,h2,j1,j2,j3,j4,j5
c-- amplitude helicities are amps(quark,lepton,gluon)
      
      amps(1,1,1)=za(j2,j3)/za(j1,j5)/za(j2,j5)
     .             *(za(j2,j1)*zb(j4,j1)+za(j2,j5)*zb(j4,j5))
           
      amps(1,1,2)=zb(j4,j1)/zb(j1,j5)/zb(j2,j5)
     .             *(za(j2,j3)*zb(j2,j1)+za(j3,j5)*zb(j1,j5))
          
      amps(1,2,1)=za(j2,j4)/za(j1,j5)/za(j2,j5)
     .             *(za(j2,j1)*zb(j3,j1)+za(j2,j5)*zb(j3,j5))
     
      amps(1,2,2)=zb(j3,j1)/zb(j1,j5)/zb(j2,j5)
     .             *(za(j2,j4)*zb(j2,j1)+za(j4,j5)*zb(j1,j5))
     
      do h1=1,2
      do h2=1,2
        amps(2,h1,h2)=-dconjg(amps(1,3-h1,3-h2))
      enddo
      enddo
      
      
      return
      end
      
