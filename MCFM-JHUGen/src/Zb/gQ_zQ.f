      subroutine gQ_zQ(p,msq)
      implicit none
************************************************************************
*    Authors: R.K. Ellis and John Campbell                             *
*    July, 2003.                                                       *
*    Matrix element for Z + heavy quark (of flavour "flav") production *
*    in order alpha_s^2                                                *
*    averaged over initial colours and spins                           *
*     g(-p1)+Q(-p2)-->Z^+(l(p3)+a(p4))+Q(p5)                           *
************************************************************************
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'heavyflav.f'
      integer j,k,hq,hl
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      double complex prop
      double precision
     . AqgZq2(2,2), AqbgZqb2(2,2),AgqbZqb2(2,2),AgqZq2(2,2)
      integer,parameter::swap(2)=(/2,1/)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call dotem(5,p,s)
      call spinoru(5,p,za,zb)

      prop=s(3,4)/Dcmplx((s(3,4)-zmass**2),zmass*zwidth)
      fac=4d0*V*esq**2*gsq

      call zgamps2(1,5,3,4,2,za,zb,AqgZq2)
      call zgamps2(2,5,3,4,1,za,zb,AgqZq2)

      do hq=1,2
      do hl=1,2
        AqbgZqb2(hq,hl)=AqgZq2(hq,swap(hl))
        AgqbZqb2(hq,hl)=AgqZq2(hq,swap(hl))
      enddo
      enddo

      do j=-flav,flav,flav
      do k=-flav,flav,flav

      if( abs(j+k) .ne. flav) goto 20

      if     ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=cdabs(Q(j)*q1+L(j)*l1*prop)**2*AqgZq2(1,1)
     .            +cdabs(Q(j)*q1+L(j)*r1*prop)**2*AqgZq2(1,2)
     .            +cdabs(Q(j)*q1+R(j)*l1*prop)**2*AqgZq2(2,1)
     .            +cdabs(Q(j)*q1+R(j)*r1*prop)**2*AqgZq2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=cdabs(Q(-j)*q1+L(-j)*l1*prop)**2*AqbgZqb2(1,1)
     .            +cdabs(Q(-j)*q1+L(-j)*r1*prop)**2*AqbgZqb2(1,2)
     .            +cdabs(Q(-j)*q1+R(-j)*l1*prop)**2*AqbgZqb2(2,1)
     .            +cdabs(Q(-j)*q1+R(-j)*r1*prop)**2*AqbgZqb2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=cdabs(Q(k)*q1+L(k)*l1*prop)**2*AgqZq2(1,1)
     .            +cdabs(Q(k)*q1+L(k)*r1*prop)**2*AgqZq2(1,2)
     .            +cdabs(Q(k)*q1+R(k)*l1*prop)**2*AgqZq2(2,1)
     .            +cdabs(Q(k)*q1+R(k)*r1*prop)**2*AgqZq2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=cdabs(Q(-k)*q1+L(-k)*l1*prop)**2*AgqbZqb2(1,1)
     .            +cdabs(Q(-k)*q1+L(-k)*r1*prop)**2*AgqbZqb2(1,2)
     .            +cdabs(Q(-k)*q1+R(-k)*l1*prop)**2*AgqbZqb2(2,1)
     .            +cdabs(Q(-k)*q1+R(-k)*r1*prop)**2*AgqbZqb2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      endif

 20   continue
      enddo
      enddo

      return
      end
