      subroutine gggghn_amp(p1,p2,p3,p4,p,n,
     & Hgggg_1256,Hgggg_1265,Hgggg_1625)
      implicit none
      include 'types.f'
C     The function gggghn calculates g(p1)+g(p2)-->H+g(p3)+g(p4)
C     with p1 contracted with n
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      integer::j,p1,p2,p3,p4,h1,h2,h3,h4
      real(dp)::Hgggg_1256,Hgggg_1265,Hgggg_1625,
     & n(4),p(mxpart,4)
      complex(dp)::amp(3,2,2,2,2),
     &  amppp(3),apmpp(3),appmp(3),apppm(3),
     &  apppp(3),
     &  ammpp(3),ampmp(3),amppm(3),apmmp(3),apmpm(3),appmm(3),
     &  ndep,ndem,vecm(mxpart,mxpart)

      call checkndotp(p,n,p1)
      
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do j=1,3
      amp(j,h1,h2,h3,h4)=czip
      enddo
      enddo
      enddo
      enddo
      enddo

      call makepppp(p1,p2,p3,p4,za,apppp)
      call makemppp(p1,p2,p3,p4,za,zb,amppp,apmpp,appmp,apppm)
      call makemmpp(p1,p2,p3,p4,za,zb,
     & ammpp,ampmp,amppm,apmmp,apmpm,appmm)



      do j=1,3
      amp(j,2,2,2,2)=apppp(j)
      amp(j,1,2,2,2)=amppp(j)
      amp(j,2,1,2,2)=apmpp(j)
      amp(j,2,2,1,2)=appmp(j)
      amp(j,2,2,2,1)=apppm(j)

      amp(j,1,1,2,2)=ammpp(j)
      amp(j,1,2,1,2)=ampmp(j)
      amp(j,1,2,2,1)=amppm(j)
      amp(j,2,1,1,2)=apmmp(j)
      amp(j,2,1,2,1)=apmpm(j)
      amp(j,2,2,1,1)=appmm(j)
      enddo
      

      do j=1,3
      amp(j,1,1,1,1)=conjg(amp(j,2,2,2,2))
      amp(j,2,1,1,1)=conjg(amp(j,1,2,2,2))
      amp(j,1,2,1,1)=conjg(amp(j,2,1,2,2))
      amp(j,1,1,2,1)=conjg(amp(j,2,2,1,2))
      amp(j,1,1,1,2)=conjg(amp(j,2,2,2,1))

      amp(j,1,1,2,2)=conjg(amp(j,2,2,1,1))
      amp(j,1,2,1,2)=conjg(amp(j,2,1,2,1))
      amp(j,1,2,2,1)=conjg(amp(j,2,1,1,2))
      amp(j,2,1,1,2)=conjg(amp(j,1,2,2,1))
      amp(j,2,1,2,1)=conjg(amp(j,1,2,1,2))
      amp(j,2,2,1,1)=conjg(amp(j,1,1,2,2))
      enddo

c      call ggggh1amp(p1,p2,p3,p4,p,diag)      

c      do h1=1,2
c      do h2=1,2
c      do h3=1,2
c      do h4=1,2
c      write(6,*) h1,h2,h3,h4,diag(h1,h2,h3,h4)/amp(1,h1,h2,h3,h4)/eight
c      enddo
c      enddo
c      enddo
c      enddo

c      pause
      call ndveccur(p1,p2,n,p,vecm)
      ndep=vecm(p2,p1)/rt2/za(p2,p1)
      ndem=vecm(p1,p2)/rt2/zb(p1,p2)
 
      Hgggg_1256=zip
      Hgggg_1265=zip
      Hgggg_1625=zip

      do h2=1,2
      do h3=1,2
      do h4=1,2
C         write(*,*) 'h4g:  ',h1,h2,h3,h4
C         write(*,*) 'h4g:  ',amp(1,h1,h2,h3,h4),amp(2,h1,h2,h3,h4),
C     &        amp(3,h1,h2,h3,h4)
C         write(*,*) 'h4gsq:',xn**2*V/two*abs(amp(1,h1,h2,h3,h4))**2
C     &        ,xn**2*V/two*abs(amp(2,h1,h2,h3,h4))**2,
C     &        xn**2*V/two*abs(amp(3,h1,h2,h3,h4))**2
      Hgggg_1256=Hgggg_1256
     & +abs(ndep*amp(1,1,h2,h3,h4)+ndem*amp(1,2,h2,h3,h4))**2
      Hgggg_1265=Hgggg_1265
     & +abs(ndep*amp(2,1,h2,h3,h4)+ndem*amp(2,2,h2,h3,h4))**2
      Hgggg_1625=Hgggg_1625
     & +abs(ndep*amp(3,1,h2,h3,h4)+ndem*amp(3,2,h2,h3,h4))**2
      enddo
      enddo
      enddo

c===  (1/4 ---> 1/2) because only three orderings)
      Hgggg_1256=xn**2*V/two*Hgggg_1256
      Hgggg_1265=xn**2*V/two*Hgggg_1265
      Hgggg_1625=xn**2*V/two*Hgggg_1625

c      Hgggg=Hgggg_1256+Hgggg_1265+Hgggg_1625

      return
      end


