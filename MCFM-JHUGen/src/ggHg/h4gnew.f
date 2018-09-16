      subroutine h4gnew(p1,p2,p3,p4,
     &                  Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j,p1,p2,p3,p4,h1,h2,h3,h4
      real(dp):: Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625
      complex(dp):: amp(3,2,2,2,2)

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

      call Amplo(p1,p2,p3,p4,amp)
      Hgggg_1256=zip
      Hgggg_1265=zip
      Hgggg_1625=zip

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      Hgggg_1256=Hgggg_1256+abs(amp(1,h1,h2,h3,h4))**2
      Hgggg_1265=Hgggg_1265+abs(amp(2,h1,h2,h3,h4))**2
      Hgggg_1625=Hgggg_1625+abs(amp(3,h1,h2,h3,h4))**2
      enddo
      enddo
      enddo
      enddo

C===  (1/4 ---> 1/2) because only three orderings)
      Hgggg_1256=xn**2*V/2._dp*Hgggg_1256
      Hgggg_1265=xn**2*V/2._dp*Hgggg_1265
      Hgggg_1625=xn**2*V/2._dp*Hgggg_1625

      Hgggg=Hgggg_1256+Hgggg_1265+Hgggg_1625

      return
      end


c--- this routine is a wrapper to the new versions of the leading
c--- order amplitudes that are calculated in the same way as the
c--- new virtual amplitudes
      subroutine Amplo(j1,j2,j3,j4,A)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      complex(dp):: A(3,2,2,2,2)
      complex(dp):: A0Hggggpppp,A0Hggggpmmm,A0Hggggmmpp,A0Hggggmpmp
      integer:: j1,j2,j3,j4,i1(3),i2(3),i3(3),i4(4),j

      i1(1)=j1
      i2(1)=j2
      i3(1)=j3
      i4(1)=j4

      i1(2)=j1
      i2(2)=j2
      i3(2)=j4
      i4(2)=j3

      i1(3)=j1
      i2(3)=j4
      i3(3)=j2
      i4(3)=j3

c--- NB: faster to do complex conjugation rather than permute,
c---     valid for these Born amplitudes only
      do j=1,3
      A(j,2,2,2,2)=A0Hggggpppp(i1(j),i2(j),i3(j),i4(j),za,zb)

      A(j,2,1,1,1)=A0Hggggpmmm(i1(j),i2(j),i3(j),i4(j),za,zb)
      A(j,1,2,1,1)=A0Hggggpmmm(i2(j),i3(j),i4(j),i1(j),za,zb)
      A(j,1,1,2,1)=A0Hggggpmmm(i3(j),i4(j),i1(j),i2(j),za,zb)
      A(j,1,1,1,2)=A0Hggggpmmm(i4(j),i1(j),i2(j),i3(j),za,zb)

      A(j,1,1,2,2)=A0Hggggmmpp(i1(j),i2(j),i3(j),i4(j),za,zb)
      A(j,1,2,2,1)=A0Hggggmmpp(i4(j),i1(j),i2(j),i3(j),za,zb)

      A(j,1,2,1,2)=A0Hggggmpmp(i1(j),i2(j),i3(j),i4(j),za,zb)

c      A(j,1,1,1,1)=A0Hggggpppp(i1(j),i2(j),i3(j),i4(j),zb,za)
c      A(j,1,2,2,2)=A0Hggggpmmm(i1(j),i2(j),i3(j),i4(j),zb,za)
c      A(j,2,1,2,2)=A0Hggggpmmm(i2(j),i3(j),i4(j),i1(j),zb,za)
c      A(j,2,2,1,2)=A0Hggggpmmm(i3(j),i4(j),i1(j),i2(j),zb,za)
c      A(j,2,2,2,1)=A0Hggggpmmm(i4(j),i1(j),i2(j),i3(j),zb,za)
c      A(j,2,2,1,1)=A0Hggggmmpp(i3(j),i4(j),i1(j),i2(j),za,zb)
c      A(j,2,1,1,2)=A0Hggggmmpp(i2(j),i3(j),i4(j),i1(j),za,zb)
c      A(j,2,1,2,1)=A0Hggggmpmp(i2(j),i3(j),i4(j),i1(j),za,zb)
      A(j,1,1,1,1)=conjg(A(j,2,2,2,2))

      A(j,1,2,2,2)=conjg(A(j,2,1,1,1))
      A(j,2,1,2,2)=conjg(A(j,1,2,1,1))
      A(j,2,2,1,2)=conjg(A(j,1,1,2,1))
      A(j,2,2,2,1)=conjg(A(j,1,1,1,2))

      A(j,2,2,1,1)=conjg(A(j,1,1,2,2))
      A(j,2,1,1,2)=conjg(A(j,1,2,2,1))

      A(j,2,1,2,1)=conjg(A(j,1,2,1,2))
      enddo

      return
      end


c--- this routine is a wrapper to the old leading order amplitudes
c--- that are based on the calculations of Kauffman, Desai and Risal
      subroutine Amplo1(p1,p2,p3,p4,amp)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      integer:: j,p1,p2,p3,p4
      complex(dp):: amp(3,2,2,2,2),
     &  amppp(3),apmpp(3),appmp(3),apppm(3),
     &  apppp(3),
     &  ammpp(3),ampmp(3),amppm(3),apmmp(3),apmpm(3),appmm(3)

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
      return
      end

