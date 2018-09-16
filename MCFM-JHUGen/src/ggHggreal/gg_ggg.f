      function GG_GGG(IHEL)
C
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR PROCESS : g g -> g g g (h)
C
      IMPLICIT NONE
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      real(dp) GG_GGG
      integer, parameter::ngluons=5

      integer::IHEL(NGLUONS)
      complex(dp)::amp_h5g

      integer::PERM(NGLUONS),IHELX(NGLUONS)
      complex(dp)::Z(6),ZTEMP
      integer::I,J

c      logical::CHECKS
c      DATA CHECKS/.false./

c--   permutations
      integer,parameter::IP(ngluons,6)=reshape((/
     & 1,2,3,4,5,
     & 1,2,4,3,5,
     & 1,3,2,4,5,
     & 1,3,4,2,5,
     & 1,4,2,3,5,
     % 1,4,3,2,5/),(/ngluons,6/))

c--   color matrix
c     c_ij= (N^2-1)*N^3/4*CIJ/2 in the GM normalization

      integer,parameter::CIJ(6,6)=reshape((/
     & 8,4,4,2,2,0,
     & 4,8,2,0,4,2,
     & 4,2,8,4,0,2,
     & 2,0,4,8,2,4,
     & 2,4,0,2,8,4,
     & 0,2,2,4,4,8/),(/6,6/))
      include 'cplx.h'


C-----
C  BEGIN CODE
C-----

      GG_GGG = zip

         DO I=1,6

         PERM(1)=IP(1,I)
         PERM(2)=IP(2,I)
         PERM(3)=IP(3,I)
         PERM(4)=IP(4,I)
         PERM(5)=IP(5,I)


         call iperm(IHEL,PERM,IHELX,NGLUONS)


         Z(I)=amp_h5g(PERM,IHELX)

c
c    check on the properties of the amplitudes
c
c         if (checks) then
c         call checker(z,i,perm,ihelx)
c         endif
c=============================================================

         ENDDO !sum over permutations

c
c   sum over color
c
         DO J = 1, 6
            ZTEMP = czip
            DO I = 1, 6
               ZTEMP = ZTEMP + Z(I)*CIJ(I,J)
            ENDDO
            GG_GGG =GG_GGG+real(ZTEMP*conjg(Z(J)))
         ENDDO



c     Overall normalization of CIJ
         GG_GGG=GG_GGG*XN**3*(XN**2-one)/four

      RETURN
      END


