      subroutine decide_flavour(pflav,pbarflav)
      implicit none
      include 'types.f'
c     ------------------------------------------------------------------
c --- Decide which flavour combination to keep for this event based
c --- on the array of relative flavour constribution weights that 
c --- was formed in lowint.
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'flavours.f'

      integer:: pflav,pbarflav

c --- To use VEGAS random number sequence :
c      integer:: idum
c      COMMON/ranno/idum
      
      real(dp):: ran2nr
      
      integer:: j,k
      real(dp):: weight_sum,weight_int
      real(dp):: pointer
      
c     ------------------------------------------------------------------

c --- First add up the weights :
      weight_sum = 0.0_dp
      do j=-nf,nf
        do k=-nf,nf
          weight_sum = weight_sum + ppbar_flavours(j,k)
        enddo
      enddo

c --- Now find a random number between zero and this integral :
      pointer = ran2nr()*weight_sum

c --- Find where this falls in the integral distribution to 
c --- discover the combination :
      weight_int = 0.0_dp
      do j=-nf,nf
        do k=-nf,nf
          weight_int = weight_int + ppbar_flavours(j,k)
          if (weight_int >= pointer) then
            pflav = j
            pbarflav = k
            goto 10
          endif
        enddo
      enddo

 10   return
      end
c

      subroutine storeevent(p,wt,pflav,pbarflav)
      implicit none
      include 'types.f'
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'eventbuffer.f'

      real(dp):: p(mxpart,4)
      real(dp):: wt 
      integer:: pflav,pbarflav

      integer:: i,j

c     ------------------------------------------------------------------

      numstored=numstored+1
      
      if (numstored > buffersize) then
        write(6,*) 'ERROR : storeevent : internal buffer size exceeded'
        numstored=buffersize
        return
      else
        do i=1,mxpart
          do j=1,4
            eventbuffer(numstored,i,j)=p(i,j)
          enddo
        enddo
        wtbuffer(numstored)=wt
        pflavbuffer(numstored)=pflav
        pbarflavbuffer(numstored)=pbarflav
      endif
  
      return
      end

c

      subroutine mcfm_getevent(p,wt,pflav,pbarflav)
      implicit none
      include 'types.f'
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'maxwt.f'
      include 'eventbuffer.f'
      include 'iterat.f'

      
      real(dp):: p(mxpart,4)
      real(dp):: wt 
      integer:: pflav,pbarflav

c --- This common block ensures we are calling VEGAS with the same
c --- parameters as during the weight scan :
      real(dp):: integ,integ_err

c --- To use VEGAS random number sequence :
c      integer:: idum
c      COMMON/ranno/idum

      integer:: i,j,position

      real(dp):: ran2nr

      real*4 rannum,randomList(buffersize)

c     ------------------------------------------------------------------

c --- Are there any events left ?
      if (numused==numstored .or. numstored==0) then
c ---   Call integration loop again, this time unweighting :
        write(6,*) 'Trying to unweight events ...'
        unweight = .true.
        numstored=0
        numused=0
 10     call mcfm_vegas(1,itmx2,ncall2,.true.,integ,integ_err)
        write(6,*) 'After event generation, numstored = ',numstored
        if (numstored==0) goto 10
        do i=1,numstored
           rannum = sngl(real(ran2nr(),kind=8))
           randomList(i) = rannum
           indexList(i) = 0
        enddo
        call SORTZV(randomList,indexList,numstored,1,0,0)
      endif
      
c --- Now pick an event from the list :
      numused = numused+1
      position = indexList(numused)

      do i=1,mxpart
        do j=1,4
          p(i,j)=eventbuffer(position,i,j)
        enddo
      enddo
      
      wt=wtbuffer(position)
      pflav=pflavbuffer(position)
      pbarflav=pbarflavbuffer(position)

      return
      end

c


c# 10 "sortzv.F" 2
      SUBROUTINE SORTZV (A,INDEX,N1,MODE,NWAY,NSORT)
C
C CERN PROGLIB# M101    SORTZV          .VERSION KERNFOR  3.15  820113
C ORIG. 02/10/75
C
      DIMENSION A(N1),INDEX(N1)
C
C
      N = N1
      IF (N.LE.0)            RETURN
      IF (NSORT.NE.0) GO TO 2
      DO 1 I=1,N
    1 INDEX(I)=I
C
    2 IF (N.EQ.1)            RETURN
      IF (MODE)    10,20,30
   10 CALL SORTTI (A,INDEX,N)
      GO TO 40
C
   20 CALL SORTTC(A,INDEX,N)
      GO TO 40
C
   30 CALL SORTTF (A,INDEX,N)
C
   40 IF (NWAY.EQ.0) GO TO 50
      N2 = N/2
      DO 41 I=1,N2
      ISWAP = INDEX(I)
      K = N+1-I
      INDEX(I) = INDEX(K)
   41 INDEX(K) = ISWAP
   50 RETURN
      END
*     ========================================
      SUBROUTINE SORTTF (A,INDEX,N1)
C
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF (AI.LE.A (I22)) GO TO 3
      INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (A(I22)-A(I222)) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (AI-A(I22)) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
*     ========================================
      SUBROUTINE SORTTI (A,INDEX,N1)
C
      INTEGER A,AI
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF (AI.LE.A (I22)) GO TO 3
      INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (A(I22)-A(I222)) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (AI-A(I22)) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
*     ========================================
      SUBROUTINE SORTTC (A,INDEX,N1)
C
      INTEGER A,AI
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF(ICMPCH(AI,A(I22)))3,3,21
   21 INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (ICMPCH(A(I22),A(I222))) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (ICMPCH(AI,A(I22))) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
*     ========================================
      FUNCTION ICMPCH(IC1,IC2)
C     FUNCTION TO COMPARE TWO 4 CHARACTER EBCDIC STRINGS - IC1,IC2
C     ICMPCH=-1 IF HEX VALUE OF IC1 IS LESS THAN IC2
C     ICMPCH=0  IF HEX VALUES OF IC1 AND IC2 ARE THE SAME
C     ICMPCH=+1 IF HEX VALUES OF IC1 IS GREATER THAN IC2
      I1=IC1
      I2=IC2
      IF(I1.GE.0.AND.I2.GE.0)GOTO 40
      IF(I1.GE.0)GOTO 60
      IF(I2.GE.0)GOTO 80
      I1=-I1
      I2=-I2
      IF(I1-I2)80,70,60
 40   IF(I1-I2)60,70,80
 60   ICMPCH=-1
      RETURN
 70   ICMPCH=0
      RETURN
 80   ICMPCH=1
      RETURN
      END


