      subroutine decide_flavour(pflav,pbarflav)
c     ------------------------------------------------------------------
c --- Decide which flavour combination to keep for this event based
c --- on the array of relative flavour constribution weights that 
c --- was formed in lowint.
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      implicit none
      include 'constants.f'
      include 'flavours.f'

      integer pflav,pbarflav

c --- To use VEGAS random number sequence :
c      integer idum
c      COMMON/ranno/idum
      
      double precision ran2
      
      integer j,k
      double precision weight_sum,weight_int
      double precision pointer
      
c     ------------------------------------------------------------------

c --- First add up the weights :
      weight_sum = 0.0d0
      do j=-nf,nf
        do k=-nf,nf
          weight_sum = weight_sum + ppbar_flavours(j,k)
        enddo
      enddo

c --- Now find a random number between zero and this integral :
      pointer = ran2()*weight_sum

c --- Find where this falls in the integral distribution to 
c --- discover the combination :
      weight_int = 0.0d0
      do j=-nf,nf
        do k=-nf,nf
          weight_int = weight_int + ppbar_flavours(j,k)
          if (weight_int .ge. pointer) then
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
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      implicit none
      include 'constants.f'
      include 'eventbuffer.f'

      double precision p(mxpart,4)
      double precision wt 
      integer pflav,pbarflav

      integer i,j

c     ------------------------------------------------------------------

      numstored=numstored+1
      
      if (numstored .gt. buffersize) then
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
c     ------------------------------------------------------------------
c        1         2         3         4         5         6         7
      implicit none
      include 'constants.f'
      include 'maxwt.f'
      include 'eventbuffer.f'

      
      double precision p(mxpart,4)
      double precision wt 
      integer pflav,pbarflav

c --- This common block ensures we are calling VEGAS with the same
c --- parameters as during the weight scan :
      integer itmx1,ncall1,itmx2,ncall2
      double precision integ,integ_err
      common/iterat/itmx1,ncall1,itmx2,ncall2

c --- To use VEGAS random number sequence :
c      integer idum
c      COMMON/ranno/idum

      integer i,j,position

      double precision ran2

      real*4 rannum,randomList(buffersize)

c     ------------------------------------------------------------------

c --- Are there any events left ?
      if (numused.eq.numstored .or. numstored.eq.0) then
c ---   Call integration loop again, this time unweighting :
        write(6,*) 'Trying to unweight events ...'
        unweight = .true.
        numstored=0
        numused=0
 10     call mcfm_vegas(1,itmx2,ncall2,.true.,integ,integ_err)
        write(6,*) 'After event generation, numstored = ',numstored
        if (numstored.eq.0) goto 10
        do i=1,numstored
           rannum = sngl(ran2())
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
