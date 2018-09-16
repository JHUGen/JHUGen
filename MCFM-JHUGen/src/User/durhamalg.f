      subroutine durhamalg(p,npartons,y32,y43,z3,z4,z5,z6)
      implicit none
      include 'types.f'
      
c--- Given a set of 4-vectors, determine the values of the parameters
c--- y32 and y43 at which the event is described in terms of a 2-jet
c--- and 3-jet event, respectively.
c--- Also returns the energy fractions of the jets (z3,z4,z5,z6)
c--- If there are at most 3 jets, returns <0 for y43 and for z6
c--- If there are at most 2 jets, returns <0 for y32 and for z5
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),q(mxpart,4),r(mxpart,4),
     & y32,y43,z3,z4,z5,z6,yij,M
      integer:: npartons,i,j,imin,jmin,icount
      
      if     (npartons > 4) then
        write(6,*) 'Too many partons for this version of the Durham'
        write(6,*) 'algorithm; at most 4 are allowed.'
        stop
      elseif (npartons < 2) then
        write(6,*) 'Too few partons for this version of the Durham'
        write(6,*) 'algorithm; need at least 2.'
        stop
      endif
      
      y43=-1._dp
      y32=-1._dp
      z6=-1._dp
      z5=-1._dp
      
c--- skip to 2-parton section if there are only 2 partons
      if (npartons == 2) then
        do i=3,4
        do j=1,4
        r(i,j)=p(i,j)
        enddo
        enddo
        goto 33
      endif
      
c--- skip to 3-parton section if there are only 3 partons
      if (npartons == 3) then
        do i=3,5
        do j=1,4
        q(i,j)=p(i,j)
        enddo
        enddo
        goto 22
      endif
      
c--- find lowest yij amongst 4 jets
      y43=0._dp
      
      do i=3,5
      do j=i+1,6
        yij=(p(i,1)*p(j,1)+p(i,2)*p(j,2)+p(i,3)*p(j,3))
     &      /sqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)
     &      /sqrt(p(j,1)**2+p(j,2)**2+p(j,3)**2)
        yij=2._dp*min(p(i,4),p(j,4))**2*(1._dp-yij)
        if ((yij < y43) .or. ((i==3) .and. (j==4))) then
          imin=i
          jmin=j
          y43=yij
        endif
      enddo
      enddo
      
      z6=min(p(3,4),p(4,4),p(5,4),p(6,4))

c--- make momenta for 3-jet event by adding 4-vectors
c---  (NO LONGER use P scheme, i.e. add 3 vectors
c---  and rescale energy to make massless)
      icount=3
      do i=3,6
      do j=1,4
        if     (i == imin) then
        q(icount,j)=p(imin,j)+p(jmin,j)
        elseif (i == jmin) then
        continue
        else
          q(icount,j)=p(i,j)
        endif
      enddo
      if (i .ne. jmin) icount=icount+1
      enddo
      
c      do i=3,5
c        q(i,4)=sqrt(q(i,1)**2+q(i,2)**2+q(i,3)**2)
c      enddo
      
c--- cluster into 3 jets      
   22 continue

c--- find lowest yij amongst 3 jets
      y32=0._dp
      
      do i=3,4
      do j=i+1,5
        yij=(q(i,1)*q(j,1)+q(i,2)*q(j,2)+q(i,3)*q(j,3))
     &      /sqrt(q(i,1)**2+q(i,2)**2+q(i,3)**2)
     &      /sqrt(q(j,1)**2+q(j,2)**2+q(j,3)**2)
        yij=2._dp*min(q(i,4),q(j,4))**2*(1._dp-yij)
        if ((yij < y32) .or. ((i==3) .and. (j==4))) then
          imin=i
          jmin=j
          y32=yij
        endif
      enddo
      enddo
      
      z5=min(q(3,4),q(4,4),q(5,4))

c--- make momenta for 3-jet event by adding 4-vectors
c---  (NO LONGER use P scheme, i.e. add 3 vectors
c---  and rescale energy to make massless)
      icount=3
      do i=3,5
      do j=1,4
        if     (i == imin) then
        r(icount,j)=q(imin,j)+q(jmin,j)
      elseif (i == jmin) then
        continue
      else
        r(icount,j)=q(i,j)
      endif
      enddo
      if (i .ne. jmin) icount=icount+1
      enddo
      
c      do i=3,4
c        r(i,4)=sqrt(r(i,1)**2+r(i,2)**2+r(i,3)**2)
c      enddo
      
c--- now we have 2 jets      
   33 continue

      z4=min(r(3,4),r(4,4))
      z3=max(r(3,4),r(4,4))

c--- Now rescale all properties by onium mass (M) appropriately
      M=-p(1,4)-p(2,4)
      
      y43=y43/M**2
      y32=y32/M**2
      
      z3=2._dp*z3/M
      z4=2._dp*z4/M
      z5=2._dp*z5/M
      z6=2._dp*z6/M
      
      return
      end
      
