
      subroutine transform(p,q,x,ip,jp,kp)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     September, 1999.                                                 *
*     Given p (-p1 + -p2 --> p3 ... px .. p_(npart+2))                 *
*     produce q (-q1 + -q2 --> q3 ... qx .. q_(npart+1))               *
*     by Lorentz transformation with jp denoting the vector            *
*     which is removed (ie all components if q(jp) set to zero)        *
*     ip is the emitter, kp is the specatator                          *
*     Correct branch chosen automatically                              *
*     x is x for ii,if,fi and y for ff                                 *
************************************************************************
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'npart.f'
      real(dp):: p(mxpart,4),q(mxpart,4),x,omx,y,omy,
     & k(4),kt(4),ks(4),kDk,ksDks,kDp(3:mxpart),ksDp(3:mxpart)
      integer:: ip,kp,j,nu,jp,ipart,inp,fip

c      do j=1,npart+2
c      do nu=1,4
c        q(j,nu)=0._dp
c      enddo
c      enddo
      q(:,:)=0._dp

      if ((ip <= 2) .and. (kp <= 2)) then
c---initial-initial
        do nu=1,4
        q(ip,nu)=x*p(ip,nu)
        q(kp,nu)=p(kp,nu)
        k(nu) =-p(ip,nu)-p(kp,nu)-p(jp,nu)
        kt(nu) =-x*p(ip,nu)-p(kp,nu)
        ks(nu)=k(nu)+kt(nu)
        enddo

        kDk=k(4)**2-k(1)**2-k(2)**2-k(3)**2
        ksDks=ks(4)**2-ks(1)**2-ks(2)**2-ks(3)**2

        ipart=3
        do j=3,npart+2
        if (j == jp) then
            go to 19
        else
            kDp(j)=k(4)*p(j,4)-k(1)*p(j,1)-k(2)*p(j,2)-k(3)*p(j,3)
            ksDp(j)=ks(4)*p(j,4)-ks(1)*p(j,1)-ks(2)*p(j,2)-ks(3)*p(j,3)
            do nu=1,4
              q(ipart,nu)=p(j,nu)-two*ksDp(j)*ks(nu)/ksDks
     &        +two*kDp(j)*kt(nu)/kDk
            enddo
        ipart=ipart+1
        endif
 19     continue
        enddo
        return
      elseif (((ip <= 2) .and. (kp > 2)) .or.
     &        ((ip > 2) .and. (kp <= 2))) then
c---initial-final or final-initial
        if (ip <= 2) then
          inp=ip
          fip=kp  
        else
          inp=kp
          fip=ip
        endif
        ipart=1
        omx=one-x
           do j=1,npart+2
                do nu=1,4
                   if (j==inp) then
                   q(ipart,nu)=x*p(inp,nu)
                   elseif (j==jp) then
                   goto 20
                   elseif (j==fip) then
                   q(ipart,nu)=p(jp,nu)+p(fip,nu)+omx*p(inp,nu)
                   else
                   q(ipart,nu)=p(j,nu)
                   endif
                enddo
                ipart=ipart+1
 20             continue
           enddo
        return
      elseif ((ip > 2) .and. (kp > 2)) then
c---final-final
        ipart=1
        y=x
        omy=one-y
           do j=1,npart+2
                do nu=1,4
                   if (j==ip) then
                   q(ipart,nu)=p(jp,nu)+p(ip,nu)-y/omy*p(kp,nu)
                   elseif (j==jp) then
                   goto 21
                   elseif (j==kp) then
                   q(ipart,nu)=p(kp,nu)/omy
                   else
                   q(ipart,nu)=p(j,nu)
                   endif
                enddo
                ipart=ipart+1
 21             continue
           enddo
        return
      endif
      end
