      subroutine transform_mass(p,q,x,ip,jp,kp,misq,mjsq,mksq,mijsq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: R.K. Ellis                                               *
*     June, 2002.                                                      *
*     Given p (-p1 + -p2 --> p3 ... px .. p_(npart+2))                 *
*     produce q (-q1 + -q2 --> q3 ... qx .. q_(npart+1))               *
*     by Lorentz transformation with jp denoting the vector            *
*     which is removed (ie all components if q(jp) set to zero)        *
*     ip is the emitter, kp is the spectator                           *
*     Correct branch chosen automatically                              *
*     x is x for ii,if,fi and y for ff                                 *
************************************************************************
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'npart.f'
      real(dp):: p(mxpart,4),q(mxpart,4),BigQ(4),pij(4),
     & x,omx,k(4),kt(4),ks(4),kDk,ksDks,kDp(3:mxpart),
     & ksDp(3:mxpart),Qsq,rat,misq,mjsq,mksq,mijsq
      real(dp):: pijsq,QDpk,pktilde(4)
      integer:: ip,kp,j,nu,jp,ipart

      do j=1,mxpart
      do nu=1,4
        q(j,nu)=0._dp
      enddo
      enddo

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
            go to 18
        else
            kDp(j)=k(4)*p(j,4)-k(1)*p(j,1)-k(2)*p(j,2)-k(3)*p(j,3)
            ksDp(j)=ks(4)*p(j,4)-ks(1)*p(j,1)-ks(2)*p(j,2)-ks(3)*p(j,3)
            do nu=1,4
              q(ipart,nu)=p(j,nu)-two*ksDp(j)*ks(nu)/ksDks
     &        +two*kDp(j)*kt(nu)/kDk
            enddo
        ipart=ipart+1
        endif
 18     continue
        enddo
        return
        elseif ((ip <= 2) .and. (kp > 2)) then
c---initial-final 
        ipart=1
        omx=one-x
           do j=1,npart+2
                do nu=1,4
                   if (j==ip) then
                   q(ipart,nu)=x*p(ip,nu)
                   elseif (j==jp) then
                   goto 19
                   elseif (j==kp) then
                   q(ipart,nu)=p(jp,nu)+p(kp,nu)+omx*p(ip,nu)
                   else
                   q(ipart,nu)=p(j,nu)
                   endif
                enddo
                ipart=ipart+1
 19             continue
           enddo
        return

        elseif ((ip > 2) .and. (kp <= 2)) then
c---final-initial
        ipart=1
        omx=one-x
           do j=1,npart+2
                do nu=1,4
                   if (j==kp) then
                   q(ipart,nu)=x*p(kp,nu)
                   elseif (j==jp) then
                   goto 20
                   elseif (j==ip) then
                   q(ipart,nu)=p(ip,nu)+p(jp,nu)+omx*p(kp,nu)
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
      do nu=1,4
      BigQ(nu)=p(ip,nu)+p(jp,nu)+p(kp,nu)
      pij(nu)=p(ip,nu)+p(jp,nu)
      enddo   
      Qsq=BigQ(4)**2-BigQ(1)**2-BigQ(2)**2-BigQ(3)**2
      QDpk=
     & +BigQ(4)*p(kp,4)-BigQ(1)*p(kp,1)-BigQ(2)*p(kp,2)-BigQ(3)*p(kp,3)
      pijsq=pij(4)**2-pij(1)**2-pij(2)**2-pij(3)**2
      rat=sqrt((Qsq-mijsq-mksq)**2-4._dp*mijsq*mksq)
      rat=rat/sqrt((Qsq-pijsq-mksq)**2-4._dp*pijsq*mksq)
        ipart=1
           do j=1,npart+2
                do nu=1,4
                   pktilde(nu)=rat*(p(kp,nu)-QDpk/Qsq*BigQ(nu))
     &                +(Qsq+mksq-mijsq)/(2._dp*Qsq)*BigQ(nu)
                   if (j==ip) then
                   q(ipart,nu)=BigQ(nu)-Pktilde(nu)
                   elseif (j==jp) then
                   goto 21
                   elseif (j==kp) then
                   q(ipart,nu)=pktilde(nu)
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
