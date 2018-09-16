      subroutine qqb_wbb_v(P,msqv)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     July, 1998.                                                      *
*     Calculate the virtual matrix element squared and subtraction     *
*     terms for the process                                            *
*     q(-p1) +Q(-p6)+ l(-p4) -->   q(p2)+Q(p5) +l(p3)                  *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scheme.f'
      include 'masses.f'
      include 'noglue.f'
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),
     & p(mxpart,4),q(mxpart,4),faclo,fac,
     & qqb,qbq
      complex(dp):: atrLLL,atrLRL,a61LLL,a61LRL
      complex(dp):: tLLL,tLRL,fLLL,fLRL
      integer:: nu,j,k

      scheme='dred'

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      enddo
      enddo

c--- shortcut if we're doing gqonly
      if (gqonly) return

c---calculate the lowest order matrix element and fill the common block
c---twopij with s_{ij}
      call qqb_wbb(p,msq)
      if (
     &      (s(5,6) < four*mbsq)
     & .or. (s(1,5)*s(2,5)/s(1,2) < mbsq)
     & .or. (s(1,6)*s(2,6)/s(1,2) < mbsq) ) return

c---  Now transform momenta into a notation
c---  suitable for calling the BDKW function with notation which is
c---  q-(-p4)+Q+(-p2)+l-(-p5) ---> q+(p1)+Q-(p3)+l+(p6)
      do nu=1,4
      q(1,nu)=p(2,nu)
      q(2,nu)=p(6,nu)
      q(3,nu)=p(5,nu)
      q(4,nu)=p(1,nu)
      q(5,nu)=p(4,nu)
      q(6,nu)=p(3,nu)
      enddo

      call spinoru(6,q,za,zb)
      faclo=V*aveqq*gw**4*gsq**2
      fac=faclo*xn*0.5_dp*ason2pi

c----do whatever needs to be done q-qb case
      tLLL=atrLLL(1,2,3,4,5,6,za,zb)
      fLLL=a61LLL(1,2,3,4,5,6,za,zb)
      tLRL=atrLRL(1,2,3,4,5,6,za,zb)
      fLRL=a61LRL(1,2,3,4,5,6,za,zb)

      qqb=fac*real(tLLL*conjg(fLLL)+fLLL*conjg(tLLL)
     &            +tLRL*conjg(fLRL)+fLRL*conjg(tLRL))

c----now look at qb-q case, swap the momenta
      tLLL=atrLLL(4,2,3,1,5,6,za,zb)
      fLLL=a61LLL(4,2,3,1,5,6,za,zb)
      tLRL=atrLRL(4,2,3,1,5,6,za,zb)
      fLRL=a61LRL(4,2,3,1,5,6,za,zb)

      qbq=fac*real(tLLL*conjg(fLLL)+fLLL*conjg(tLLL)
     &            +tLRL*conjg(fLRL)+fLRL*conjg(tLRL))

      do j=-(nf-1),(nf-1)
      do k=-(nf-1),(nf-1)
      if (Vsq(j,k) == 0._dp) goto 20
            if     ((j > 0) .and. (k < 0)) then
               msqv(j,k)=Vsq(j,k)*qqb
            elseif ((j < 0) .and. (k > 0)) then
               msqv(j,k)=Vsq(j,k)*qbq
            else
               msqv(j,k)=0._dp
            endif
 20   continue
      enddo
      enddo

      return
      end

