      subroutine msq_qqQQg(i1,i2,i3,i4,i5,i6,i7,xsq)
      implicit none
      include 'types.f'
c--- this will compute Z -> q qb Q Qb g amplitudes according to NT
c--- note that we simplify the situation by assuming that q != Q
c--- in which case all (1<-->3) and (2<-->4) amps vanish
c--- notation:  quark pairs (q,qb) are (1,2) and (3,4)
c---            lepton-antilepton (6,7) and gluon 5
c--- note that this is now updated so that the routine returns
c--- the matrix element separated according to coupling
c--- returned is msq(hq,Qh,hg,lh,j) where
c---   j=1 : Z couples to the (1,2) line in M and M*
c---   j=2 : Z couples to the (3,4) line in M and M*
c---   j=3 : Z couples to the (1,2) line in M and the (3,4) in M* (RE)
c---   j=4 : Z couples to the (1,2) line in M and the (3,4) in M* (IM)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: i1,i2,i3,i4,i5,i6,i7,j,hq,Qh,hg,lh
      complex(dp)::a1(2,2,2,2),a2(2,2,2,2),a3(2,2,2,2),a4(2,2,2,2),
     &               b1(2,2,2,2),b2(2,2,2,2),b3(2,2,2,2),b4(2,2,2,2)
      complex(dp)::temp
      real(dp)::xA(4),xD(4),xF(4),xG(4)
      real(dp)::xy(4),xz(4),xxy(4),xsq(2,2,2,2,4)

c--- note that 'a' corresponds to the perm 1,2,3,4,5
      call nagyqqQQg(i1,i2,i3,i4,i5,i6,i7,a1,a2,a3,a4)
c---           'b' corresponds to the perm 3,4,1,2,5
      call nagyqqQQg(i3,i4,i1,i2,i5,i6,i7,b1,b2,b3,b4)

      do hq=1,2
      do Qh=1,2
      do hg=1,2
      do lh=1,2

c--- eqn (B57)
c        mbar1=a1(hq,Qh,hg,lh)+b3(Qh,hq,hg,lh)
c        mbar2=a2(hq,Qh,hg,lh)+b4(Qh,hq,hg,lh)
c--- introduce new notation for eqn (B57) with (1<-->3, 2<-->4)
c--- mbar3 = m1(3,4,1,2) + m3(1,2,3,4)
c--- mbar4 = m2(3,4,1,2) + m4(1,2,3,4)
c        mbar3=b1(Qh,hq,hg,lh)+a3(hq,Qh,hg,lh)
c        mbar4=b2(Qh,hq,hg,lh)+a4(hq,Qh,hg,lh)
c--- eqn (B50)
c        mA=abs(mbar1)**2+abs(mbar2)**2
c     &    +abs(mbar3)**2+abs(mbar4)**2
        xA(1)=abs(a1(hq,Qh,hg,lh))**2+abs(a2(hq,Qh,hg,lh))**2
     &       +abs(a3(hq,Qh,hg,lh))**2+abs(a4(hq,Qh,hg,lh))**2
        xA(2)=abs(b1(Qh,hq,hg,lh))**2+abs(b2(Qh,hq,hg,lh))**2
     &       +abs(b3(Qh,hq,hg,lh))**2+abs(b4(Qh,hq,hg,lh))**2
        temp =2._dp*(a1(hq,Qh,hg,lh)*conjg(b3(Qh,hq,hg,lh))
     &            +a2(hq,Qh,hg,lh)*conjg(b4(Qh,hq,hg,lh))
     &            +a3(hq,Qh,hg,lh)*conjg(b1(Qh,hq,hg,lh))
     &            +a4(hq,Qh,hg,lh)*conjg(b2(Qh,hq,hg,lh)))
        xA(3)=real(temp)
        xA(4)=aimag(temp)
c--- eqn (B51) - vanishes because different quarks
c        mB=0._dp
c--- eqn (B52) - vanishes because different quarks
c        mC=0._dp
c--- eqn (B53)
c        mD=2._dp*real(mbar1*conjg(mbar2)+mbar3*conjg(mbar4))
        xD(1)=2._dp*real(a1(hq,Qh,hg,lh)*conjg(a2(hq,Qh,hg,lh))
     &                 +a3(hq,Qh,hg,lh)*conjg(a4(hq,Qh,hg,lh)))
        xD(2)=2._dp*real(b1(Qh,hq,hg,lh)*conjg(b2(Qh,hq,hg,lh))
     &                 +b3(Qh,hq,hg,lh)*conjg(b4(Qh,hq,hg,lh)))
        temp =2._dp*(a1(hq,Qh,hg,lh)*conjg(b4(Qh,hq,hg,lh))
     &            +conjg(b3(Qh,hq,hg,lh))*a2(hq,Qh,hg,lh)
     &            +a3(hq,Qh,hg,lh)*conjg(b2(Qh,hq,hg,lh))
     &            +conjg(b1(Qh,hq,hg,lh))*a4(hq,Qh,hg,lh))
        xD(3)=real(temp)
        xD(4)=aimag(temp)
c--- eqn (B54) - vanishes because different quarks
c        mE=0._dp
c--- eqn (B55)
c        mF=2._dp*real(mbar1*conjg(mbar3)+mbar2*conjg(mbar4))
        xF(1)=2._dp*real(a1(hq,Qh,hg,lh)*conjg(a3(hq,Qh,hg,lh))
     &                 +a2(hq,Qh,hg,lh)*conjg(a4(hq,Qh,hg,lh)))
        xF(2)=2._dp*real(b1(Qh,hq,hg,lh)*conjg(b3(Qh,hq,hg,lh))
     &                 +b2(Qh,hq,hg,lh)*conjg(b4(Qh,hq,hg,lh)))
        temp =2._dp*(a1(hq,Qh,hg,lh)*conjg(b1(Qh,hq,hg,lh))
     &            +conjg(b3(Qh,hq,hg,lh))*a3(hq,Qh,hg,lh)
     &            +a2(hq,Qh,hg,lh)*conjg(b2(Qh,hq,hg,lh))
     &            +conjg(b4(Qh,hq,hg,lh))*a4(hq,Qh,hg,lh))
        xF(3)=real(temp)
        xF(4)=aimag(temp)
c--- eqn (B56)
c        mG=2._dp*real(mbar1*conjg(mbar4)+mbar2*conjg(mbar3))
        xG(1)=2._dp*real(a1(hq,Qh,hg,lh)*conjg(a4(hq,Qh,hg,lh))
     &                 +a2(hq,Qh,hg,lh)*conjg(a3(hq,Qh,hg,lh)))
        xG(2)=2._dp*real(b1(Qh,hq,hg,lh)*conjg(b4(Qh,hq,hg,lh))
     &                 +b2(Qh,hq,hg,lh)*conjg(b3(Qh,hq,hg,lh)))
        temp=2._dp*(a1(hq,Qh,hg,lh)*conjg(b2(Qh,hq,hg,lh))
     &           +conjg(b3(Qh,hq,hg,lh))*a4(hq,Qh,hg,lh)
     &           +a2(hq,Qh,hg,lh)*conjg(b1(Qh,hq,hg,lh))
     &           +conjg(b4(Qh,hq,hg,lh))*a3(hq,Qh,hg,lh))
        xG(3)=real(temp)
        xG(4)=aimag(temp)
c--- eqns (B44)-(B49)
c        m0=mB+mC+mE
c        mx=-0.5_dp*(3._dp*mC+2._dp*mE+mB)
c        my=mA+mD
c        mz=mF+mG
c        mxx=0.25_dp*(2._dp*mC+mE)
c        mxy=-0.5_dp*(mF+mD)
c--- eqn (B42) translated
c--- note that C3 = (N*CF**2-Tr*CF)/2 = CF/4*(N**2-2)
c--- in my notation, where Tr=1/2
c        msq(hq,Qh,hg,lh)=cf/2._dp*(xn*cf*my+xn**2*mxy+(xn**2-2._dp)/2._dp*mz)
        do j=1,4
        xy(j)=xA(j)+xD(j)
        xz(j)=xF(j)+xG(j)
        xxy(j)=-0.5_dp*(xF(j)+xD(j))
        xsq(hq,Qh,hg,lh,j)=cf/2._dp*(xn*cf*xy(j)+xn**2*xxy(j)
     &                            +(xn**2-2._dp)/2._dp*xz(j))
        enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
