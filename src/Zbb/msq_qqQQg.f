      subroutine msq_qqQQg(i1,i2,i3,i4,i5,i6,i7,xsq)
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
      integer i1,i2,i3,i4,i5,i6,i7,j,hq,Qh,hg,lh
      double complex a1(2,2,2,2),a2(2,2,2,2),a3(2,2,2,2),a4(2,2,2,2),
     .               b1(2,2,2,2),b2(2,2,2,2),b3(2,2,2,2),b4(2,2,2,2)
      double complex temp
      double precision xA(4),xD(4),xF(4),xG(4)
      double precision xy(4),xz(4),xxy(4),xsq(2,2,2,2,4)
      
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
c        mA=cdabs(mbar1)**2+cdabs(mbar2)**2
c     .    +cdabs(mbar3)**2+cdabs(mbar4)**2
        xA(1)=cdabs(a1(hq,Qh,hg,lh))**2+cdabs(a2(hq,Qh,hg,lh))**2
     .       +cdabs(a3(hq,Qh,hg,lh))**2+cdabs(a4(hq,Qh,hg,lh))**2
        xA(2)=cdabs(b1(Qh,hq,hg,lh))**2+cdabs(b2(Qh,hq,hg,lh))**2
     .       +cdabs(b3(Qh,hq,hg,lh))**2+cdabs(b4(Qh,hq,hg,lh))**2
        temp =2d0*(a1(hq,Qh,hg,lh)*Dconjg(b3(Qh,hq,hg,lh))
     .            +a2(hq,Qh,hg,lh)*Dconjg(b4(Qh,hq,hg,lh))
     .            +a3(hq,Qh,hg,lh)*Dconjg(b1(Qh,hq,hg,lh))
     .            +a4(hq,Qh,hg,lh)*Dconjg(b2(Qh,hq,hg,lh)))
        xA(3)=dble(temp)
        xA(4)=dimag(temp)
c--- eqn (B51) - vanishes because different quarks
c        mB=0d0
c--- eqn (B52) - vanishes because different quarks
c        mC=0d0
c--- eqn (B53)
c        mD=2d0*dble(mbar1*Dconjg(mbar2)+mbar3*Dconjg(mbar4))
        xD(1)=2d0*dble(a1(hq,Qh,hg,lh)*Dconjg(a2(hq,Qh,hg,lh))
     .                 +a3(hq,Qh,hg,lh)*Dconjg(a4(hq,Qh,hg,lh)))
        xD(2)=2d0*dble(b1(Qh,hq,hg,lh)*Dconjg(b2(Qh,hq,hg,lh))
     .                 +b3(Qh,hq,hg,lh)*Dconjg(b4(Qh,hq,hg,lh)))
        temp =2d0*(a1(hq,Qh,hg,lh)*Dconjg(b4(Qh,hq,hg,lh))
     .            +Dconjg(b3(Qh,hq,hg,lh))*a2(hq,Qh,hg,lh)
     .            +a3(hq,Qh,hg,lh)*Dconjg(b2(Qh,hq,hg,lh))
     .            +Dconjg(b1(Qh,hq,hg,lh))*a4(hq,Qh,hg,lh))
        xD(3)=dble(temp)
        xD(4)=dimag(temp)
c--- eqn (B54) - vanishes because different quarks
c        mE=0d0
c--- eqn (B55)
c        mF=2d0*dble(mbar1*Dconjg(mbar3)+mbar2*Dconjg(mbar4))
        xF(1)=2d0*dble(a1(hq,Qh,hg,lh)*Dconjg(a3(hq,Qh,hg,lh))
     .                 +a2(hq,Qh,hg,lh)*Dconjg(a4(hq,Qh,hg,lh)))
        xF(2)=2d0*dble(b1(Qh,hq,hg,lh)*Dconjg(b3(Qh,hq,hg,lh))
     .                 +b2(Qh,hq,hg,lh)*Dconjg(b4(Qh,hq,hg,lh)))
        temp =2d0*(a1(hq,Qh,hg,lh)*Dconjg(b1(Qh,hq,hg,lh))
     .            +Dconjg(b3(Qh,hq,hg,lh))*a3(hq,Qh,hg,lh)
     .            +a2(hq,Qh,hg,lh)*Dconjg(b2(Qh,hq,hg,lh))
     .            +Dconjg(b4(Qh,hq,hg,lh))*a4(hq,Qh,hg,lh))
        xF(3)=dble(temp)
        xF(4)=dimag(temp)
c--- eqn (B56)
c        mG=2d0*dble(mbar1*Dconjg(mbar4)+mbar2*Dconjg(mbar3))
        xG(1)=2d0*dble(a1(hq,Qh,hg,lh)*Dconjg(a4(hq,Qh,hg,lh))
     .                 +a2(hq,Qh,hg,lh)*Dconjg(a3(hq,Qh,hg,lh)))
        xG(2)=2d0*dble(b1(Qh,hq,hg,lh)*Dconjg(b4(Qh,hq,hg,lh))
     .                 +b2(Qh,hq,hg,lh)*Dconjg(b3(Qh,hq,hg,lh)))
        temp=2d0*(a1(hq,Qh,hg,lh)*Dconjg(b2(Qh,hq,hg,lh))
     .           +Dconjg(b3(Qh,hq,hg,lh))*a4(hq,Qh,hg,lh)
     .           +a2(hq,Qh,hg,lh)*Dconjg(b1(Qh,hq,hg,lh))
     .           +Dconjg(b4(Qh,hq,hg,lh))*a3(hq,Qh,hg,lh))
        xG(3)=dble(temp)
        xG(4)=dimag(temp)
c--- eqns (B44)-(B49)
c        m0=mB+mC+mE
c        mx=-0.5d0*(3d0*mC+2d0*mE+mB)
c        my=mA+mD
c        mz=mF+mG
c        mxx=0.25d0*(2d0*mC+mE)
c        mxy=-0.5d0*(mF+mD)
c--- eqn (B42) translated
c--- note that C3 = (N*CF**2-Tr*CF)/2 = CF/4*(N**2-2)
c--- in my notation, where Tr=1/2
c        msq(hq,Qh,hg,lh)=cf/2d0*(xn*cf*my+xn**2*mxy+(xn**2-2d0)/2d0*mz)
        do j=1,4
        xy(j)=xA(j)+xD(j)
        xz(j)=xF(j)+xG(j)
        xxy(j)=-0.5d0*(xF(j)+xD(j))
        xsq(hq,Qh,hg,lh,j)=cf/2d0*(xn*cf*xy(j)+xn**2*xxy(j)
     .                            +(xn**2-2d0)/2d0*xz(j))
        enddo
      enddo
      enddo
      enddo
      enddo

      return
      end
