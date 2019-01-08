      subroutine Cfill_recur2(p1,p2,p1p2,m1,m2,m3,N0,exceptional)
      implicit none
C     Implements the calculation of the formfactors
C     for small Gram Determinant and small Y, as in DD Eq.5.54-5.61 etc
C     N0 is the offset in the common block

C--- Currently: calculates up to rank 3 with at least one recursion
c---            calculates ranks 4 and 5 with no recursion
c---            calculates metric tensor components of ranks 6 and 7

c--- JC: 11/22/2012 added an extra level of recursion. No additional
c---     identities are used, but the extra loop improves the
c---     numerical precision

      include 'TRconstants.f'
      include 'pvBnames.f'
      include 'pvBv.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'pvverbose.f'
      integer B12,B23,B13,np,ep,N0,pvBcache,
     , i,j,k,l,m,n,i1,i2,i3,i4,i5,step,jx,kgt,lgt,ixt,jxt
      parameter(np=2)
      double precision p1,p2,p1p2,m1,m2,m3,f(np),
     . Gtwiddle(np,np),Xtwiddle0(np),Gr(np,np),DetGr,Gtt(np,np,np,np),
     . Xtwiddle(0:np,0:np),Xtmax
      double complex
     . Shat3zz(np,-2:0),Shat4zz(np,z1max,-2:0),
     . Shat5zz(np,z2max,-2:0),Shat6zz(np,z3max,-2:0),
     . Shat5zzzz(np,-2:0),Shat6zzzz(np,z1max,-2:0),
     . Shat7zz(np,z4max,-2:0),Shat7zzzz(np,z2max,-2:0),
     . Shat1(np,-2:0),Shat2(np,z1max,-2:0),
     . Shat3(np,z2max,-2:0),Shat4(np,z3max,-2:0),
     . Shat5(np,z4max,-2:0),Shat6(np,z5max,-2:0),Shat7(np,z6max,-2:0)
      double complex bsum1(-2:0),
     . bsum0(-2:0),bsum11(-2:0),bsum00(-2:0),
     . bsum111(-2:0),bsum1111(-2:0),bsum001(-2:0),
     . bsum0011(-2:0),bsum0000(-2:0),
     . bsum00001(-2:0),bsum00111(-2:0),bsum11111(-2:0),
     . Bzero5(z5max,-2:0),Bzero4(z4max,-2:0),
     . Bzero3(z3max,-2:0),Bzero2(z2max,-2:0),Bzero1(z1max,-2:0),
     . Bzero0(-2:0)
      logical exceptional
      logical,save :: first=.true.
!$omp threadprivate(first)

      exceptional=.false.

      if (first) then
        first=.false.
        call Array2dim
        call CArraysetup
      endif

c--- Not necessary, routine upgraded now
c      if ((m1 .ne. 0d0).or.(m2 .ne. 0d0).or.(m3 .ne. 0d0)) then
c      write(6,*) 'nonzero internal masses'
c      stop
c      endif

      B12=pvBcache(p1,m1,m2)
      B23=pvBcache(p2,m2,m3)
      B13=pvBcache(p1p2,m1,m3)

C----We have changed the sign of fi (different from Dfill) to agree
C----with notation of Denner-Dittmaier
      f(1) = -m2 + m1 + p1
      f(2) = -m3 + m1 + p1p2

      Gr(1,1)=2*p1
      Gr(2,2)=2*p1p2
      Gr(1,2)=p1+p1p2-p2
      Gr(2,1)=Gr(1,2)

      call determinant(2,np,Gr,DetGr)
      if (pvverbose) write(6,*) 'small Y: 2x2 DetGr = ',DetGr

      Gtwiddle(1,1)=Gr(2,2)
      Gtwiddle(2,2)=Gr(1,1)
      Gtwiddle(1,2)=-Gr(1,2)
      Gtwiddle(2,1)=-Gr(2,1)

C----setup Gtt
      do i=1,2
      do k=1,2
      do j=1,2
      do l=1,2
      Gtt(i,k,j,l)=delta(i,l)*delta(k,j)-delta(i,j)*delta(k,l)
      enddo
      enddo
      enddo
      enddo

c--- setup Xtwiddle
      do j=1,2
      Xtwiddle0(j)=-Gtwiddle(j,1)*f(1)-Gtwiddle(j,2)*f(2)
      Xtwiddle(0,j)=Xtwiddle0(j)
      Xtwiddle(j,0)=Xtwiddle(0,j)
      enddo
      Xtwiddle(0,0)=DetGr

      do i=1,2
      do j=1,2
      Xtwiddle(i,j)=2d0*m1*Gtwiddle(i,j)
      do n=1,2
      do m=1,2
      Xtwiddle(i,j)=Xtwiddle(i,j)+Gtt(i,n,j,m)*f(n)*f(m)
      enddo
      enddo
      enddo
      enddo


      do ep=-2,0
      Bsum0(ep)=Bv(bb0+B23,ep)+Bv(bb1+B23,ep)
      Bsum1(ep)=Bv(bb1+B23,ep)+Bv(bb11+B23,ep)
      Bsum00(ep)=Bv(bb00+B23,ep)+Bv(bb001+B23,ep)
      Bsum11(ep)=Bv(bb11+B23,ep)+Bv(bb111+B23,ep)
      Bsum001(ep)=Bv(bb001+B23,ep)+Bv(bb0011+B23,ep)
      Bsum111(ep)=Bv(bb111+B23,ep)+Bv(bb1111+B23,ep)
      Bsum0000(ep)=Bv(bb0000+B23,ep)+Bv(bb00001+B23,ep)
      Bsum0011(ep)=Bv(bb0011+B23,ep)+Bv(bb00111+B23,ep)
      Bsum1111(ep)=Bv(bb1111+B23,ep)+Bv(bb11111+B23,ep)
      Bsum00001(ep)=Bv(bb00001+B23,ep)+Bv(bb000011+B23,ep)
      Bsum00111(ep)=Bv(bb00111+B23,ep)+Bv(bb001111+B23,ep)
      Bsum11111(ep)=Bv(bb11111+B23,ep)+Bv(bb111111+B23,ep)
      enddo

c      write(6,'(a9,2f20.15)') 'Bsum0',Bsum0(-1)
c      write(6,'(a9,2f20.15)') 'Bsum1',Bsum1(-1)
c      write(6,'(a9,2f20.15)') 'Bsum00',Bsum00(-1)
c      write(6,'(a9,2f20.15)') 'Bsum11',Bsum11(-1)
c      write(6,'(a9,2f20.15)') 'Bsum001',Bsum001(-1)
c      write(6,'(a9,2f20.15)') 'Bsum111',Bsum111(-1)
c      write(6,'(a9,2f20.15)') 'Bsum0000',Bsum0000(-1)
c      write(6,'(a9,2f20.15)') 'Bsum0011',Bsum0011(-1)
c      write(6,'(a9,2f20.15)') 'Bsum1111',Bsum1111(-1)

c--- new implementation, in the same style as Dfill_alt.f
c---  (except ShatC.f also includes the zz definitions)
      do ep=-2,0
      include 'ShatC.f'
      enddo


c--- note: these are the triangle parts of the S00 functions that
c---       are defined above (and commented out), except that these
c---       are a factor of two smaller

      do ep=-2,0
      include 'Bzero.f'
      enddo


      jx=1
      do j=2,np
      if (abs(Xtwiddle0(j)) .ge. abs(Xtwiddle0(jx))) jx=j
      enddo

      kgt=1
      lgt=1
      do k=1,np
      do l=k,np
      if (abs(Gtwiddle(k,l)) .ge. abs(Gtwiddle(kgt,lgt))) then
      kgt=k
      lgt=l
      endif
      enddo
      enddo

      ixt=1
      jxt=1
      do i=1,np
      do j=i,np
      if (abs(Xtwiddle(i,j)) .ge. abs(Xtwiddle(ixt,jxt))) then
      ixt=i
      jxt=j
      endif
      enddo
      enddo
      
c      do k=1,np
c      do l=k,np
c      write(6,*) k,l,Gtwiddle(k,l),Gtwiddle(kgt,lgt)
c      enddo
c      enddo
c      pause
c      write(6,*) '   Xtwiddle(0,0)', Xtwiddle(0,0)
c      do j=1,np
c      write(6,*) 'j, Xtwiddle(0,j)',j, Xtwiddle(0,j)
c      enddo
c      do j=1,np
c      do k=1,np
c      write(6,*) 'j, k, Xtwiddle(j,k)',j,k, Xtwiddle(j,k)
c      enddo
c      enddo
      
c--- calculate maximum entry in Xtwiddle
      Xtmax=abs(Gr(1,1))
      do j=1,np
      do k=1,np
      Xtmax=max(abs(Xtwiddle(j,k)),Xtmax)
      enddo
      enddo
c      write(6,*) 'Xtmax=',Xtmax

c--- check for exceptional case where none of the Xtwiddle(i,j) elements
c---  are large compared to DetGr (Xtwiddle(0,0)) or Xtwiddle(0,j)
c---  [see note at end of Sec.5.5 in DD].
      if ( (Xtmax/abs(Xtwiddle(0,0)) .lt. 1d1)
     . .or.(Xtmax/abs(Xtwiddle0(jx)) .lt. 1d1)) then
        if (pvverbose) then
          write(6,*) 'EXCEPTIONAL CASE'
        write(6,*) 'Maximum Xtwiddle(i,j) = ',Xtmax
        write(6,*) '        Xtwiddle(0,0) = ',Xtwiddle(0,0)
        write(6,*) 'maximum Xtwiddle(0,j) = ',Xtwiddle(0,jx)
      endif
      exceptional=.true.
      return
      endif
      
C----Begin the iteration scheme

C set all the Cv to zero
      do ep=-2,0
      do j=1,Ncc
      Cv(j+N0,ep)=czip
      enddo
      enddo

      do step=0,3
      if (step .eq. 3) goto 103
      if (step .eq. 2) goto 102
      if (step .eq. 1) goto 101
      if (step .eq. 0) goto 100

C--- step 3
 103  continue

C--- step 2: calculate C00iiii, C00iiiii, Ciiii, Ciiiii,
c---                   C0000ii, C0000iii, C000000,C000000i
 102  continue

C--- a) Calculate C00iiii
C---    Small terms of order Xtwiddle(0,k)*Ciiiii,Xtwiddle(0,0)*Ciiiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00llll(kgt,lgt,Xtwiddle,Gtwiddle,Shat6,N0)

      do i1=1,np
      if (i1 .ne. lgt) then
C---    Calculate C00llli, requires C00llll
C---    Small terms of order Xtwiddle(0,k)*Ciiiii,Xtwiddle(0,0)*Ciiiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00llli(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat6,N0)
      endif
      enddo

      do i1=1,np
      do i2=i1,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt)) then
C---    Calculate C00lli1i2, requires C00llli1
C---    Small terms of order Xtwiddle(0,k)*Ciiiii,Xtwiddle(0,0)*Ciiiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00lli1i2(kgt,lgt,i1,i2,Xtwiddle,Gtwiddle,Shat6,N0)
      endif
      enddo
      enddo

      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt) .and. (i3 .ne. lgt)) then
C---    Calculate C00li1i2i3, requires C00lli1i2
C---    Small terms of order Xtwiddle(0,k)*Ciiiii,Xtwiddle(0,0)*Ciiiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00li1i2i3(kgt,lgt,i1,i2,i3,Xtwiddle,Gtwiddle,Shat6,N0)
      endif
      enddo
      enddo
      enddo

      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      if (   (i1 .ne. lgt) .and. (i2 .ne. lgt)
     . .and. (i3 .ne. lgt) .and. (i4 .ne. lgt)) then
C---    Calculate C00i1i2i3i4, requires C00li1i2i3
C---    Small terms of order Xtwiddle(0,k)*Ciiiii,Xtwiddle(0,0)*Ciiiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00i1i2i3i4(kgt,lgt,i1,i2,i3,i4,
     .                     Xtwiddle,Gtwiddle,Shat6,N0)
      endif
      enddo
      enddo
      enddo
      enddo

C--- b) Calculate C00iiiii
C---    Small terms of order Xtwiddle(0,k)*Ciiiiii,Xtwiddle(0,0)*Ciiiiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00lllll(kgt,lgt,Xtwiddle,Gtwiddle,Shat7,N0)

      do i1=1,np
      if (i1 .ne. lgt) then
C---    Calculate C00lllli, requires C00lllll
C---    Small terms of order Xtwiddle(0,k)*Ciiiiii,Xtwiddle(0,0)*Ciiiiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00lllli(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat7,N0)
      endif
      enddo

      do i1=1,np
      do i2=i1,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt)) then
C---    Calculate C00llli1i2, requires C00lllli1
C---    Small terms of order Xtwiddle(0,k)*Ciiiiii,Xtwiddle(0,0)*Ciiiiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00llli1i2(kgt,lgt,i1,i2,Xtwiddle,Gtwiddle,Shat7,N0)
      endif
      enddo
      enddo

      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt) .and. (i3 .ne. lgt)) then
C---    Calculate C00lli1i2i3, requires C00llli1i2
C---    Small terms of order Xtwiddle(0,k)*Ciiiiii,Xtwiddle(0,0)*Ciiiiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00lli1i2i3(kgt,lgt,i1,i2,i3,Xtwiddle,Gtwiddle,Shat7,N0)
      endif
      enddo
      enddo
      enddo

      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      if (   (i1 .ne. lgt) .and. (i2 .ne. lgt)
     . .and. (i3 .ne. lgt) .and. (i4 .ne. lgt)) then
C---    Calculate C00li1i2i3i4, requires C00lli1i2i3
C---    Small terms of order Xtwiddle(0,k)*Ciiiiii,Xtwiddle(0,0)*Ciiiiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00li1i2i3i4(kgt,lgt,i1,i2,i3,i4,
     .                      Xtwiddle,Gtwiddle,Shat7,N0)
      endif
      enddo
      enddo
      enddo
      enddo
      
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      do i5=i4,np
      if (   (i1 .ne. lgt) .and. (i2 .ne. lgt)
     . .and. (i3 .ne. lgt) .and. (i4 .ne. lgt) .and. (i5 .ne. lgt)) then
C---    Calculate C00i1i2i3i4i5, requires C00li1i2i3i4
C---    Small terms of order Xtwiddle(0,k)*Ciiiiii,Xtwiddle(0,0)*Ciiiiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00i1i2i3i4i5(kgt,lgt,i1,i2,i3,i4,i5,
     .                       Xtwiddle,Gtwiddle,Shat7,N0)
      endif
      enddo
      enddo
      enddo
      enddo
      enddo

C--- c) Calculate Ciiiii, requires C00iiii,C00iiiii
C---    Small terms of order Xtwiddle(0,j)*Ciiiiii
C---    Denominator Xtwiddle(i,j)
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      do i5=i4,np
      call runCY_i1i2i3i4i5(ixt,jxt,i1,i2,i3,i4,i5,
     .                     f,Xtwiddle,Gtt,Gtwiddle,Shat6,Bzero5,N0)
      enddo
      enddo
      enddo
      enddo
      enddo

C--- d) Calculate Ciiii, requires C00iii,C00iiii
C---    Small terms of order Xtwiddle(0,j)*Ciiiii
C---    Denominator Xtwiddle(i,j)
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      call runCY_i1i2i3i4(ixt,jxt,i1,i2,i3,i4,
     . f,Xtwiddle,Gtt,Gtwiddle,Shat5,Bzero4,N0)
      enddo
      enddo
      enddo
      enddo

C--- e) Calculate C0000ii
C---    Small terms of order Xtwiddle(0,k)*Czziii,Xtwiddle(0,0)*Czziiii
C---    Denominator Gtwiddle(k,l)
      call runCY_0000ll(kgt,lgt,Xtwiddle,Gtwiddle,Shat6zz,N0)

      do i1=1,np
      if (i1 .ne. lgt) then
C---    Calculate C0000li, requires C0000ll
C---    Small terms of order Xtwiddle(0,k)*C00iii,Xtwiddle(0,0)*C00iiii
C---    Denominator Gtwiddle(k,l)
      call runCY_0000li(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat6zz,N0)
      endif 
      enddo

      do i1=1,np
      do i2=i1,np
      if ((i1.ne.lgt) .and. (i2 .ne. lgt)) then  
C---    Calculate C0000i1i2, requires C0000li1,C0000li2
C---    Small terms of order Xtwiddle(0,k)*C00iii,Xtwiddle(0,0)*C00iiii
C---    Denominator Gtwiddle(k,l)
      call runCY_0000i1i2(kgt,lgt,i1,i2,Xtwiddle,Gtwiddle,Shat6zz,N0)
      endif
      enddo
      enddo
      
C--- f) Calculate C0000iii
C---    Small terms of order Xtwiddle(0,k)*Czziiii,Xtwiddle(0,0)*C00iiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_0000lll(kgt,lgt,Xtwiddle,Gtwiddle,Shat7zz,N0)

      do i1=1,np
      if (i1 .ne. lgt) then
C---    Calculate C0000lli, requires C0000lll
C---    Small terms of order Xtwiddle(0,k)*C00iiii,Xtwiddle(0,0)*C00iiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_0000lli(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat7zz,N0)
      endif
      enddo

      do i1=1,np
      do i2=i1,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt)) then
C---    Calculate C0000li1i2, requires C0000lli1
C---    Small terms of order Xtwiddle(0,k)*C00iiii,Xtwiddle(0,0)*C00iiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_0000li1i2(kgt,lgt,i1,i2,Xtwiddle,Gtwiddle,Shat7zz,N0)
      endif
      enddo
      enddo

C---    Calculate C0000i1i2i3, requires C0000li1i2,C0000li2i3,C0000li3i1
C---    Small terms of order Xtwiddle(0,k)*C00iiii,Xtwiddle(0,0)*C00iiiii
C---    Denominator Gtwiddle(k,l)
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      if ((i1 .ne. lgt) .and.(i2 .ne. lgt) .and.(i3 .ne. lgt)) then
      call runCY_0000i1i2i3(kgt,lgt,i1,i2,i3,
     .                     Xtwiddle,Gtwiddle,Shat7zz,N0)
      endif
      enddo
      enddo
      enddo

C--- g) Calculate C000000
C---    Small terms of order Xtwiddle(0,k)*C0000i,Xtwiddle(0,0)*C0000ii
C---    Denominator Gtwiddle(k,l)
      call runCY_000000(kgt,lgt,Xtwiddle,Gtwiddle,Shat6zzzz,N0)

C--- h) Calculate C000000i
C---    Small terms of order Xtwiddle(0,k)*C0000ii,Xtwiddle(0,0)*C0000iii
C---    Denominator Gtwiddle(k,l)
      call runCY_000000l(kgt,lgt,Xtwiddle,Gtwiddle,Shat7zzzz,N0)

      do i1=1,np
      if (i1 .ne. lgt)
     . call runCY_000000i(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat7zzzz,N0)
      enddo

C--- step 1: calculate C00ii, C00iii, Cii, Ciii, C0000, C0000i
 101  continue

C--- a) Calculate C00ii
C---    Small terms of order Xtwiddle(0,k)*Ciii,Xtwiddle(0,0)*Ciiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00ll(kgt,lgt,Xtwiddle,Gtwiddle,Shat4,N0)

      do i1=1,np
      if (i1 .ne. lgt) then
C---    Calculate C00li, requires C00ll
C---    Small terms of order Xtwiddle(0,k)*Ciii,Xtwiddle(0,0)*Ciiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00li(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat4,N0)
      endif 
      enddo

      do i1=1,np
      do i2=i1,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt)) then  
C---    Calculate C00i1i2, requires C00li1,C00li2
C---    Small terms of order Xtwiddle(0,k)*Ciii,Xtwiddle(0,0)*Ciiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00i1i2(kgt,lgt,i1,i2,Xtwiddle,Gtwiddle,Shat4,N0)
      endif
      enddo
      enddo
      
c--- b) Calculate C00iii
C---    Small terms of order Xtwiddle(0,k)*Ciiii,Xtwiddle(0,0)*Ciiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00lll(kgt,lgt,Xtwiddle,Gtwiddle,Shat5,N0)

      do i1=1,np
      if (i1 .ne. lgt) then
C---    Calculate C00lli, requires C00lll
C---    Small terms of order Xtwiddle(0,k)*Ciiii,Xtwiddle(0,0)*Ciiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00lli(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat5,N0)
      endif
      enddo

      do i1=1,np
      do i2=i1,np
      if ((i1 .ne. lgt) .and. (i2 .ne. lgt)) then
C---    Calculate C00li1i2, requires C00lli1
C---    Small terms of order Xtwiddle(0,k)*Ciiii,Xtwiddle(0,0)*Ciiiii
C---    Denominator Gtwiddle(k,l)
      call runCY_00li1i2(kgt,lgt,i1,i2,Xtwiddle,Gtwiddle,Shat5,N0)
      endif
      enddo
      enddo

C---    Calculate C00i1i2i3, requires C00li1i2,C00li2i3,C00li3i1
C---    Small terms of order Xtwiddle(0,k)*Ciiii,Xtwiddle(0,0)*Ciiiii
C---    Denominator Gtwiddle(k,l)
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      if ((i1 .ne. lgt) .and.(i2 .ne. lgt) .and.(i3 .ne. lgt)) then
      call runCY_00i1i2i3(kgt,lgt,i1,i2,i3,Xtwiddle,Gtwiddle,Shat5,N0)
      endif
      enddo
      enddo
      enddo

C--- c) Calculate Ciii, requires C00ii,C00iii
C---    Small terms of order Xtwiddle(0,j)*Ciiii
C---    Denominator Xtwiddle(i,j)
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      call runCY_i1i2i3(ixt,jxt,i1,i2,i3,f,Xtwiddle,Gtt,Gtwiddle,Shat4,
     . Bzero3,N0)
      enddo
      enddo
      enddo

C--- d) Calculate Cii, requires C00i,C00ii
C---    Small terms of order Xtwiddle(0,j)*Ciii
C---    Denominator Xtwiddle(i,j)
      do i1=1,np
      do i2=i1,np
      call runCY_i1i2(ixt,jxt,i1,i2,f,Xtwiddle,Gtt,Gtwiddle,Shat3,
     . Bzero2,N0)
      enddo
      enddo

C--- e) Calculate C0000
C---    Small terms of order Xtwiddle(0,k)*C00i,Xtwiddle(0,0)*C00ii
C---    Denominator Gtwiddle(k,l)
      call runCY_0000(kgt,lgt,Xtwiddle,Gtwiddle,Shat4zz,N0)

C--- f) Calculate C0000l
C---    Small terms of order Xtwiddle(0,k)*C00ii,Xtwiddle(0,0)*C00iii
C---    Denominator Gtwiddle(k,l)
      call runCY_0000l(kgt,lgt,Xtwiddle,Gtwiddle,Shat5zz,N0)

      do i1=1,np
      if (i1 .ne. lgt)
     . call runCY_0000i(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat5zz,N0)
      enddo


C--- step 0: calculate C00,C00i,C0 and Ci
 100  continue
C--- a) Calculate C00
C---    Small terms of order Xtwiddle(0,k)*Ci,Xtwiddle(0,0)*Cii
C---    Denominator Gtwiddle(kgt,lgt)
      call runCY_00(kgt,lgt,Xtwiddle,Gtwiddle,Shat2,N0)

C--- b) Calculate C00l, C00i
C---    Small terms of order Xtwiddle(0,k)*Cii,Xtwiddle(0,0)*Ciii
C---    Denominator Gtwiddle(k,l)
      call runCY_00l(kgt,lgt,Xtwiddle,Gtwiddle,Shat3,N0)
C---    Calculate C00i1, requires C00l
C---    Small terms of order Xtwiddle(0,k)*Dli1,Xtwiddle(0,0)*Dkli1
C---    Denominator Gtwiddle(k,l)
      do i1=1,np
      if (i1 .ne. lgt) then
        call runCY_00i(kgt,lgt,i1,Xtwiddle,Gtwiddle,Shat3,N0)
      endif
      enddo

C--- c) Calculate Ci, requires C00i
C---    Small terms of order Xtwiddle(0,j)*Cii
C---    Denominator Xtwiddle(i,j)
      do i1=1,np
      call runCY_i(ixt,jxt,i1,f,Xtwiddle,Gtt,Gtwiddle,Shat2,Bzero1,N0)
      enddo

C--- d) Calculates C0
C---    Requires C00, small terms of order Xtwiddle(0,j)*Ci
C---    Denominator Xtwiddle(i,j)
      call runCY_0(ixt,jxt,f,Xtwiddle,Gtwiddle,Gtt,Shat1,Bzero0,N0)
     
c--- check the contents of triangle array    
c      write(6,*) 'recur2: C array'
c      do ip=1,13
c      write(6,'(i3,2e20.12)') ip,Cv(ip+N0,-1)
c      enddo
c      pause

      enddo

    
c--- check the contents of triangle array    
c      write(6,*) 'C array'
c      write(6,*) p1p2,p1,p2,m1,m2,m3
c      do ip=1,Ncc
c      if (abs(Csing(ip,p1p2,p1,p2,m1,m2,m3)) .ne. 0d0) then
c      write(6,'(i3,4f20.15)') ip,Cv(ip+N0,-1),Cv(ip+N0,-1)
c     .                    /Csing(ip,p1p2,p1,p2,m1,m2,m3)
c      endif
c      enddo
c      pause
      
c--- check the contents of bubble arrays    
c      write(6,*) 'B12 array'
c      do ip=1,Nbb
c      if (abs(Bsing(ip,p1,m1,m2)) .ne. 0d0) then
c      write(6,'(i3,2f20.15)') ip,Bv(ip+B12,-1)/Bsing(ip,p1,m1,m2)
c      endif
c      enddo
   
c      write(6,*) 'B13 array'
c      do ip=1,Nbb
c      if (abs(Bsing(ip,p1p2,m1,m3)) .ne. 0d0) then
c      write(6,'(i3,2f20.15)') ip,Bv(ip+B13,-1)/Bsing(ip,p1p2,m1,m3)
c      endif
c      enddo
   
c      write(6,*) 'B23 array'
c      do ip=1,Nbb
c      if (abs(Bsing(ip,p2,m2,m3)) .ne. 0d0) then
c      write(6,'(i3,2f20.15)') ip,Bv(ip+B23,-1)/Bsing(ip,p2,m2,m3)
c      endif
c      enddo
   
c      pause

c   77 format(a3,i2,a5,3('(',e13.6,',',e13.6,') '))
    
      end
