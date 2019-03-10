      subroutine Cfill_recur3(p1,p2,p1p2,m1,m2,m3,N0)
      implicit none
C     Implements the calculation of the formfactors
C     for small momenta, as in DD Eq.5.64-5.70 etc
C     N0 is the offset in the common block

C--- Currently: calculates up to rank 4 with at least one recursion
c---            calculates rank 5 with no recursion
c---            calculates C00iiii components of rank 6

      include 'TRconstants.f'
      include 'pvBnames.f'
      include 'pvBv.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      integer B12,B23,B13,np,ep,N0,pvBcache,
     , j,k,i1,i2,i3,i4,i5,step,kmax
      parameter(np=2)
      double precision p1,p2,p1p2,m1,m2,m3,f(np),Gr(np,np),DetGr
      double complex S0000(-2:0),S0000i(np,-2:0),
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
      logical,save:: first=.true.
!$omp threadprivate(first)

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
      write(6,*) 'small P: 2x2 DetGr = ',DetGr

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


c--- find the largest f(k) for recursion relations
      kmax=1
      do k=2,np
      if (abs(f(k)) .ge. abs(f(kmax))) kmax=k
      enddo

      write(6,*) 'f(kmax) =',f(kmax)
      
C----Begin the iteration scheme

C set all the Cv to zero
      do ep=-2,0
      do j=1,Ncc
      Cv(j+N0,ep)=czip
      enddo
      enddo

      do step=0,5
      if (step .eq. 5) goto 105
      if (step .eq. 4) goto 104
      if (step .eq. 3) goto 103
      if (step .eq. 2) goto 102
      if (step .eq. 1) goto 101
      if (step .eq. 0) goto 100

C----step 5: calculate C00iiii, Ciiiii 
 105  continue

C--- Fixes C00iiii using extension of 5.69
C--- knowing Ciiii, with corrections of order Gr(i,j)*Ciiiiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      call runCP_00iiii(i1,i2,i3,i4,m1,Gr,Bzero4,N0)
      enddo
      enddo
      enddo
      enddo

C--- Fixes Ciiiii using extension of 5.70
c--- knowing C00iiii, with a correction of order Gr(i,j)*Ciiiiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      do i5=i4,np
      call runCP_iiiii(kmax,i1,i2,i3,i4,i5,f,Gr,Shat6,N0)
      enddo
      enddo
      enddo
      enddo
      enddo

C----step 4: calculate C00iii, Ciiii, C0000i
 104  continue

C--- Fixes C00iii using extension of 5.69
C--- knowing Ciii, with corrections of order Gr(i,j)*Ciiiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      call runCP_00iii(i1,i2,i3,m1,Gr,Bzero3,N0)
      enddo
      enddo
      enddo

C--- Fixes Ciiii using extension of 5.70
c--- knowing C00iii, with a correction of order Gr(i,j)*Ciiiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      call runCP_iiii(kmax,i1,i2,i3,i4,f,Gr,Shat5,N0)
      enddo
      enddo
      enddo
      enddo

c--- calculate S0000i (needs C00i) - required for C0000i
      include 'S0000iC_def.f'

C--- Fixes C0000i, with corrections of order Gr(i,j)*C00iii
      do i1=1,np
      call runCP_0000i(i1,Gr,S0000i,N0)
      enddo
      
C----step 3: calculate C00ii, Ciii, C0000 
 103  continue

C--- Fixes C00ii using 5.69
C--- knowing Cii, with corrections of order Gr(i,j)*Ciiii
      do i1=1,np
      do i2=i1,np
      call runCP_00ii(i1,i2,m1,Gr,Bzero2,N0)
      enddo
      enddo

C--- Fixes Ciii using 5.70
c--- knowing C00ii, with a correction of order Gr(i,j)*Ciiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      call runCP_iii(kmax,i1,i2,i3,f,Gr,Shat4,N0)
      enddo
      enddo
      enddo

c--- calculate S0000 (needs C00) - required for C0000
      include 'S0000C_def.f'

C--- Fixes C0000, with corrections of order Gr(i,j)*C00ii
      call runCP_0000(Gr,S0000,N0)
     
C---- step 2: calculate C00i and Cii
 102  continue
C--- Fixes C00i according to 5.67 Denner-Dittmaier
C--- knowing Ci, with correction of order Gr(i,j)*Ciii
      do i1=1,np
      call runCP_00i(i1,m1,Gr,Bzero1,N0)
      enddo

C--- Fixes Cii using Eq. 5.68 Denner-Dittmaier
C--- knowing C00i with correction of order Gr(i,j)*Ciii
      do i1=1,np
      do i2=i1,np
      call runCP_ii(kmax,i1,i2,f,Gr,Shat3,N0)
      enddo
      enddo

C----step 1: calculate C00 and Ci
 101  continue
C--- Fixes C00 according to 5.65 Denner-Dittmaier
C--- knowing C0 with corrections of order Gr(i,j)*Cii  
      call runCP_00(m1,Gr,Bzero0,N0)

C--- Fixes Ci according to 5.66 Denner-Dittmaier
C--- knowing C00 with corrections of order Gr(i,j)*Cii
      do i1=1,np
      call runCP_i(kmax,i1,f,Gr,Shat2,N0)
      enddo

C----step 0: calculate C0
 100  continue
C--- Fixes C0 according to 5.64 Denner-Dittmaier
C--- with corrections of order Gr(i,j)*Ci
      call runCP_0(kmax,f,Gr,Shat1,N0)

      enddo

    
c--- check the contents of triangle array    
c      write(6,*) 'C array'
c      do ip=1,Ncc
c      if (abs(Csing(ip,p1p2,p1,p2,m1,m2,m3)) .ne. 0d0) then
c      write(6,'(i3,2f20.15)') ip,Cv(ip+N0,-1)
c     .                           /Csing(ip,p1p2,p1,p2,m1,m2,m3)
c      endif
c      enddo
c      pause
      
c--- check the contents of bubble arrays    
c      write(6,*) 'B12 array'
c      do ip=1,Nbb
c      write(6,'(i3,2f20.15)') ip,Bv(ip+B12,-1)/Bsing(ip,p1,m1,m2)
c      enddo
    
c      write(6,*) 'B13 array'
c      do ip=1,Nbb
c      write(6,'(i3,2f20.15)') ip,Bv(ip+B13,-1)/Bsing(ip,p1p2,m1,m3)
c      enddo
    
c      write(6,*) 'B23 array'
c      do ip=1,Nbb
c      write(6,'(i3,2f20.15)') ip,Bv(ip+B23,-1)/Bsing(ip,p2,m2,m3)
c      enddo
    
c      pause

      end


