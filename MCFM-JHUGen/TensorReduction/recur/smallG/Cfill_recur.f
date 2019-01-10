      subroutine Cfill_recur(p1,p2,p1p2,m1,m2,m3,N0)
      implicit none
C     Implements the calculation of the formfactors
C     for small Gram Determinant as in DD Eq.5.41-5.48 etc
C     N0 is the offset in the common block

C--- Currently: calculates up to rank 4 with at least one recursion
c---            calculates rank 5 with no recursion
c---            calculates metric tensor components of rank 6

      include 'TRconstants.f'
      include 'pvBnames.f'
      include 'pvBv.f'
      include 'pvCnames.f'
      include 'pvCv.f'
      include 'Carraydef.f'
      include 'pvverbose.f'
      integer B12,B23,B13,np,ep,N0,pvBcache,
     , i,j,k,l,i1,i2,i3,i4,i5,step,jmax,kmax,lmax
      parameter(np=2)
      double precision p1,p2,p1p2,m1,m2,m3,f(np),
     . Gtwiddle(np,np),Xtwiddle0(np),Gr(np,np),DetGr,Gtt(np,np,np,np)
      double complex S00(-2:0),S0000(-2:0),S000000(-2:0),
     . S0000i(np,-2:0),S0000ii(z2max,-2:0),
     . S00i(np,-2:0),S00ii(z2max,-2:0),
     . S00iii(z3max,-2:0),S00iiii(z4max,-2:0),
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
     . bsum00001(-2:0),bsum00111(-2:0),bsum11111(-2:0)
 
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

      GR(1,1)=2*p1
      GR(2,2)=2*p1p2
      GR(1,2)=p1+p1p2-p2
      GR(2,1)=GR(1,2)

c--- use recursion in singular regions only
c      do i=1,2
c      do j=1,2
c        Grc(i,j)=dcmplx(GR(i,j))
c      enddo
c      enddo      
c      Gsing=pvGramsing(GRc,2)
cc--- note that logic is inverted by local version of PVGramsing
c      if (Gsing .eqv. .true.) then
cc--- not in singular region: use normal PV reduction  
c        call pvCfill(p1,p2,p1p2,m1,m2,m3,N0)    
c        return
c      else
c        write(6,*) 'singular triangle Gram'
c      endif
      
      call determinant(2,np,GR,DetGR)
      Gtwiddle(1,1)=GR(2,2)
      Gtwiddle(2,2)=GR(1,1)
      Gtwiddle(1,2)=-GR(1,2)
      Gtwiddle(2,1)=-GR(2,1)

      if (pvverbose) write(6,*) 'small G: 2x2 DetGr = ',DetGr

      do j=1,2
      Xtwiddle0(j)=-Gtwiddle(j,1)*f(1)-Gtwiddle(j,2)*f(2)
      enddo


      do i=1,2
      do k=1,2
      do j=1,2
      do l=1,2
      Gtt(i,k,j,l)=delta(i,l)*delta(k,j)-delta(i,j)*delta(k,l)
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


c--- These are now calculated in the recursion, which is required when
c--- the first propagator has a non-zero mass
c      do ep=-2,0
c      include 'S00C.f'
c      enddo


      jmax=1
      do j=2,np
      if (abs(Xtwiddle0(j)) .ge. abs(Xtwiddle0(jmax))) jmax=j
      enddo

      kmax=1
      lmax=1
      do k=1,np
      do l=k,np
      if (abs(Gtwiddle(k,l)) .ge. abs(Gtwiddle(kmax,lmax))) then
      kmax=k
      lmax=l
      endif
      enddo
      enddo

C----Begin the iteration scheme

C set all the Cv to zero
      do ep=-2,0
      do j=1,Ncc
      Cv(j+N0,ep)=czip
      enddo
      enddo


      do step=0,5
c      if (step .eq. 6) goto 106
      if (step .eq. 5) goto 105
      if (step .eq. 4) goto 104
      if (step .eq. 3) goto 103
      if (step .eq. 2) goto 102
      if (step .eq. 1) goto 101
      if (step .eq. 0) goto 100

C----step 6 
c 106  continue

cC---Fixes Ciiiii according to extension of Denner-Dittmaier
c---knowing C00iiii with a correction of order Delta*Ciiiiii
c      do i1=1,np
c      do i2=i1,np
c      do i3=i2,np
c      do i4=i3,np
c      do i5=i4,np
c      call runC_iiiii(jmax,i1,i2,i3,i4,i5,DetGr,Xtwiddle0,Gtwiddle,
c     . Shat6,N0)

c      enddo
c      enddo
c      enddo
c      enddo
c      enddo

C---step 5
 105  continue
C--- Fixes C000000 according to extension of Denner-Dittmaier
C--- knowing C0000 with correction of order Delta*C0000ii

c--- calculate S000000 (needs C0000)
      do ep=-2,0
      S000000(ep)=2d0*Bv(bb0000+B23,ep)+2d0*m1*Cv(cc0000+N0,ep)
      enddo

      call runC_000000(kmax,lmax,DetGr,f,Gtwiddle,Gtt,
     . Shat5zzzz,Shat6zzzz,S000000,N0)
C--- Fixes C0000ii according to extension of Denner-Dittmaier 
C--- knowing C00ii,C0000i,C000000 with correction of order Delta*C00iiii

c--- calculate S0000ii (needs C00ii)
      do ep=-2,0
      S0000ii(z2(1,1),ep)=+2d0*(Bsum00(ep)+Bsum001(ep))
     .                    +2d0*m1*Cv(cc0011+N0,ep)
      S0000ii(z2(1,2),ep)=-2d0*Bsum001(ep)
     .                    +2d0*m1*Cv(cc0012+N0,ep)
      S0000ii(z2(2,2),ep)=+2d0*Bv(bb0011 + B23,ep)
     .                    +2d0*m1*Cv(cc0022+N0,ep)
      enddo

      do i1=1,np
      do i2=i1,np
      call runC_0000ii(kmax,lmax,i1,i2,DetGr,f,Gtwiddle,Gtt, 
     . Shat5zz,Shat6zzzz,S0000ii,Shat6zz,N0) 
      enddo
      enddo
C--- Fixes C00iiii according to extension of Denner-Dittmaier
C--- knowing Ciiii,C00iii,C0000ii, corrections of order Delta*Ciiiiii

c--- calculate S00iiii (needs Ciiii)
      do ep=-2,0
      S00iiii(z4(1,1,1,1),ep)=+2d0*(Bsum0(ep)+3*Bsum1(ep)+3*Bsum11(ep)
     .           +Bsum111(ep))+2d0*m1*Cv(cc1111+N0,ep)
      S00iiii(z4(2,2,2,2),ep)=+2d0*Bv(bb1111 + B23,ep)
     .                        +2d0*m1*Cv(cc2222+N0,ep)
      S00iiii(z4(1,1,1,2),ep)=-2d0*(Bsum1(ep)+2*Bsum11(ep)+Bsum111(ep))
     .                        +2d0*m1*Cv(cc1112+N0,ep)
      S00iiii(z4(1,1,2,2),ep)=2d0*(Bsum11(ep)+Bsum111(ep))
     .                        +2d0*m1*Cv(cc1122+N0,ep)
      S00iiii(z4(1,2,2,2),ep)=-2d0*Bsum111(ep)
     .                        +2d0*m1*Cv(cc1222+N0,ep)
      enddo

      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      call runC_00iiii(kmax,lmax,i1,i2,i3,i4,DetGr,f,Gtwiddle,Gtt,
     . Shat5,Shat6,S00iiii,Shat6zz,N0)
      enddo
      enddo
      enddo
      enddo
C--- Fixes Ciiiii according to extension of Denner-Dittmaier
c--- knowing C00iiii with a correction of order Delta*Ciiiiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      do i5=i4,np
      call runC_iiiii(jmax,i1,i2,i3,i4,i5,DetGr,Xtwiddle0,Gtwiddle,
     . Shat6,N0)
      enddo
      enddo
      enddo
      enddo
      enddo

C----step 4 Calculate Ciiii
 104  continue
C--- Fixes C0000i according to extension of Denner-Dittmaier
C--- knowing C00i,C0000 with corrections of order Delta*C00iii 

c--- calculate S0000i (needs C00i)
      do ep=-2,0
      S0000i(1,ep)=-2d0*Bsum00(ep)+2d0*m1*Cv(cc001+N0,ep)
      S0000i(2,ep)=2d0*Bv(bb001 + B23,ep)+2d0*m1*Cv(cc002+N0,ep)
      enddo
  
      do i1=1,np
      call runC_0000i(kmax,lmax,i1,DetGr,f,Gtwiddle,Gtt, 
     . Shat4zz,Shat5zzzz,S0000i,Shat5zz,N0) 
      enddo
C--- Fixes C00iii according to extension of Denner-Dittmaier
c--- knowing Ciii,C00ii,C0000i with corrections of order Delta*Ciiiii

c--- calculate S00iii (needs Ciii)
      do ep=-2,0
      S00iii(z3(1,1,1),ep)=-2d0*(Bsum0(ep)+2*Bsum1(ep)+Bsum11(ep))
     .                     +2d0*m1*Cv(cc111+N0,ep)
      S00iii(z3(1,1,2),ep)=+2d0*(Bsum1(ep)+Bsum11(ep))
     .                     +2d0*m1*Cv(cc112+N0,ep)
      S00iii(z3(1,2,2),ep)=-2d0*Bsum11(ep)
     .                     +2d0*m1*Cv(cc122+N0,ep)
      S00iii(z3(2,2,2),ep)=+2d0*Bv(bb111 + B23,ep)
     .                     +2d0*m1*Cv(cc222+N0,ep)
      enddo
  
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      call runC_00iii(kmax,lmax,i1,i2,i3,DetGr,f,Gtwiddle,Gtt,
     . Shat4,Shat5,S00iii,Shat5zz,N0)
      enddo
      enddo
      enddo
C--- Fixes Ciiii according to extension of Denner-Dittmaier
C--- knowing C00iii and dropping terms of order Delta*Ciiiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      do i4=i3,np
      call runC_iiii(jmax,i1,i2,i3,i4,DetGr,Xtwiddle0,Gtwiddle,Shat5,N0)
      enddo
      enddo
      enddo
      enddo

C----step 3 
 103  continue
C--- Fixes C0000 using 5.46
C--- knowing C00 with correction of order Delta*C00ii

c--- calculate S0000 (needs C00)
      do ep=-2,0
      S0000(ep)=2d0*Bv(bb00+B23,ep)+2d0*m1*Cv(cc00+N0,ep)
      enddo

      call runC_0000(kmax,lmax,DetGr,f,Gtwiddle,Gtt,Shat3zz,Shat4zz,
     . S0000,N0)
C--- Fixes C00ii using 5.47
C--- knowing Cii,C00i,C0000 with corrections of order Delta*Ciiii

c--- calculate S00ii (needs Cii)
      do ep=-2,0
      S00ii(z2(1,1),ep)=+2d0*(Bsum0(ep)+Bsum1(ep))+2d0*m1*Cv(cc11+N0,ep)
      S00ii(z2(1,2),ep)=-2d0*Bsum1(ep)+2d0*m1*Cv(cc12+N0,ep)
      S00ii(z2(2,2),ep)=+2d0*Bv(bb11 + B23,ep)+2d0*m1*Cv(cc22+N0,ep)
      enddo
  
      do i1=1,np
      do i2=i1,np
      call runC_00ii(kmax,lmax,i1,i2,DetGr,f,Gtwiddle,Gtt,Shat3,Shat4,
     . S00ii,Shat4zz,N0)
      enddo
      enddo
C--- Fixes Ciii using 5.48
c--- knowing C00ii with a correction of order Delta*Ciiii
      do i1=1,np
      do i2=i1,np
      do i3=i2,np
      call runC_iii(jmax,i1,i2,i3,DetGr,Xtwiddle0,Gtwiddle,Shat4,N0)
      enddo
      enddo
      enddo

C---- step 2 calculate C00i and Cii
 102  continue
C--- Fixes C00i according to 5.44 Denner-Dittmaier
C--- knowing Ci, C00 with correction of order Delta*Ciii

c--- calculate S00i (needs Ci)
      do ep=-2,0
      S00i(1,ep)=-2d0*Bsum0(ep)+2d0*m1*Cv(cc1+N0,ep)
      S00i(2,ep)=2d0*Bv(bb1 + B23,ep)+2d0*m1*Cv(cc2+N0,ep)
      enddo
  
      do i3=1,np
      call runC_00i(kmax,lmax,i3,DetGr,f,Gtwiddle,Gtt,Shat2,Shat3,
     . Shat3zz,S00i,N0)
      enddo
C--- Fixes Cii using Eq. 5.45 Denner-Dittmaier
C--- knowing C00i with correction of order Delta*Ciii
      do i2=1,np
      do i3=i2,np
      call runC_ii(jmax,i2,i3,DetGr,Xtwiddle0,Gtwiddle,Shat3,N0)
      enddo
      enddo

C----step 1 calculate C00
 101  continue
C--- Fixes C00 according to 5.42 Denner-Dittmaier
C--- knowing C0 with corrections of order Delta*Cii

c--- calculate S00 (needs C0)
      do ep=-2,0
      S00(ep)=2d0*Bv(bb0+B23,ep)+2d0*m1*Cv(cc0+N0,ep)
      enddo
  
      call runC_00(kmax,lmax,DetGr,f,Gtwiddle,Gtt,Shat1,Shat2,S00,N0)
C--- Fixes Ci according to 5.43 Denner-Dittmaier
C--- knowing C00 with corrections of order Delta*Cii
      do i2=1,np
      call runC_i(jmax,i2,DetGr,Xtwiddle0,Gtwiddle,Shat2,N0)
      enddo

C----step 0 calculate C0
 100  continue
      call runC_0(jmax,DetGr,Xtwiddle0,Gtwiddle,Shat1,N0)

c--- Write out checking information to fort.* files     
c      do ip=1,Ncc
c      write(10+ip,77)'Cv(',ip,'+N0) ',
c     . Cv(ip+N0,-2),Cv(ip+N0,-1),Cv(ip+N0,0)
c      enddo
            
      enddo
    
c--- check the contents of triangle array    
c      write(6,*) 'C array'
c      do ip=1,Ncc
c      if (abs(Csing(ip,p1p2,p1,p2,m1,m2,m3)) .ne. 0d0) then
c      write(6,'(i3,2f20.15)') ip,Cv(ip+N0,-1)
c     .  			 /Csing(ip,p1p2,p1,p2,m1,m2,m3)
c      endif
c      enddo

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

c   77 format(a3,i2,a5,3('(',e13.6,',',e13.6,') '))
    
      end


