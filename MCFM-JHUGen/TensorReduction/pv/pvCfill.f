      subroutine pvCfill(p1,p2,p1p2,m1s,m2s,m3s,N)
      implicit none
C     Calculate the form factors for massless triangle diagrams
C     p1=p1sq,p2=p2sq,p1p2=(p1+p2)^2
C     N is the offset in the common block
      include 'TRconstants.f'
      include 'TRscale.f'
      include 'pvCnames.f'
      include 'pvBnames.f'
      include 'pvBv.f'
      include 'pvCv.f'
      include 'pvverbose.f'
      include 'pvrecurflags.f'
      integer B12,B23,B13,np,ep,epmj,N,j,perm(2),pvBcache
      parameter(np=2)
      double complex G(np,np),in(2,-2:0),trI3,
     . bsum(-2:0),b0sum(-2:0),b1sum(-2:0),b11sum(-2:0),b111sum(-2:0),
     . b1111sum(-2:0),b11111sum(-2:0),
     . b00sum(-2:0),b001sum(-2:0),b0011sum(-2:0),b0000sum(-2:0)
      double precision p1,p2,p1p2,m1s,m2s,m3s,f1,f2
      logical exceptional
      integer,save:: icall,irecur,irecur2,irecur3,irecur4
      double precision,save::idp3(0:2),idp2(0:2),idp1(0:2),id(0:2),
     . idm1(0:2),idm2(0:2)

      logical,save:: first=.true.
      logical,save:: scaleset=.false.
!$omp threadprivate(first,idp3,idp2,idp1,id,idm1,idm2)
!$omp threadprivate(icall,irecur,irecur2,irecur3,irecur4)

      if (first) then
      first=.false.
C--idp3=1/[D+3]
      idp3(0)=1d0/7d0
      idp3(1)=idp3(0)*2d0/7d0
      idp3(2)=idp3(1)*2d0/7d0
C--idp2=1/[D+2]
      idp2(0)=1d0/6d0
      idp2(1)=idp2(0)/3d0
      idp2(2)=idp2(1)/3d0
C--idp1=1/[D+1]
      idp1(0)=0.2d0
      idp1(1)=idp1(0)*0.4d0
      idp1(2)=idp1(1)*0.4d0
C--id=1/D
      id(0)=0.25d0
      id(1)=id(0)*0.5d0
      id(2)=id(1)*0.5d0
C--idm1=1/[D-1]
      idm1(0)=1d0/3d0
      idm1(1)=idm1(0)*2d0/3d0
      idm1(2)=idm1(1)*2d0/3d0
C--idm2=1/[D-2]
      idm2(0)=0.5d0
      idm2(1)=idm2(0)
      idm2(2)=idm2(1)
c--- variables for statistics reporting
      irecur=0
      irecur2=0
      irecur3=0
      irecur4=0
      icall=0      
c--- print out flags for recursion
c      write(6,*) 'pvCfill recursion flags:'
c      write(6,*) '  doGsing  ',doGsing
c      write(6,*) '  doGYsing ',doGYsing
c      write(6,*) '  doPsing  ',doPsing
c      write(6,*) '  doPFsing ',doPFsing
      endif

c--- statistics accounting and reporting      
      icall=icall+1
      if (pvverbose) then
      if (mod(icall,100000) .eq. 0) then
        write(6,77) icall,
     & 1d2*dfloat(icall-irecur-irecur2-irecur3-irecur4)/dfloat(icall),
     & 1d2*dfloat(irecur)/dfloat(icall),
     & 1d2*dfloat(irecur2)/dfloat(icall),
     & 1d2*dfloat(irecur3)/dfloat(icall),
     & 1d2*dfloat(irecur4)/dfloat(icall)
      endif
      endif
   77 format(' Cfill ',i9,': ',5(f6.2,'% : '))
      
      B12=pvBcache(p1,m1s,m2s)
      B23=pvBcache(p2,m2s,m3s)
      B13=pvBcache(p1p2,m1s,m3s)

      f1=m2s-m1s-p1
      f2=m3s-m1s-p1p2

      G(1,1)=dcmplx(2d0*p1)
      G(2,2)=dcmplx(2d0*p1p2)
      G(1,2)=dcmplx(p1+p1p2-p2)
      G(2,1)=G(1,2)

c      if (pvverbose) write(6,*) 'Check triangle Gsing'
c      Gsing=pvGramsing(G,2)

C     Y(i,j)=mi^2+mj^2-(q_i-q_j)^2
C     where q_1=0,  q_2=p1,  q_3=p_1+p_2;

c      Y(1,1) = dcmplx(2d0*m1s)
c      Y(1,2) = dcmplx(m1s + m2s - p1)
c      Y(2,1) = Y(1,2)
c      Y(1,3) = dcmplx(m1s + m3s - p1p2)
c      Y(3,1) = Y(1,3)
c      Y(2,2) = dcmplx(2d0*m2s)
c      Y(2,3) = dcmplx(m2s + m3s - p2)
c      Y(3,2) = Y(2,3)
c      Y(3,3) = dcmplx(2d0*m3s)
      
c      if (pvverbose) write(6,*) 'Check triangle Ysing'
c      Ysing=pvGramsing(Y,3)

c-- find maximum entry in Gram matrix
c      Gmax=zip
c      do j=1,2
c      do k=j,2
c      if (abs(G(j,k)) .gt. Gmax) Gmax=abs(G(j,k))
c      enddo
c      enddo
      
c      if (pvverbose) write(6,*) 'Gmax=',Gmax
c      Psing=.false.
c--- criterion for small momenta recursion
c      if (Gmax .lt. weenumber) Psing=.true.
      
c-- find maximum of f1 and f2
c      fmax=max(abs(f1),abs(f2))
      
c      if (pvverbose) write(6,*) 'fmax=',fmax
c      Fsing=.false.
c--- criterion for small momenta and small f(k) recursion
c      if (fmax .lt. weenumber) Fsing=.true.
      
c      write(6,*) 'Gsing,Ysing,Psing,Fsing',Gsing,Ysing,Psing,Fsing
      
      exceptional=.false.
      
      if     (doPFsing) then
c--- for small momenta and small f(k)
        if (pvverbose) then
          write(6,*) 'USING TRIANGLE SMALL MOMENTA AND f(k) RECURSION'
      endif
        call Cfill_recur4(p1,p2,p1p2,m1s,m2s,m3s,N)
      irecur4=irecur4+1
c        write(98,'(4(l2),11(f21.15))') Gsing,Ysing,Psing,Fsing,
c     &   q1save(1),q1save(2),q1save(3),q1save(4),
c     &   q2save(1),q2save(2),q2save(3),q2save(4),
c     &   m1s,m2s,m3s
c        write(98,'(6(f21.15))') p1,p2,p1p2,m1s,m2s,m3s
      return
      elseif (doPsing) then
        if (pvverbose) then
        write(6,*) 'USING TRIANGLE SMALL MOMENTA RECURSION'
      endif
        call Cfill_recur3(p1,p2,p1p2,m1s,m2s,m3s,N)
      irecur3=irecur3+1
c        write(98,'(4(l2),11(f21.15))') Gsing,Ysing,Psing,Fsing,
c     &   q1save(1),q1save(2),q1save(3),q1save(4),
c     &   q2save(1),q2save(2),q2save(3),q2save(4),
c     &   m1s,m2s,m3s
c        write(98,'(6(f21.15))') p1,p2,p1p2,m1s,m2s,m3s
      return
      elseif (doGYsing) then
c--- for small Gram and small Y
        if (pvverbose) then
          write(6,*) 'USING TRIANGLE SMALL Y AND SMALL G RECURSION'
      endif
        call Cfill_recur2(p1,p2,p1p2,m1s,m2s,m3s,N,exceptional)
      irecur2=irecur2+1
c        write(98,'(4(l2),11(f21.15))') Gsing,Ysing,Psing,Fsing,
c     &   q1save(1),q1save(2),q1save(3),q1save(4),
c     &   q2save(1),q2save(2),q2save(3),q2save(4),
c     &   m1s,m2s,m3s
c        write(98,'(6(f21.15))') p1,p2,p1p2,m1s,m2s,m3s
        if (exceptional) then
c------ for exceptional configurations, fall through to normal PV
        continue
      else
c------ otherwise, we're done
        return
      endif
      elseif (doGsing) then
c--- for small Gram only 
        if (pvverbose) then  
          write(6,*) 'USING TRIANGLE SMALL G RECURSION'
      endif
        call Cfill_recur (p1,p2,p1p2,m1s,m2s,m3s,N)
      irecur=irecur+1
c        write(98,'(4(l2),11(f21.15))') Gsing,Ysing,Psing,Fsing,
c     &   q1save(1),q1save(2),q1save(3),q1save(4),
c     &   q2save(1),q2save(2),q2save(3),q2save(4),
c     &   m1s,m2s,m3s
c        write(98,'(6(f21.15))') p1,p2,p1p2,m1s,m2s,m3s
      return
      endif
c--- otherwise, usual PV is fine

c      if (exceptional) write(6,*) 'WARNING: EXCEPTIONAL POINT'

c        write(98,'(4(l2),11(f21.15))') Gsing,Ysing,Psing,Fsing,
c     &   q1save(1),q1save(2),q1save(3),q1save(4),
c     &   q2save(1),q2save(2),q2save(3),q2save(4),
c     &   m1s,m2s,m3s
c        write(98,'(6(f21.15))') p1,p2,p1p2,m1s,m2s,m3s

c--- initialize integrals      
      do ep=-2,0
      do j=1,Ncc
      Cv(N+j,ep)=dcmplx(1d5,-1d5)
      enddo
      enddo

      call XLUDecomp(G, 2, perm)

C---one index form factors
c'B'Id,Bsum111(P?,K?,m1?,m2?)=MM(B111,P,K,m1,m2)+MM(B1111,P,K,m1,m2);
c'B'Id,Bsum001(P?,K?,m1?,m2?)=MM(B001,P,K,m1,m2)+MM(B0011,P,K,m1,m2);
c'B'Id,Bsum00(P?,K?,m1?,m2?)=MM(B00,P,K,m1,m2)+MM(B001,P,K,m1,m2);
c'B'Id,Bsum11(P?,K?,m1?,m2?)=MM(B11,P,K,m1,m2)+MM(B111,P,K,m1,m2);
c'B'Id,Bsum1(P?,K?,m1?,m2?)=MM(B1,P,K,m1,m2)+MM(B11,P,K,m1,m2);
c'B'Id,Bsum0(P?,K?,m1?,m2?)=MM(B0,P,K,m1,m2)+MM(B1,P,K,m1,m2);

      do ep=-2,0
      Cv(N+cc0,ep)=trI3(p1,p2,p1p2,m1s,m2s,m3s,musq,ep)
      bsum(ep)=Bv(bb0+B23,ep)+Bv(bb1+B23,ep)

      b0sum(ep)=Bv(bb0+B23,ep)+Bv(bb1+B23,ep)
      b00sum(ep)=Bv(bb00+B23,ep)+Bv(bb001+B23,ep)
      b001sum(ep)=Bv(bb001+B23,ep)+Bv(bb0011+B23,ep)
      b0011sum(ep)=Bv(bb0011+B23,ep)+Bv(bb00111+B23,ep)

      b0000sum(ep)=Bv(bb0000+B23,ep)+Bv(bb00001+B23,ep)

      b1sum(ep)=Bv(bb1+B23,ep)+Bv(bb11+B23,ep)
      b11sum(ep)=Bv(bb11+B23,ep)+Bv(bb111+B23,ep)
      b111sum(ep)=Bv(bb111+B23,ep)+Bv(bb1111+B23,ep)
      b1111sum(ep)=Bv(bb1111+B23,ep)+Bv(bb11111+B23,ep)
      b11111sum(ep)=Bv(bb11111+B23,ep)+Bv(bb111111+B23,ep)

      in(1,ep)=f1*Cv(N+cc0,ep)-Bv(bb0+B23,ep)+Bv(bb0+B13,ep)
      in(2,ep)=f2*Cv(N+cc0,ep)-Bv(bb0+B23,ep)+Bv(bb0+B12,ep)
      enddo
      call pvBackSubst(G,2,perm,in)

      do ep=-2,0
      Cv(N+cc1,ep)=in(1,ep)
      Cv(N+cc2,ep)=in(2,ep)
c      write(66,*) 'ox',cc1,ep,Cv(N+cc1,ep)
c      write(66,*) 'ox',cc2,ep,Cv(N+cc2,ep)
      enddo
      
C---two index form factors
      do ep=-2,0
      Cv(N+cc00,ep)=czip
      if (ep .eq. -2) goto 20
      do j=0,ep+2
      epmj=ep-j
      Cv(N+cc00,ep)=Cv(N+cc00,ep)+idm2(j)*(m1s*Cv(N+cc0,epmj)
     & -half*(f1*Cv(N+cc1,epmj)+f2*Cv(N+cc2,epmj)-Bv(bb0+B23,epmj)))
      enddo
 20   continue

      in(1,ep)=f1*Cv(N+cc1,ep)+bsum(ep)-2d0*Cv(N+cc00,ep)
      in(2,ep)=f2*Cv(N+cc1,ep)+bsum(ep)+Bv(bb1+B12,ep)
      enddo
      call pvBackSubst(G,2,perm,in)
      do ep=-2,0
      Cv(N+cc11,ep)=in(1,ep)
      Cv(N+cc12,ep)=in(2,ep)
c      write(66,*) 'ox',cc11,ep,Cv(N+cc11,ep)
c      write(66,*) 'ox',cc12,ep,Cv(N+cc12,ep)
      
      in(1,ep)=f1*Cv(N+cc2,ep)-Bv(bb1+B23,ep)+Bv(bb1+B13,ep)
      in(2,ep)=f2*Cv(N+cc2,ep)-Bv(bb1+B23,ep)-2d0*Cv(N+cc00,ep)
      enddo
      call pvBackSubst(G,2,perm,in)
      do ep=-2,0
      Cv(N+cc12,ep)=0.5D0*(Cv(N+cc12,ep)+in(1,ep))
      Cv(N+cc22,ep)=in(2,ep)
c      write(66,*) 'ox',cc12,ep,Cv(N+cc12,ep)
c      write(66,*) 'ox',cc22,ep,Cv(N+cc22,ep)
      enddo
       
c      if ((maxcindex .eq. 2) .and. (pvRespectmaxcindex)) return

C---three index form factors
      do ep=-2,0
      Cv(N+cc001,ep)=czip
      Cv(N+cc002,ep)=czip
      if (ep .eq. -2) goto 30
      do j=0,ep+2
      epmj=ep-j
      Cv(N+cc001,ep)=Cv(N+cc001,ep)+idm1(j)*(m1s*Cv(N+cc1,epmj)
     & -half*(f1*Cv(N+cc11,epmj)+f2*Cv(N+cc12,epmj)+bsum(epmj)))
      Cv(N+cc002,ep)=Cv(N+cc002,ep)+idm1(j)*(m1s*Cv(N+cc2,epmj)
     & -half*(f1*Cv(N+cc12,epmj)+f2*Cv(N+cc22,epmj)-Bv(bb1+B23,epmj)))
      enddo
 30   continue
      enddo

      do ep=-2,0
      bsum(ep)=bsum(ep)+b1sum(ep)
C--- bsum is now equal to 
c--- Bv(bb1+B23,ep)+2*Bv(bb1+B23,ep)+Bv(bb11+B23,ep)
      in(1,ep)=f1*Cv(N+cc11,ep)-bsum(ep)-4d0*Cv(N+cc001,ep)
      in(2,ep)=f2*Cv(N+cc11,ep)-bsum(ep)+Bv(bb11+B12,ep)
      enddo
      call pvBackSubst(G,2,perm,in)
      do ep=-2,0
      Cv(N+cc111,ep)=in(1,ep)
      Cv(N+cc112,ep)=in(2,ep)
c      write(66,*) 'ox',cc111,ep,Cv(N+cc111,ep)
c      write(66,*) 'ox',cc112,ep,Cv(N+cc112,ep)

      in(1,ep)=f1*Cv(N+cc22,ep)-Bv(bb11+B23,ep)+Bv(bb11+B13,ep)
      in(2,ep)=f2*Cv(N+cc22,ep)-Bv(bb11+B23,ep)-4d0*Cv(N+cc002,ep)
      enddo
      call pvBackSubst(G,2,perm,in)
      do ep=-2,0
      Cv(N+cc122,ep)=in(1,ep)
      Cv(N+cc222,ep)=in(2,ep)
c      write(66,*) 'ox',cc122,ep,Cv(N+cc122,ep)
c      write(66,*) 'ox',cc222,ep,Cv(N+cc222,ep)

c      b1sum(ep)=Bv(bb1+B23,ep)+Bv(bb11+B23,ep)
      in(1,ep)=f1*Cv(N+cc12,ep)+b1sum(ep)-2d0*Cv(N+cc002,ep)
      in(2,ep)=f2*Cv(N+cc12,ep)+b1sum(ep)-2d0*Cv(N+cc001,ep)
      enddo
      call pvBackSubst(G,2,perm,in)
      do ep=-2,0
      Cv(N+cc112,ep)=in(1,ep)
      Cv(N+cc122,ep)=in(2,ep)
c      write(66,*) 'ox',cc112,ep,Cv(N+cc112,ep)
c      write(66,*) 'ox',cc122,ep,Cv(N+cc122,ep)
      enddo

c      if ((maxcindex .eq. 3) .and. (pvRespectmaxcindex)) return

C---four index form factors
      do ep=-2,0
      do j=cc0000,cc0022
      Cv(N+j,ep)=czip
      enddo
      do j=0,ep+2
      epmj=ep-j
      Cv(N+cc0000,ep)=Cv(N+cc0000,ep)+id(j)*(
     & +m1s*Cv(N+cc00,epmj)
     & -half*(f1*Cv(N+cc001,epmj)+f2*Cv(N+cc002,epmj)
     & -Bv(bb00+B23,epmj)))
      Cv(N+cc0011,ep)=Cv(N+cc0011,ep)+id(j)*(
     & +m1s*Cv(N+cc11,epmj)
     & -half*(f1*Cv(N+cc111,epmj)+f2*Cv(N+cc112,epmj)
     & -b0sum(epmj)-b1sum(epmj)))
      Cv(N+cc0012,ep)=Cv(N+cc0012,ep)+id(j)*(
     & +m1s*Cv(N+cc12,epmj)
     & -half*(f1*Cv(N+cc112,epmj)+f2*Cv(N+cc122,epmj)
     & +b1sum(epmj)))
      Cv(N+cc0022,ep)=Cv(N+cc0022,ep)+id(j)*(
     & +m1s*Cv(N+cc22,epmj)
     & -half*(f1*Cv(N+cc122,epmj)+f2*Cv(N+cc222,epmj)
     & -Bv(bb11+B23,epmj)))

      enddo
c      write(66,*) 'ox',cc0000,ep,Cv(N+cc0000,ep)
c      write(66,*) 'ox',cc0011,ep,Cv(N+cc0011,ep)
c      write(66,*) 'ox',cc0012,ep,Cv(N+cc0012,ep)
c      write(66,*) 'ox',cc0022,ep,Cv(N+cc0022,ep)
      enddo

      do ep=-2,0
      bsum(ep)=bsum(ep)+b1sum(ep)+b11sum(ep)
C--- bsum is now equal to 
c--- Bv(bb1+B23,ep)+3*Bv(bb1+B23,ep)+3*Bv(bb11+B23,ep)+Bv(bb111+B23,ep)

      in(1,ep)=f1*Cv(N+cc111,ep)+bsum(ep)-6d0*Cv(N+cc0011,ep)
      in(2,ep)=f2*Cv(N+cc111,ep)+bsum(ep)+Bv(bb111+B12,ep)
      enddo
      call pvBackSubst(G,2,perm,in)
      do ep=-2,0
      Cv(N+cc1111,ep)=in(1,ep)
      Cv(N+cc1112,ep)=in(2,ep)
c      write(66,*) 'ox',cc1111,ep,Cv(N+cc1111,ep)
c      write(66,*) 'ox',cc1112,ep,Cv(N+cc1112,ep)

      in(1,ep)=f1*Cv(N+cc112,ep)-b1sum(ep)-b11sum(ep)
     . -4d0*Cv(N+cc0012,ep)
      in(2,ep)=f2*Cv(N+cc112,ep)-b1sum(ep)-b11sum(ep)
     . -2d0*Cv(N+cc0011,ep)
      enddo
      call pvBackSubst(G,2,perm,in)
      do ep=-2,0
      Cv(N+cc1112,ep)=in(1,ep)
      Cv(N+cc1122,ep)=in(2,ep)
c      write(66,*) 'ox',cc1112,ep,Cv(N+cc1112,ep)
c      write(66,*) 'ox',cc1122,ep,Cv(N+cc1122,ep)

      in(1,ep)=f1*Cv(N+cc222,ep)-Bv(bb111+B23,ep)+Bv(bb111+B13,ep)
      in(2,ep)=f2*Cv(N+cc222,ep)-Bv(bb111+B23,ep)-6d0*Cv(N+cc0022,ep)
      enddo
      call pvBackSubst(G,2,perm,in)
      do ep=-2,0
      Cv(N+cc1222,ep)=in(1,ep)
      Cv(N+cc2222,ep)=in(2,ep)
c      write(66,*) 'ox',cc1222,ep,Cv(N+cc1222,ep)
c      write(66,*) 'ox',cc2222,ep,Cv(N+cc2222,ep)

      in(1,ep)=f1*Cv(N+cc122,ep)+b11sum(ep)-2d0*Cv(N+cc0022,ep)
      in(2,ep)=f2*Cv(N+cc122,ep)+b11sum(ep)-4d0*Cv(N+cc0012,ep)
      enddo
      call pvBackSubst(G,2,perm,in)
      do ep=-2,0
      Cv(N+cc1122,ep)=in(1,ep)
      Cv(N+cc1222,ep)=half*(Cv(N+cc1222,ep)+in(2,ep))
      enddo 

c      if ((maxcindex .eq. 4) .and. (pvRespectmaxcindex)) return

C---five index form factors
      do ep=-2,0
      do j=cc00001,cc00222
      Cv(N+j,ep)=czip
      enddo
      do j=0,ep+2
      epmj=ep-j
      Cv(N+cc00001,ep)=Cv(N+cc00001,ep)+idp1(j)*(
     & +m1s*Cv(N+cc001,epmj)
     & -half*(f1*Cv(N+cc0011,epmj)+f2*Cv(N+cc0012,epmj)
     & +b00sum(epmj)))
      Cv(N+cc00002,ep)=Cv(N+cc00002,ep)+idp1(j)*(
     & +m1s*Cv(N+cc002,epmj)
     & -half*(f1*Cv(N+cc0012,epmj)+f2*Cv(N+cc0022,epmj)
     & -Bv(bb001+B23,epmj)))
      Cv(N+cc00111,ep)=Cv(N+cc00111,ep)+idp1(j)*(
     & +m1s*Cv(N+cc111,epmj)
     & -half*(f1*Cv(N+cc1111,epmj)+f2*Cv(N+cc1112,epmj)
     & +b0sum(epmj)+2D0*b1sum(epmj)+b11sum(epmj)))
      Cv(N+cc00112,ep)=Cv(N+cc00112,ep)+idp1(j)*(
     & +m1s*Cv(N+cc112,epmj)
     & -half*(f1*Cv(N+cc1112,epmj)+f2*Cv(N+cc1122,epmj)
     & -b1sum(epmj)-b11sum(epmj)))
      Cv(N+cc00122,ep)=Cv(N+cc00122,ep)+idp1(j)*(
     & +m1s*Cv(N+cc122,epmj)
     & -half*(f1*Cv(N+cc1122,epmj)+f2*Cv(N+cc1222,epmj)
     & +b11sum(epmj)))
      Cv(N+cc00222,ep)=Cv(N+cc00222,ep)+idp1(j)*(
     & +m1s*Cv(N+cc222,epmj)
     & -half*(f1*Cv(N+cc1222,epmj)+f2*Cv(N+cc2222,epmj)
     & -Bv(bb111+B23,epmj)))
      enddo

c      write(66,*) 'ox',cc00001,ep,Cv(N+cc00001,ep)
c      write(66,*) 'ox',cc00002,ep,Cv(N+cc00002,ep)
c      write(66,*) 'ox',cc00111,ep,Cv(N+cc00111,ep)
c      write(66,*) 'ox',cc00112,ep,Cv(N+cc00112,ep)
c      write(66,*) 'ox',cc00122,ep,Cv(N+cc00122,ep)
c      write(66,*) 'ox',cc00222,ep,Cv(N+cc00222,ep)
      enddo

      do ep=-2,0
C Cv(pppp
      in(1,ep) = f1*Cv(N+cc1111,ep)-8.D0*Cv(N+cc00111,ep)- Bv(bb0+B23,
     & ep)-4.D0*Bv(bb1+B23,ep)-6.D0*Bv(bb11+B23,ep)-4.D0*Bv(
     & bb111+B23,ep)-Bv(bb1111+B23,ep)
      in(2,ep) = f2*Cv(N+cc1111,ep)-Bv(bb0+B23,ep)-4.D0*Bv(bb1+
     & B23,ep)-6.D0*Bv(bb11+B23,ep)-4.D0*Bv(bb111+B23,ep)+Bv(
     & bb1111+B12,ep)-Bv(bb1111+B23,ep)
      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc11111,ep)=in(1,ep)
      Cv(N+cc11112,ep)=in(2,ep)
c      write(66,*) 'ox',cc11111,ep,Cv(N+cc11111,ep)
c      write(66,*) 'ox',cc11112,ep,Cv(N+cc11112,ep)
     
C Cv(pppk
      in(1,ep) = f1*Cv(N+cc1112,ep)-6.D0*Cv(N+cc00112,ep)+ Bv(bb1+B23,
     & ep)+3.D0*Bv(bb11+B23,ep)+3.D0*Bv(bb111+B23,ep)+Bv(
     & bb1111+B23,ep)
      in(2,ep) = f2*Cv(N+cc1112,ep)-2.D0*Cv(N+cc00111,ep) +Bv(bb1+B23,
     & ep)+3.D0*Bv(bb11+B23,ep)+3.D0*Bv(bb111+B23,ep)+Bv(
     & bb1111+B23,ep)
      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc11112,ep)=in(1,ep)
      Cv(N+cc11122,ep)=in(2,ep)
c      write(66,*) 'ox',cc11112,ep,Cv(N+cc11112,ep)
c      write(66,*) 'ox',cc11122,ep,Cv(N+cc11122,ep)
     
C Cv(pkkk
      in(1,ep) = f1*Cv(N+cc1222,ep)-2.D0*Cv(N+cc00222,ep)+Bv(bb111+
     & B23,ep)+Bv(bb1111+B23,ep)

      in(2,ep) = f2*Cv(N+cc1222,ep)-6.D0*Cv(N+cc00122,ep)+Bv(bb111+
     & B23,ep)+Bv(bb1111+B23,ep)

      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc11222,ep)=in(1,ep)
      Cv(N+cc12222,ep)=in(2,ep)
c      write(66,*) 'ox',cc11222,ep,Cv(N+cc11222,ep)
c      write(66,*) 'ox',cc12222,ep,Cv(N+cc12222,ep)
     
C Cv(kkkk
      in(1,ep) = f1*Cv(N+cc2222,ep)+Bv(bb1111+B13,ep)-Bv(bb1111+B23,ep)
      in(2,ep) = f2*Cv(N+cc2222,ep)-8.D0*Cv(N+cc00222,ep)
     & -Bv(bb1111+B23,ep)
      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc12222,ep)=in(1,ep)
      Cv(N+cc22222,ep)=in(2,ep)
c      write(66,*) 'ox',cc12222,ep,Cv(N+cc12222,ep)
c      write(66,*) 'ox',cc22222,ep,Cv(N+cc22222,ep)
      enddo

c      if ((maxcindex .eq. 5) .and. (pvRespectmaxcindex)) return

C---six index form factors
      do ep=-2,0
      do j=cc000000,cc002222
      Cv(N+j,ep)=czip
      enddo
      if (ep .eq. -2) goto 60
      do j=0,ep+2
      epmj=ep-j
      Cv(N+cc000000,ep)=Cv(N+cc000000,ep)+idp2(j)*(
     &  +m1s*Cv(N+cc0000,epmj)
     &  -half*(f1*Cv(N+cc00001,epmj)+f2*Cv(N+cc00002,epmj)
     &  -Bv(bb0000+B23,epmj)))

      Cv(N+cc000011,ep)=Cv(N+cc000011,ep)+idp2(j)*(
     &  +m1s*Cv(N+cc0011,epmj)
     &  -half*(f1*Cv(N+cc00111,epmj)+f2*Cv(N+cc00112,epmj)
     &  -B00sum(epmj)-B001sum(epmj)))

      Cv(N+cc000012,ep)=Cv(N+cc000012,ep)+idp2(j)*(
     &  +m1s*Cv(N+cc0012,epmj)
     &  -half*(f1*Cv(N+cc00112,epmj)+f2*Cv(N+cc00122,epmj)
     &  +B001sum(epmj)))

      Cv(N+cc000022,ep)=Cv(N+cc000022,ep)+idp2(j)*(
     &  +m1s*Cv(N+cc0022,epmj)
     &  -half*(f1*Cv(N+cc00122,epmj)+ f2*Cv(N+cc00222,epmj)
     &  -Bv(bb0011+B23,epmj)))

      Cv(N+cc001111,ep)=Cv(N+cc001111,ep)+idp2(j)*(
     &  +m1s*Cv(N+cc1111,epmj)
     &  -half*(f1*Cv(N+cc11111,epmj)+f2*Cv(N+cc11112,epmj)
     &  -B111sum(epmj)-B0sum(epmj)-3d0*B1sum(epmj)-3d0*B11sum(epmj)))

      Cv(N+cc001112,ep)=Cv(N+cc001112,ep)+idp2(j)*(
     &  +m1s*Cv(N+cc1112,epmj)
     &  -half*(f1*Cv(N+cc11112,epmj)+f2*Cv(N+cc11122,epmj)
     &  +B1sum(epmj)+2d0*B11sum(epmj)+B111sum(epmj)))
      Cv(N+cc001122,ep)=Cv(N+cc001122,ep)+idp2(j)*(
     &  +m1s*Cv(N+cc1122,epmj)
     &  -half*(f1*Cv(N+cc11122,epmj)+f2*Cv(N+cc11222,epmj)
     &  -B11sum(epmj)-B111sum(epmj))) 
      Cv(N+cc001222,ep)=Cv(N+cc001222,ep)+idp2(j)*(
     &  +m1s*Cv(N+cc1222,epmj)
     &  -half*(f1*Cv(N+cc11222,epmj)+f2*Cv(N+cc12222,epmj)
     &  +B111sum(epmj)))

      Cv(N+cc002222,ep)=Cv(N+cc002222,ep)+idp2(j)*(
     &  +m1s*Cv(N+cc2222,epmj)
     &  -half*(f1*Cv(N+cc12222,epmj)+f2*Cv(N+cc22222,epmj)
     &  -Bv(bb1111+B23,epmj)))

      enddo
 60   continue

c      write(66,*) 'ox',cc000000,ep,Cv(N+cc000000,ep)
c      write(66,*) 'ox',cc000011,ep,Cv(N+cc000011,ep)
c      write(66,*) 'ox',cc000012,ep,Cv(N+cc000012,ep)
c      write(66,*) 'ox',cc000022,ep,Cv(N+cc000022,ep)

c      write(66,*) 'ox',cc001111,ep,Cv(N+cc001111,ep)
c      write(66,*) 'ox',cc001112,ep,Cv(N+cc001112,ep)
c      write(66,*) 'ox',cc001122,ep,Cv(N+cc001122,ep)
c      write(66,*) 'ox',cc001222,ep,Cv(N+cc001222,ep)
c      write(66,*) 'ox',cc002222,ep,Cv(N+cc002222,ep)

      enddo
      
C    (Cv(N+cc000011,Cv(N+cc000012,zzzzp)
      do ep=-2,0
      in(1,ep)=
     .   + f1*Cv(N+cc00001,ep)
     .   -2*Cv(N+cc000000,ep)
     .   + B0000sum(ep)
      in(2,ep)=
     .   + f2*Cv(N+cc00001,ep)
     .   + Bv(bb00001+B12,ep)
     .   + B0000sum(ep)

      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc000011,ep)=in(1,ep)
      Cv(N+cc000012,ep)=in(2,ep)
c      write(66,*) 'ox',cc000011,ep,Cv(N+cc000011,ep)
c      write(66,*) 'ox',cc000012,ep,Cv(N+cc000012,ep)
     
C    (Cv(N+cc000012,Cv(N+cc000022,zzzzk)
      in(1,ep)=
     .   + f1*Cv(N+cc00002,ep)
     .   + Bv(bb00001+B13,ep)
     .   - Bv(bb00001+B23,ep)
      in(2,ep)=
     .   + f2*Cv(N+cc00002,ep)
     .   -2*Cv(N+cc000000,ep)
     .   - Bv(bb00001+B23,ep)
      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc000012,ep)=in(1,ep)
      Cv(N+cc000022,ep)=in(2,ep)

c      write(66,*) 'ox',cc000012,ep,Cv(N+cc000012,ep)
c      write(66,*) 'ox',cc000022,ep,Cv(N+cc000022,ep)

C    (Cv(N+cc001111,Cv(N+cc001112,zzppp)
      in(1,ep)=
     .   + f1*Cv(N+cc00111,ep)-6*Cv(N+cc000011,ep)
     .   + B0011sum(ep)+ 2d0*B001sum(ep)+ B00sum(ep)
      in(2,ep)=
     .  + f2*Cv(N+cc00111,ep)+ Bv(bb00111+B12,ep)
     .   + B0011sum(ep)+ 2d0*B001sum(ep)+ B00sum(ep)
      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc001111,ep)=in(1,ep)
      Cv(N+cc001112,ep)=in(2,ep)
c      write(66,*) 'ox',cc001111,ep,Cv(N+cc001111,ep)
c      write(66,*) 'ox',cc001112,ep,Cv(N+cc001112,ep)
     
C    (Cv(N+cc001112,Cv(N+cc001122,zzppk)
      in(1,ep)=
     .   + f1*Cv(N+cc00112,ep)-4*Cv(N+cc000012,ep)
     .    -B001sum(ep)-B0011sum(ep)
      in(2,ep)=
     .   + f2*Cv(N+cc00112,ep)
     .   -2*Cv(N+cc000011,ep)
     .    -B001sum(ep)-B0011sum(ep)

      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc001112,ep)=in(1,ep)
      Cv(N+cc001122,ep)=in(2,ep)
c      write(66,*) 'ox',cc001112,ep,Cv(N+cc001112,ep)
c      write(66,*) 'ox',cc001122,ep,Cv(N+cc001122,ep)
     
C    (Cv(N+cc001122,Cv(N+cc001222,zzpkk)
      in(1,ep)=
     .   + f1*Cv(N+cc00122,ep)-2*Cv(N+cc000022,ep)
     .   + B0011sum(ep)

      in(2,ep)=
     .   + f2*Cv(N+cc00122,ep)-4*Cv(N+cc000012,ep)
     .   + B0011sum(ep)

      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc001122,ep)=in(1,ep)
      Cv(N+cc001222,ep)=in(2,ep)
c      write(66,*) 'ox',cc001122,ep,Cv(N+cc001122,ep)
c      write(66,*) 'ox',cc001222,ep,Cv(N+cc001222,ep)
     
C    (Cv(N+cc001222,Cv(N+cc002222,zzkkk)
      in(1,ep)=
     .   + f1*Cv(N+cc00222,ep)
     .   + Bv(bb00111+B13,ep)- Bv(bb00111+B23,ep)
      in(2,ep)=
     .   + f2*Cv(N+cc00222,ep)-6*Cv(N+cc000022,ep)
     .   - Bv(bb00111+B23,ep)


      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc001222,ep)=in(1,ep)
      Cv(N+cc002222,ep)=in(2,ep)
c      write(66,*) 'ox',cc001222,ep,Cv(N+cc001222,ep)
c      write(66,*) 'ox',cc002222,ep,Cv(N+cc002222,ep)
     


C    (Cv(N+cc111111,Cv(N+cc111112,ppppp)
      in(1,ep)=
     .   +f1*Cv(N+cc11111,ep)-10*Cv(N+cc001111,ep)
     .   +B0sum(ep)+4*B1sum(ep)+4*B111sum(ep)
     .   +6*B11sum(ep)+B1111sum(ep)
      in(2,ep)=
     .   +f2*Cv(N+cc11111,ep)+Bv(bb11111+B12,ep)
     .   +B0sum(ep)+4*B1sum(ep)
     .   +4*B111sum(ep)+6*B11sum(ep)+B1111sum(ep)
      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc111111,ep)=in(1,ep)
      Cv(N+cc111112,ep)=in(2,ep)
c      write(66,*) 'ox',cc111111,ep,Cv(N+cc111111,ep)
c      write(66,*) 'ox',cc111112,ep,Cv(N+cc111112,ep)
     
C    (Cv(N+cc111112,Cv(N+cc111122,ppppk)
      in(1,ep)=
     .    +f1*Cv(N+cc11112,ep)-8*Cv(N+cc001112,ep)
     .   -B1sum(ep)-3*B11sum(ep)-3*B111sum(ep)-B1111sum(ep)

      in(2,ep)=
     .   + f2*Cv(N+cc11112,ep)-2*Cv(N+cc001111,ep)
     .   -B1sum(ep)-3*B11sum(ep)-3*B111sum(ep)-B1111sum(ep)

      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc111112,ep)=in(1,ep)
      Cv(N+cc111122,ep)=in(2,ep)
c      write(66,*) 'ox',cc111112,ep,Cv(N+cc111112,ep)
c      write(66,*) 'ox',cc111122,ep,Cv(N+cc111122,ep)
     
C    (Cv(N+cc111122,Cv(N+cc111222,pppkk)
      in(1,ep)=
     .   + f1*Cv(N+cc11122,ep)-6*Cv(N+cc001122,ep)
     .   +B11sum(ep)+2*B111sum(ep)+B1111sum(ep)
      in(2,ep)=
     .   + f2*Cv(N+cc11122,ep)-4*Cv(N+cc001112,ep)
     .   +B11sum(ep)+2*B111sum(ep)+B1111sum(ep)

      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc111122,ep)=in(1,ep)
      Cv(N+cc111222,ep)=in(2,ep)
c      write(66,*) 'ox',cc111122,ep,Cv(N+cc111122,ep)
c      write(66,*) 'ox',cc111222,ep,Cv(N+cc111222,ep)
     
C    (Cv(N+cc111222,Cv(N+cc112222,ppkkk)
      in(1,ep)=
     .   + f1*Cv(N+cc11222,ep)-4*Cv(N+cc001222,ep)
     .   -B111sum(ep)-B1111sum(ep)
      in(2,ep)=
     .   + f2*Cv(N+cc11222,ep)-6*Cv(N+cc001122,ep)
     .   -B111sum(ep)-B1111sum(ep)

      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc111222,ep)=in(1,ep)
      Cv(N+cc112222,ep)=in(2,ep)
c      write(66,*) 'ox',cc111222,ep,Cv(N+cc111222,ep)
c      write(66,*) 'ox',cc112222,ep,Cv(N+cc112222,ep)

C    (Cv(N+cc112222,Cv(N+cc122222,pkkkk)
      in(1,ep)=f1*Cv(N+cc12222,ep)-2*Cv(N+cc002222,ep)
     .   + B1111sum(ep)

      in(2,ep)=f2*Cv(N+cc12222,ep)-8*Cv(N+cc001222,ep)
     .   + B1111sum(ep)

      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc112222,ep)=in(1,ep)
      Cv(N+cc122222,ep)=in(2,ep)
c      write(66,*) 'ox',cc112222,ep,Cv(N+cc112222,ep)
c      write(66,*) 'ox',cc122222,ep,Cv(N+cc122222,ep)
     
C    (Cv(N+cc122222,Cv(N+cc222222,kkkkk)
      in(1,ep)=
     .   + f1*Cv(N+cc22222,ep)
     .   + Bv(bb11111+B13,ep)-Bv(bb11111+B23,ep)
 
      in(2,ep)=
     .   + f2*Cv(N+cc22222,ep)-10* Cv(N+cc002222,ep)
     .   - Bv(bb11111+B23,ep)
      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc122222,ep)=in(1,ep)
      Cv(N+cc222222,ep)=in(2,ep)
c      write(66,*) 'ox',cc122222,ep,Cv(N+cc122222,ep)
c      write(66,*) 'ox',cc222222,ep,Cv(N+cc222222,ep)

      enddo    

c      if ((maxcindex .eq. 6) .and. (pvRespectmaxcindex)) return

C---seven index form factors
      do ep=-2,0
      do j=cc0000001,cc0022222
      Cv(N+j,ep)=czip
      enddo
      if (ep .eq. -2) goto 61
      do j=0,ep+2
      epmj=ep-j

      Cv(N+cc0011111,ep)=Cv(N+cc0011111,ep)-idp3(j)*(
     & +half*b1111sum(epmj)
     & +2*b111sum(epmj)
     & +3*b11sum(epmj)
     & +2*b1sum(epmj)
     & +half*b0sum(epmj)
     & +half*f1*Cv(N+cc111111,epmj)
     & +half*f2*Cv(N+cc111112,epmj)
     & -Cv(N+cc11111,epmj)*m1s)


      Cv(N+cc0011112,ep)=Cv(N+cc0011112,ep)-idp3(j)*(
     &  -half*b1sum(epmj)
     &  -3d0*half*b11sum(epmj)
     &  -3d0*half*b111sum(epmj)
     &  -half*b1111sum(epmj)
     &  +half*f1*Cv(N+cc111112,epmj)
     &  +half*f2*Cv(N+cc111122,epmj)
     &  -Cv(N+cc11112,epmj)*m1s)


      Cv(N+cc0012222,ep)=Cv(N+cc0012222,ep)-idp3(j)*(
     &  + half*Bv(bb1111+B23,epmj)
     &  + half*Bv(bb11111+B23,epmj)
     &  + half*f1*Cv(N+cc112222,epmj)
     &  + half*f2*Cv(N+cc122222,epmj)
     &  - Cv(N+cc12222,epmj)*m1s)


      Cv(N+cc0022222,ep)=Cv(N+cc0022222,ep)-idp3(j)*(
     & - half*Bv(bb11111+B23,epmj)
     &  + half*f1*Cv(N+cc122222,epmj)
     &  + half*f2*Cv(N+cc222222,epmj)
     &  - Cv(N+cc22222,epmj)*m1s)


      Cv(N+cc0011222,ep)=Cv(N+cc0011222,ep)-idp3(j)*(
     & - half*Bv(bb111+B23,epmj)
     &  - Bv(bb1111+B23,epmj)
     &  - half*Bv(bb11111+B23,epmj)
     &  + half*f1*Cv(N+cc111222,epmj)
     &  + half*f2*Cv(N+cc112222,epmj)
     &  - Cv(N+cc11222,epmj)*m1s)


      Cv(N+cc0011122,ep)=Cv(N+cc0011122,ep)-idp3(j)*(
     &  + half*Bv(bb11+B23,epmj)
     &  + 3d0*half*Bv(bb111+B23,epmj)
     &  + 3d0*half*Bv(bb1111+B23,epmj)
     &  + half*Bv(bb11111+B23,epmj)
     &  + half*f1*Cv(N+cc111122,epmj)
     &  + half*f2*Cv(N+cc111222,epmj)
     &  - Cv(N+cc11122,epmj)*m1s)


      Cv(N+cc0000001,ep)=Cv(N+cc0000001,ep)-idp3(j)*(
     &  + half*Bv(bb0000+B23,epmj)
     &  + half*Bv(bb00001+B23,epmj)
     &  + half*f1*Cv(N+cc000011,epmj)
     &  + half*f2*Cv(N+cc000012,epmj)
     &  - Cv(N+cc00001,epmj)*m1s)

      Cv(N+cc0000002,ep)=Cv(N+cc0000002,ep)-idp3(j)*(
     &  - half*Bv(bb00001+B23,epmj)
     &  + half*f1*Cv(N+cc000012,epmj)
     &  + half*f2*Cv(N+cc000022,epmj)
     &  - Cv(N+cc00002,epmj)*m1s)

      Cv(N+cc0000111,ep)=Cv(N+cc0000111,ep)-idp3(j)*(
     &  + half*Bv(bb00+B23,epmj)
     &  + 3d0*half*Bv(bb001+B23,epmj)
     &  + 3d0*half*Bv(bb0011+B23,epmj)
     &  + half*Bv(bb00111+B23,epmj)
     &  + half*f1*Cv(N+cc001111,epmj)
     &  + half*f2*Cv(N+cc001112,epmj)
     &  - Cv(N+cc00111,epmj)*m1s)

      Cv(N+cc0000112,ep)=Cv(N+cc0000112,ep)-idp3(j)*(
     &  - half*Bv(bb001+B23,epmj)
     &  - Bv(bb0011+B23,epmj)
     &  - half*Bv(bb00111+B23,epmj)
     &  + half*f1*Cv(N+cc001112,epmj)
     &  + half*f2*Cv(N+cc001122,epmj)
     &  - Cv(N+cc00112,epmj)*m1s)

      Cv(N+cc0000122,ep)=Cv(N+cc0000122,ep)-idp3(j)*(
     & + half*Bv(bb0011+B23,epmj)
     &  + half*Bv(bb00111+B23,epmj)
     &  + half*f1*Cv(N+cc001122,epmj)
     &  + half*f2*Cv(N+cc001222,epmj)
     &  - Cv(N+cc00122,epmj)*m1s)


      Cv(N+cc0000222,ep)=Cv(N+cc0000222,ep)-idp3(j)*(
     &  - half*Bv(bb00111+B23,epmj)
     &  + half*f1*Cv(N+cc001222,epmj)
     &  + half*f2*Cv(N+cc002222,epmj)
     &  - Cv(N+cc00222,epmj)*m1s
     &  )

      enddo
 61   continue

c      write(66,*) 'zx',cc0000001,ep,Cv(N+cc0000001,ep)
c      write(66,*) 'zx',cc0000002,ep,Cv(N+cc0000002,ep)
c      write(66,*) 'zx',cc0000111,ep,Cv(N+cc0000111,ep)
c      write(66,*) 'zx',cc0000112,ep,Cv(N+cc0000112,ep)
c      write(66,*) 'zx',cc0000122,ep,Cv(N+cc0000122,ep)
c      write(66,*) 'zx',cc0000222,ep,Cv(N+cc0000222,ep)

c      write(66,*) 'zx',cc0011111,ep,Cv(N+cc0011111,ep)
c      write(66,*) 'zx',cc0011112,ep,Cv(N+cc0011112,ep)
c      write(66,*) 'zx',cc0011122,ep,Cv(N+cc0011122,ep)
c      write(66,*) 'zx',cc0011222,ep,Cv(N+cc0011222,ep)
c      write(66,*) 'zx',cc0012222,ep,Cv(N+cc0012222,ep)
c      write(66,*) 'zx',cc0022222,ep,Cv(N+cc0022222,ep)

      enddo
      

C   Cv(N+cc0000001,Cv(N+cc0000002,zzzzzz)
      do ep=-2,0
      in(1,ep)=
     &  + f1*Cv(N+cc000000,ep)
     &  + Bv(bb000000+B13,ep)
     &  - Bv(bb000000+B23,ep)
      in(2,ep)=
     &  + f2*Cv(N+cc000000,ep)
     &  + Bv(bb000000+B12,ep)
     &  - Bv(bb000000+B23,ep)
      enddo

      call pvBackSubst(G,2,perm, in)

      do ep=-2,0
      Cv(N+cc0000001,ep)=in(1,ep)
      Cv(N+cc0000002,ep)=in(2,ep)
c      write(66,*) 'nx',cc0000001,ep,Cv(N+cc0000001,ep)
c      write(66,*) 'nx',cc0000002,ep,Cv(N+cc0000002,ep)


C   Cv(N+cc0000111,Cv(N+cc0000112,zzzzpp)
      in(1,ep)=
     &  - Bv(bb0000+B23,ep)
     &  - 2*Bv(bb00001+B23,ep)
     &  + f1*Cv(N+cc000011,ep)
     &  - 4*Cv(N+cc0000001,ep)
     &  - Bv(bb000011+B23,ep)
      in(2,ep)=
     &  - Bv(bb0000+B23,ep)
     &  - 2*Bv(bb00001+B23,ep)
     &  + f2*Cv(N+cc000011,ep)
     &  + Bv(bb000011+B12,ep)
     &  - Bv(bb000011+B23,ep)
      enddo


      call pvBackSubst(G,2,perm, in)

      do ep=-2,0
      Cv(N+cc0000111,ep)=in(1,ep)
      Cv(N+cc0000112,ep)=in(2,ep)
c      write(66,*) 'nx',cc0000111,ep,Cv(N+cc0000111,ep)
c      write(66,*) 'nx',cc0000112,ep,Cv(N+cc0000112,ep)

C   Cv(N+cc0000112,Cv(N+cc0000122,zzzzpk)
      in(1,ep)=
     &  + Bv(bb00001+B23,ep)
     &  + f1*Cv(N+cc000012,ep)
     &  - 2*Cv(N+cc0000002,ep)
     &  + Bv(bb000011+B23,ep)
      in(2,ep)=
     &  + Bv(bb00001+B23,ep)
     &  + f2*Cv(N+cc000012,ep)
     &  - 2*Cv(N+cc0000001,ep)
     &  + Bv(bb000011+B23,ep)
      enddo



      call pvBackSubst(G,2,perm, in)

      do ep=-2,0
      Cv(N+cc0000112,ep)=in(1,ep)
      Cv(N+cc0000122,ep)=in(2,ep)
c      write(66,*) 'nx',cc0000112,ep,Cv(N+cc0000112,ep)
c      write(66,*) 'nx',cc0000122,ep,Cv(N+cc0000122,ep)

C   Cv(N+cc0000122,Cv(N+cc0000222,zzzzkk)
      in(1,ep)=
     &  + f1*Cv(N+cc000022,ep)
     &  + Bv(bb000011+B13,ep)
     &  - Bv(bb000011+B23,ep)
      in(2,ep)=
     &  + f2*Cv(N+cc000022,ep)
     &  - 4*Cv(N+cc0000002,ep)
     &  - Bv(bb000011+B23,ep)
      enddo

      call pvBackSubst(G,2,perm, in)

      do ep=-2,0
      Cv(N+cc0000122,ep)=in(1,ep)
      Cv(N+cc0000222,ep)=in(2,ep)
c      write(66,*) 'nx',cc0000122,ep,Cv(N+cc0000122,ep)
c      write(66,*) 'nx',cc0000222,ep,Cv(N+cc0000222,ep)

C   Cv(N+cc0011111,Cv(N+cc0011112,zzpppp)
      in(1,ep)=
     &  - Bv(bb00+B23,ep)
     &  - 4*Bv(bb001+B23,ep)
     &  - 6*Bv(bb0011+B23,ep)
     &  - 4*Bv(bb00111+B23,ep)
     &  + f1*Cv(N+cc001111,ep)
     &  - 8*Cv(N+cc0000111,ep)
     &  - Bv(bb001111+B23,ep)
Cv(N+cc(zzpppp,2,P?,K?,m1?,m2?,m3?)=(
      in(2,ep)=
     &  - Bv(bb00+B23,ep)
     &  - 4*Bv(bb001+B23,ep)
     &  - 6*Bv(bb0011+B23,ep)
     &  - 4*Bv(bb00111+B23,ep)
     &  + f2*Cv(N+cc001111,ep)
     &  + Bv(bb001111+B12,ep)
     &  - Bv(bb001111+B23,ep)
      enddo

      call pvBackSubst(G,2,perm, in)

      do ep=-2,0
      Cv(N+cc0011111,ep)=in(1,ep)
      Cv(N+cc0011112,ep)=in(2,ep)
c      write(66,*) 'nx',cc0011111,ep,Cv(N+cc0011111,ep)
c      write(66,*) 'nx',cc0011112,ep,Cv(N+cc0011112,ep)

C   Cv(N+cc0011112,Cv(N+cc0011122,zzpppk)
      in(1,ep)=
     &  + Bv(bb001+B23,ep)
     &  + 3*Bv(bb0011+B23,ep)
     &  + 3*Bv(bb00111+B23,ep)
     &  + f1*Cv(N+cc001112,ep)
     &  - 6*Cv(N+cc0000112,ep)
     &  + Bv(bb001111+B23,ep)
      in(2,ep)=
     &  + Bv(bb001+B23,ep)
     &  + 3*Bv(bb0011+B23,ep)
     &  + 3*Bv(bb00111+B23,ep)
     &  + f2*Cv(N+cc001112,ep)
     &  - 2*Cv(N+cc0000111,ep)
     &  + Bv(bb001111+B23,ep)
      enddo


      call pvBackSubst(G,2,perm, in)

      do ep=-2,0
      Cv(N+cc0011112,ep)=in(1,ep)
      Cv(N+cc0011122,ep)=in(2,ep)
c      write(66,*) 'nx',cc0011112,ep,Cv(N+cc0011112,ep)
c      write(66,*) 'nx',cc0011122,ep,Cv(N+cc0011122,ep)

C   Cv(N+cc0011122,Cv(N+cc0011222,zzppkk)
      in(1,ep)=
     &  - Bv(bb0011+B23,ep)
     &  - 2*Bv(bb00111+B23,ep)
     &  + f1*Cv(N+cc001122,ep)
     &  - 4*Cv(N+cc0000122,ep)
     &  - Bv(bb001111+B23,ep)
      in(2,ep)=
     &  - Bv(bb0011+B23,ep)
     &  - 2*Bv(bb00111+B23,ep)
     &  + f2*Cv(N+cc001122,ep)
     &  - 4*Cv(N+cc0000112,ep)
     &  - Bv(bb001111+B23,ep)
      enddo


      call pvBackSubst(G,2,perm, in)

      do ep=-2,0
      Cv(N+cc0011122,ep)=in(1,ep)
      Cv(N+cc0011222,ep)=in(2,ep)
c      write(66,*) 'nx',cc0011122,ep,Cv(N+cc0011122,ep)
c      write(66,*) 'nx',cc0011222,ep,Cv(N+cc0011222,ep)

C   Cv(N+cc0011222,Cv(N+cc0012222,zzpkkk)
      in(1,ep)=
     &  + Bv(bb00111+B23,ep)
     &  + f1*Cv(N+cc001222,ep)
     &  - 2*Cv(N+cc0000222,ep)
     &  + Bv(bb001111+B23,ep)
      in(2,ep)=
     &  + Bv(bb00111+B23,ep)
     &  + f2*Cv(N+cc001222,ep)
     &  - 6*Cv(N+cc0000122,ep)
     &  + Bv(bb001111+B23,ep)
      enddo


      call pvBackSubst(G,2,perm, in)

      do ep=-2,0
      Cv(N+cc0011222,ep)=in(1,ep)
      Cv(N+cc0012222,ep)=in(2,ep)
c      write(66,*) 'nx',cc0011222,ep,Cv(N+cc0011222,ep)
c      write(66,*) 'nx',cc0012222,ep,Cv(N+cc0012222,ep)

C   Cv(N+cc0012222,Cv(N+cc0022222,zzkkkk)
      in(1,ep)=
     &  + f1*Cv(N+cc002222,ep)
     &  + Bv(bb001111+B13,ep)
     &  - Bv(bb001111+B23,ep)
      in(2,ep)=
     &  + f2*Cv(N+cc002222,ep)
     &  - 8*Cv(N+cc0000222,ep)
     &  - Bv(bb001111+B23,ep)
      enddo


      call pvBackSubst(G,2,perm, in)

      do ep=-2,0
      Cv(N+cc0012222,ep)=in(1,ep)
      Cv(N+cc0022222,ep)=in(2,ep)
c      write(66,*) 'nx',cc0012222,ep,Cv(N+cc0012222,ep)
c      write(66,*) 'nx',cc0022222,ep,Cv(N+cc0022222,ep)


C   Cv(N+cc1111111,Cv(N+cc1111112,pppppp)
      in(1,ep)=
     &  + f1*Cv(N+cc111111,ep)
     &  - 12*Cv(N+cc0011111,ep)
     &  - b11111sum(ep)
     &  - 5*b1111sum(ep)
     &  - 10*b111sum(ep)
     &  - 10*b11sum(ep)
     &  - 5*b1sum(ep)
     &  - b0sum(ep)
      in(2,ep)=
     &  + f2*Cv(N+cc111111,ep)
     &  - b11111sum(ep)
     &  - 5*b1111sum(ep)
     &  - 10*b111sum(ep)
     &  - 10*b11sum(ep)
     &  - 5*b1sum(ep)
     &  - b0sum(ep)
     &  + Bv(bb111111+B12,ep)
      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc1111111,ep)=in(1,ep)
      Cv(N+cc1111112,ep)=in(2,ep)
c      write(66,*) 'nx',cc1111111,ep,Cv(N+cc1111111,ep)
c      write(66,*) 'nx',cc1111112,ep,Cv(N+cc1111112,ep)

      
C   Cv(N+cc1111112,Cv(N+cc1111122,pppppk)
      in(1,ep)=
     &  + b1sum(ep)
     &  + 4*b11sum(ep)
     &  + 6*b111sum(ep)
     &  + 4*b1111sum(ep)
     &  + b11111sum(ep)
     &  + f1*Cv(N+cc111112,ep)
     &  - 10*Cv(N+cc0011112,ep)
      in(2,ep)=
     &  + b1sum(ep)
     &  + 4*b11sum(ep)
     &  + 6*b111sum(ep)
     &  + 4*b1111sum(ep)
     &  + b11111sum(ep)
     &  + f2*Cv(N+cc111112,ep)
     &  - 2*Cv(N+cc0011111,ep)
      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc1111112,ep)=in(1,ep)
      Cv(N+cc1111122,ep)=in(2,ep)
c      write(66,*) 'nx',cc1111112,ep,Cv(N+cc1111112,ep)
c      write(66,*) 'nx',cc1111122,ep,Cv(N+cc1111122,ep)

C   Cv(N+cc1111122,Cv(N+cc1111222,ppppkk)
      in(1,ep)=
     &  - b11sum(ep)
     &  - 3*b111sum(ep)
     &  - 3*b1111sum(ep)
     &  - b11111sum(ep)
     &  + f1*Cv(N+cc111122,ep)
     &  - 8*Cv(N+cc0011122,ep)
      in(2,ep)=
     &  - b11sum(ep)
     &  - 3*b111sum(ep)
     &  - 3*b1111sum(ep)
     &  - b11111sum(ep)
     &  + f2*Cv(N+cc111122,ep)
     &  - 4*Cv(N+cc0011112,ep)
      enddo
      call pvBackSubst(G,2,perm, in)
      do ep=-2,0
      Cv(N+cc1111122,ep)=in(1,ep)
      Cv(N+cc1111222,ep)=in(2,ep)
c      write(66,*) 'nx',cc1111122,ep,Cv(N+cc1111122,ep)
c      write(66,*) 'nx',cc1111222,ep,Cv(N+cc1111222,ep)

C   Cv(N+cc1111222,Cv(N+cc1112222,pppkkk)
      in(1,ep)=
     &  + b111sum(ep)
     &  + 2*b1111sum(ep)
     &  + b11111sum(ep)
     &  + f1*Cv(N+cc111222,ep)
     &  - 6*Cv(N+cc0011222,ep)
      in(2,ep)=
     &  + b111sum(ep)
     &  + 2*b1111sum(ep)
     &  + b11111sum(ep)
     &  + f2*Cv(N+cc111222,ep)
     &  - 6*Cv(N+cc0011122,ep)
      enddo

 
      call pvBackSubst(G,2,perm, in)

      do ep=-2,0
      Cv(N+cc1111222,ep)=in(1,ep)
      Cv(N+cc1112222,ep)=in(2,ep)
c      write(66,*) 'nx',cc1111222,ep,Cv(N+cc1111222,ep)
c      write(66,*) 'nx',cc1112222,ep,Cv(N+cc1112222,ep)

C   Cv(N+cc1112222,Cv(N+cc1122222,ppkkkk)
      in(1,ep)=
     &  - b1111sum(ep)
     &  - b11111sum(ep)
     &  + f1*Cv(N+cc112222,ep)
     &  - 4*Cv(N+cc0012222,ep)
      in(2,ep)=
     &  - b1111sum(ep)
     &  - b11111sum(ep)
     &  + f2*Cv(N+cc112222,ep)
     &  - 8*Cv(N+cc0011222,ep)
      enddo
      call pvBackSubst(G,2,perm, in)

      do ep=-2,0
      Cv(N+cc1112222,ep)=in(1,ep)
      Cv(N+cc1122222,ep)=in(2,ep)
c      write(66,*) 'nx',cc1112222,ep,Cv(N+cc1112222,ep)
c      write(66,*) 'nx',cc1122222,ep,Cv(N+cc1122222,ep)

C   Cv(N+cc1122222,Cv(N+cc1222222,pkkkkk)
      in(1,ep)=
     &  + f1*Cv(N+cc122222,ep)
     &  - 2*Cv(N+cc0022222,ep)
     &  + b11111sum(ep)
      in(2,ep)=
     &  + f2*Cv(N+cc122222,ep)
     &  - 10*Cv(N+cc0012222,ep)
     &  + b11111sum(ep)
      enddo

      call pvBackSubst(G,2,perm, in)

      do ep=-2,0
      Cv(N+cc1122222,ep)=in(1,ep)
      Cv(N+cc1222222,ep)=in(2,ep)
c      write(66,*) 'nx',cc1122222,ep,Cv(N+cc1122222,ep)
c      write(66,*) 'nx',cc1222222,ep,Cv(N+cc1222222,ep)

C   Cv(N+cc1222222,Cv(N+cc2222222,kkkkkk)
      in(1,ep)=
     &  + f1*Cv(N+cc222222,ep)
     &  + Bv(bb111111+B13,ep)
     &  - Bv(bb111111+B23,ep)
      in(2,ep)=
     &  + f2*Cv(N+cc222222,ep)
     &  - 12*Cv(N+cc0022222,ep)
     &  - Bv(bb111111+B23,ep)
      enddo

      call pvBackSubst(G,2,perm, in)

      do ep=-2,0
      Cv(N+cc1222222,ep)=in(1,ep)
      Cv(N+cc2222222,ep)=in(2,ep)
c      write(66,*) 'nx',cc1222222,ep,Cv(N+cc1222222,ep)
c      write(66,*) 'nx',cc2222222,ep,Cv(N+cc2222222,ep)
      enddo


c--- to check recursion identities      
c      call Cfill_alt(p1,p2,p1p2,m1s,m2s,m3s,N)
       
      return
      end
