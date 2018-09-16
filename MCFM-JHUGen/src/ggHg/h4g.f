      subroutine h4g(p1,p2,p3,p4,Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      integer:: j,p1,p2,p3,p4,h1,h2,h3,h4
      real(dp):: Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625
      complex(dp):: amp(3,2,2,2,2),
     &  amppp(3),apmpp(3),appmp(3),apppm(3),
     &  apppp(3),
     &  ammpp(3),ampmp(3),amppm(3),apmmp(3),apmpm(3),appmm(3)

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do j=1,3
      amp(j,h1,h2,h3,h4)=czip
      enddo
      enddo
      enddo
      enddo
      enddo

      call makepppp(p1,p2,p3,p4,za,apppp)
      call makemppp(p1,p2,p3,p4,za,zb,amppp,apmpp,appmp,apppm)
      call makemmpp(p1,p2,p3,p4,za,zb,
     & ammpp,ampmp,amppm,apmmp,apmpm,appmm)



      do j=1,3
      amp(j,2,2,2,2)=apppp(j)
      amp(j,1,2,2,2)=amppp(j)
      amp(j,2,1,2,2)=apmpp(j)
      amp(j,2,2,1,2)=appmp(j)
      amp(j,2,2,2,1)=apppm(j)

      amp(j,1,1,2,2)=ammpp(j)
      amp(j,1,2,1,2)=ampmp(j)
      amp(j,1,2,2,1)=amppm(j)
      amp(j,2,1,1,2)=apmmp(j)
      amp(j,2,1,2,1)=apmpm(j)
      amp(j,2,2,1,1)=appmm(j)
      enddo


c      call makepppp(p1,p2,p3,p4,zb,apppp)
c      call makemppp(p1,p2,p3,p4,zb,za,amppp,apmpp,appmp,apppm)
c      call makemmpp(p1,p2,p3,p4,zb,za,
c     & ammpp,ampmp,amppm,apmmp,apmpm,appmm)

c      do j=1,3
c      amp(j,1,1,1,1)=apppp(j)
c      amp(j,2,1,1,1)=amppp(j)
c      amp(j,1,2,1,1)=apmpp(j)
c      amp(j,1,1,2,1)=appmp(j)
c      amp(j,1,1,1,2)=apppm(j)

c      amp(j,1,1,2,2)=ammpp(j)
c      amp(j,1,2,1,2)=ampmp(j)
c      amp(j,1,2,2,1)=amppm(j)
c      amp(j,2,1,1,2)=apmmp(j)
c      amp(j,2,1,2,1)=apmpm(j)
c      amp(j,2,2,1,1)=appmm(j)

c      enddo

      do j=1,3
      amp(j,1,1,1,1)=conjg(amp(j,2,2,2,2))
      amp(j,2,1,1,1)=conjg(amp(j,1,2,2,2))
      amp(j,1,2,1,1)=conjg(amp(j,2,1,2,2))
      amp(j,1,1,2,1)=conjg(amp(j,2,2,1,2))
      amp(j,1,1,1,2)=conjg(amp(j,2,2,2,1))

      amp(j,1,1,2,2)=conjg(amp(j,2,2,1,1))
      amp(j,1,2,1,2)=conjg(amp(j,2,1,2,1))
      amp(j,1,2,2,1)=conjg(amp(j,2,1,1,2))
      amp(j,2,1,1,2)=conjg(amp(j,1,2,2,1))
      amp(j,2,1,2,1)=conjg(amp(j,1,2,1,2))
      amp(j,2,2,1,1)=conjg(amp(j,1,1,2,2))
      enddo

      Hgggg_1256=0._dp
      Hgggg_1265=0._dp
      Hgggg_1625=0._dp

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
C         write(*,*) 'h4g:  ',h1,h2,h3,h4
C         write(*,*) 'h4g:  ',amp(1,h1,h2,h3,h4),amp(2,h1,h2,h3,h4),
C     &        amp(3,h1,h2,h3,h4)
C         write(*,*) 'h4gsq:',xn**2*V/2._dp*abs(amp(1,h1,h2,h3,h4))**2
C     &        ,xn**2*V/2._dp*abs(amp(2,h1,h2,h3,h4))**2,
C     &        xn**2*V/2._dp*abs(amp(3,h1,h2,h3,h4))**2
      Hgggg_1256=Hgggg_1256+abs(amp(1,h1,h2,h3,h4))**2
      Hgggg_1265=Hgggg_1265+abs(amp(2,h1,h2,h3,h4))**2
      Hgggg_1625=Hgggg_1625+abs(amp(3,h1,h2,h3,h4))**2
      enddo
      enddo
      enddo
      enddo

C===  (1/4 ---> 1/2) because only three orderings)
      Hgggg_1256=xn**2*V/2._dp*Hgggg_1256
      Hgggg_1265=xn**2*V/2._dp*Hgggg_1265
      Hgggg_1625=xn**2*V/2._dp*Hgggg_1625

      Hgggg=Hgggg_1256+Hgggg_1265+Hgggg_1625

      return
      end


      subroutine makepppp(p1,p2,p3,p4,za,apppp)
      implicit none
      include 'types.f'

C     Taken from Kauffman hep-ph/9903330
C     and (older formula)
C     %\cite{Kauffman:1996ix}
C     \bibitem{Kauffman:1996ix}
C     R.~P.~Kauffman, S.~V.~Desai and D.~Risal,
C     %``Production of a Higgs boson plus two jets in hadronic collisions,''
C     Phys.\ Rev.\ D {\bf 55}, 4005 (1997)
C     [Erratum-ibid.\ D {\bf 58}, 119901 (1998)]
C     [arXiv:hep-ph/9610541].
C     %%CITATION = HEP-PH 9610541;%%
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j,p1,p2,p3,p4,i1(4),i2(4),i3(4),i4(4)
      complex(dp):: apppp(3)
      real(dp):: hm2

C if Higgs has non-zero width hm**2 must be recalculated
      hm2 = s(p1,p2)+s(p1,p3)+s(p1,p4)+s(p2,p3)+ s(p2,p4)+s(p3,p4)
C      write(*,*) 'hm2',hm2,hmass**2
      do j=1,2
      i1(j)=p1
            if (j==1) then
            i2(j)=p2
            i3(j)=p3
            i4(j)=p4
            elseif (j==2) then
            i2(j)=p2
            i3(j)=p4
            i4(j)=p3
c            elseif (j==3) then
c            i2(j)=p4
c            i3(j)=p2
c            i4(j)=p3
            endif
C---PRD55 Eq(21)
      apppp(j)=-cplx2(hm2**2,zip)/(za(i1(j),i2(j))*za(i2(j),i3(j))
     &  *za(i3(j),i4(j))*za(i4(j),i1(j)))
      enddo
C---determine apppp(3) using sub-cyclic identity
      apppp(3)=-apppp(1)-apppp(2)
      return
      end


      subroutine makemppp(p1,p2,p3,p4,za,zb,amppp,apmpp,appmp,apppm)
      implicit none
      include 'types.f'

C     Taken from Kauffman hep-ph/9903330
C     and (older formula)
C     %\cite{Kauffman:1996ix}
C     \bibitem{Kauffman:1996ix}
C     R.~P.~Kauffman, S.~V.~Desai and D.~Risal,
C     %``Production of a Higgs boson plus two jets in hadronic collisions,''
C     Phys.\ Rev.\ D {\bf 55}, 4005 (1997)
C     [Erratum-ibid.\ D {\bf 58}, 119901 (1998)]
C     [arXiv:hep-ph/9610541].
C     %%CITATION = HEP-PH 9610541;%%

      integer:: j,k,p1,p2,p3,p4,j1,j2,jk(4),
     & i1(4),i2(4),i3(4),i4(4)
      real(dp):: s123,s124,s134,s234
      complex(dp):: z2,amppp(3),apmpp(3),appmp(3),apppm(3),temp
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer, parameter:: k1(4)=(/1,2,3,4/)
      integer, parameter:: k2(4)=(/2,3,4,1/)
      integer, parameter:: k3(4)=(/3,4,1,2/)
      integer, parameter:: k4(4)=(/4,1,2,3/)


C---statement function
      z2(j1,j2)=-za(j1,p1)*zb(p1,j2)-za(j1,p2)*zb(p2,j2)
     &          -za(j1,p3)*zb(p3,j2)-za(j1,p4)*zb(p4,j2)
C---statement function
      jk(1)=p1
      jk(2)=p2
      jk(3)=p3
      jk(4)=p4
      do k=1,4
      do j=1,2
      i1(j)=jk(k1(k))
            if (j==1) then
            i2(j)=jk(k2(k))
            i3(j)=jk(k3(k))
            i4(j)=jk(k4(k))
      elseif (j==2) then
            i2(j)=jk(k2(k))
            i3(j)=jk(k4(k))
            i4(j)=jk(k3(k))
c      elseif (j==3) then
c            i2(j)=jk(k4(k))
c            i3(j)=jk(k2(k))
c            i4(j)=jk(k3(k))
      endif
      s124=s(i1(j),i2(j))+s(i1(j),i4(j))+s(i2(j),i4(j))
      s123=s(i1(j),i2(j))+s(i1(j),i3(j))+s(i2(j),i3(j))
      s134=s(i1(j),i3(j))+s(i1(j),i4(j))+s(i3(j),i4(j))
      s234=s(i2(j),i3(j))+s(i2(j),i4(j))+s(i3(j),i4(j))
C---PRD55 Eq(22)
c      amppp=
c     & -(z2(p1,p3)*zb(p2,p4))**2
c     & /((s(p1,p2)+s(p1,p4)+s(p2,p4))*s(p1,p2)*s(p1,p4))
c     & -(z2(p1,p4)*zb(p2,p3))**2/((s(p1,p2)+s(p1,p3)+s(p2,p3))*s(p1,p2)*s(p2,p3))
c     & -(z2(p1,p2)*zb(p3,p4))**2/((s(p1,p3)+s(p1,p4)+s(p3,p4))*s(p1,p4)*s(p3,p4))
c     & +zb(p2,p4)/(zb(p1,p2)*za(p2,p3)*za(p3,p4)*zb(p4,p1))
c     & *(+s(p2,p3)*z2(p1,p2)/za(p4,p1)
c     &   +s(p3,p4)*z2(p1,p4)/za(p1,p2)-zb(p2,p4)*(s(p2,p3)+s(p2,p4)+s(p3,p4)))
C---PRD55 Eq(A8+erratum)
c      amppp=
c     & -(z2(p1,p3)*zb(p2,p4))**2/((s(p1,p2)+s(p1,p4)+s(p2,p4))*s(p1,p2)*s(p1,p4))
c     & -(z2(p1,p4)*zb(p2,p3))**2/((s(p1,p2)+s(p1,p3)+s(p2,p3))*s(p1,p2)*s(p2,p3))
c     & -(z2(p1,p2)*zb(p3,p4))**2/((s(p1,p3)+s(p1,p4)+s(p3,p4))*s(p1,p4)*s(p3,p4))
c     & -zb(p2,p4)/(zb(p1,p2)*zb(p1,p4)*za(p1,p3))
c     & *(z2(p1,p2)**2/(za(p1,p4)*za(p3,p4))
c     &  +z2(p1,p4)**2/(za(p1,p2)*za(p2,p3)))

C---hep-ph/9903330 Eq(11)
      temp=
     & -(z2(i1(j),i3(j))*zb(i2(j),i4(j)))**2
     & /(s124*s(i1(j),i2(j))*s(i1(j),i4(j)))
     & -(z2(i1(j),i4(j))*zb(i2(j),i3(j)))**2
     & /(s123*s(i1(j),i2(j))*s(i2(j),i3(j)))
     & -(z2(i1(j),i2(j))*zb(i3(j),i4(j)))**2
     & /(s134*s(i1(j),i4(j))*s(i3(j),i4(j)))
     & +zb(i2(j),i4(j))/
     & (zb(i1(j),i2(j))*za(i2(j),i3(j))*za(i3(j),i4(j))*zb(i4(j),i1(j)))
     & *(s(i2(j),i3(j))*z2(i1(j),i2(j))/za(i4(j),i1(j))
     & +s(i3(j),i4(j))*z2(i1(j),i4(j))/za(i1(j),i2(j))
     & -zb(i2(j),i4(j))*s234)
C      if (k==1) amppp(j)=temp
C      if (k==2) apmpp(j)=temp
C      if (k==3) appmp(j)=temp
C      if (k==4) apppm(j)=temp


C -- GZ
      if (k==1) amppp(j)=temp
      if (k==2 .and. j==1) then
         apmpp(j)=temp
      elseif (k==2 .and. j==2) then
         apmpp(j+1)=temp
      endif
      if (k==3) appmp(j)=temp
      if (k==4 .and. j==1) then
         apppm(j)=temp
      elseif (k==4 .and. j==2) then
         apppm(j+1)=temp
      endif

C---determine axxxx(3) using sub-cyclic identity
      if (j ==2) then
         if (k==1) amppp(3)=-amppp(1)-amppp(2)
         if (k==2) apmpp(2)=-apmpp(1)-apmpp(3)
         if (k==3) appmp(3)=-appmp(1)-appmp(2)
         if (k==4) apppm(2)=-apppm(1)-apppm(3)
      endif

      enddo
      enddo
C---determine axxxx(3) using sub-cyclic identity
C      amppp(3)=-amppp(1)-amppp(2)
C      apmpp(3)=-apmpp(1)-apmpp(2)
C      appmp(3)=-appmp(1)-appmp(2)
C      apppm(3)=-apppm(1)-apppm(2)
      return
      end



      subroutine makemmpp(p1,p2,p3,p4,za,zb,
     & ammpp,ampmp,amppm,apmmp,apmpm,appmm)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
C     Taken from Kauffman hep-ph/9903330
C     and (older formula)
C     %\cite{Kauffman:1996ix}
C     \bibitem{Kauffman:1996ix}
C     R.~P.~Kauffman, S.~V.~Desai and D.~Risal,
C     %``Production of a Higgs boson plus two jets in hadronic collisions,''
C     Phys.\ Rev.\ D {\bf 55}, 4005 (1997)
C     [Erratum-ibid.\ D {\bf 58}, 119901 (1998)]
C     [arXiv:hep-ph/9610541].
C     %%CITATION = HEP-PH 9610541;%%
      integer:: j,k,p1,p2,p3,p4,jk(4),
     & i1(6),i2(6),i3(6),i4(6)
      complex(dp):: temp,
     & ammpp(3),ampmp(3),amppm(3),apmmp(3),apmpm(3),appmm(3)
      integer,parameter:: k1(6)=(/1,1,1,2,2,3/)
      integer,parameter:: k2(6)=(/2,3,4,3,4,4/)
      integer,parameter:: k3(6)=(/3,2,2,1,1,1/)
      integer,parameter:: k4(6)=(/4,4,3,4,3,2/)


      jk(1)=p1
      jk(2)=p2
      jk(3)=p3
      jk(4)=p4
      do k=1,6
            do j=1,2
            if (j==1) then
            i1(j)=jk(k1(k))
            i2(j)=jk(k2(k))
            i3(j)=jk(k3(k))
            i4(j)=jk(k4(k))
            elseif (j==2) then
            i1(j)=jk(k2(k))
            i2(j)=jk(k1(k))
            i3(j)=jk(k3(k))
            i4(j)=jk(k4(k))
c            elseif (j==3) then
c            i1(j)=jk(k2(k))
c            i2(j)=jk(k3(k))
c            i3(j)=jk(k1(k))
c            i4(j)=jk(k4(k))
            endif
            temp=
     &      -za(jk(k1(k)),jk(k2(k)))**4/(za(i1(j),i2(j))*za(i2(j),i3(j))
     &      *za(i3(j),i4(j))*za(i4(j),i1(j)))
     &      -zb(jk(k3(k)),jk(k4(k)))**4/(zb(i1(j),i2(j))*zb(i2(j),i3(j))
     &      *zb(i3(j),i4(j))*zb(i4(j),i1(j)))
C            if (k==1) ammpp(j)=temp
C            if (k==2) ampmp(j)=temp
C            if (k==3) amppm(j)=temp
C            if (k==4) apmmp(j)=temp
C            if (k==5) apmpm(j)=temp
C            if (k==6) appmm(j)=temp

            if (k==1) then
               if (j==1) ammpp(j)=temp
               if (j==2) then
                  ammpp(2)=temp
                  ammpp(3)=-ammpp(2)-ammpp(1)
               endif
            endif

            if (k==2) then
               if (j==1) ampmp(3)=temp
               if (j==2) then
                  ampmp(2)=temp
                  ampmp(1)=-ampmp(2)-ampmp(3)
               endif
            endif

            if (k==3) then
               if (j==1) amppm(3)=temp
               if (j==2) then
                  amppm(1)=temp
                  amppm(2)=-amppm(1)-amppm(3)
               endif
            endif

            if (k==4) then
               if (j==1) apmmp(3)=temp
               if (j==2) then
                  apmmp(1)=temp
                  apmmp(2)=-apmmp(1)-apmmp(3)
               endif
            endif

            if (k==5) then
               if (j==1) apmpm(3)=temp
               if (j==2) then
                  apmpm(2)=temp
                  apmpm(1)=-apmpm(2)-apmpm(3)
               endif
            endif

            if (k==6) then
               if (j==1) appmm(1)=temp
               if (j==2) then
                  appmm(2)=temp
                  appmm(3)=-appmm(1)-appmm(2)
               endif
            endif

      enddo
      enddo

C---determine axxxx(3) using sub-cyclic identity
C      ammpp(3)=-ammpp(1)-ammpp(2)
C      ampmp(3)=-ampmp(1)-ampmp(2)
C      amppm(3)=-amppm(1)-amppm(2)
C      apmmp(3)=-apmmp(1)-apmmp(2)
C      apmpm(3)=-apmpm(1)-apmpm(2)
C      appmm(3)=-appmm(1)-appmm(2)

      return
      end
