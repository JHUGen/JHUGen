      subroutine amp_zqqQQa_ql(i1,i2,i3,i4,i5,i6,i7,ai)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: j,k,j1,j2,j3,j4,j5,j6,j7,i1,i2,i3,i4,i5,i6,i7
      integer:: hq,Qh,lh,hg,h2,h3,h4
      integer:: hqc,Qhc,lhc,hgc
      real(dp):: t
      complex(dp):: ai(2,2,2,2),A7h1,A7h2,A7h3,A7h4
      complex(dp):: xa(mxpart,mxpart),xb(mxpart,mxpart)
C-----statement function
      j1=i1
      j2=i2
      j5=i5
C----hq,Qh,hg,lh
      do hq=2,2
      do Qh=1,2
      do hg=1,2
      do lh=1,2
         if (hq == 2) then
            do j=1,mxpart
            do k=1,mxpart
            xa(j,k)=za(j,k)
            xb(j,k)=zb(j,k)
c            h1=2
            h2=Qh
            h3=hg
            h4=lh
            enddo
            enddo
         elseif (hq == 1) then
            do j=1,mxpart
            do k=1,mxpart
            xa(j,k)=zb(k,j)
            xb(j,k)=za(k,j)
            enddo
            enddo
c            h1=3-hq
            h2=3-Qh
            h3=3-hg
            h4=3-lh
         endif
         if (h2 == 2) then
            j3=i3
            j4=i4
         elseif (h2 == 1) then
            j3=i4
            j4=i3
         endif
         if (h4 == 2) then
            j6=i6
            j7=i7
         elseif (h4 == 1) then
            j6=i7
            j7=i6
         endif
c-----
      if (h3 == 2) then
c-----A(2,2,2,2)
      A7h1=czip
      A7h1=1/(xa(j5,j6)*xa(j5,j7))*( 
     .-xb(j1,j3)*xa(j2,j6)*( xa(j4,j1)*xb(j1,j5)*xa(j5,j6) +
     &  xa(j4,j1)*xb(j1,j7)*xa(j7,j6) + xa(j4,j3)*xb(j3,j5)*xa(j5,j6) +
     &  xa(j4,j3)*xb(j3,j7)*xa(j7,j6) )/t(j1,j3,j4)
     .-xa(j4,j2)*( xa(j6,j2)*xb(j2,j3) + xa(j6,j4)*xb(j4,j3) )*
     &  ( xa(j6,j5)*xb(j5,j1) + xa(j6,j7)*xb(j7,j1) )/t(j2,j3,j4) )
      ai(hq,Qh,hg,lh) = A7h1
      elseif (h3 == 1)  then
c-----A(2,2,1,2)
      A7h4=czip
      A7h4=1/(xb(j5,j6)*xb(j5,j7))*( 
     .-xb(j1,j3)*( xa(j4,j1)*xb(j1,j7) + xa(j4,j3)*xb(j3,j7) )*
     &  ( xa(j2,j5)*xb(j5,j7) + xa(j2,j6)*xb(j6,j7) )/t(j1,j3,j4) 
     .+xa(j4,j2)*xb(j7,j1)*( xb(j3,j2)*xa(j2,j5)*xb(j5,j7) +
     &  xb(j3,j2)*xa(j2,j6)*xb(j6,j7) + xb(j3,j4)*xa(j4,j5)*xb(j5,j7) +
     &  xb(j3,j4)*xa(j4,j6)*xb(j6,j7) )/t(j2,j3,j4) )
      ai(hq,Qh,hg,lh) = A7h4
      endif
c-----put in photon propagator
      ai(hq,Qh,hg,lh)=ai(hq,Qh,hg,lh)/t(j5,j6,j7)/s(j3,j4)
c--------obtain hq=1 from complex conjugation
      hqc=mod(hq+2,2)+1
      Qhc=mod(Qh+2,2)+1
      hgc=mod(hg+2,2)+1
      lhc=mod(lh+2,2)+1
      ai(hqc,Qhc,hgc,lhc)=conjg(ai(hq,Qh,hg,lh))
c-----
      enddo
      enddo
      enddo
      enddo
c-----done
      return
      end
