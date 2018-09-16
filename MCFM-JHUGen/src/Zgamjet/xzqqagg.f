      subroutine xzqqagg_qq(j1,j2,j3,j4,j5,j6,j7,a70h)
      implicit none
      include 'types.f'
      
************************************************************************
*     return the helicity amplitudes for                               *
*     0 ---> q(p1)+ph(p2)+g(p3)+g(p4)+qbar(p5)+lb(p6)+l(p7)            *
*     photon is coming from quark line                                 *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: i2(6),i3(6),i4(6),j,lh,h2,h3,h4,hq,h(7)
      integer:: j1,j2,j3,j4,j5,j6,j7,icol
      integer:: hqc,lhc,h2c,h3c,h4c
      complex(dp):: a70h(2,2,2,2,2,2)
      complex(dp):: atemp,m(6),amp_qqggg,a34,a43
c-----(2,3,4) permutation
      i2(1)=j2
      i3(1)=j3
      i4(1)=j4
      i2(2)=j2
      i3(2)=j4
      i4(2)=j3
      i2(3)=j4
      i3(3)=j2
      i4(3)=j3
      i2(4)=j3
      i3(4)=j4
      i4(4)=j2
      i2(5)=j3
      i3(5)=j2
      i4(5)=j4
      i2(6)=j4
      i3(6)=j3
      i4(6)=j2
c-----initialize helicity amplitudes to zero
c-----a70h(hq,hg2,hg3,hg4,hl)
      do hq=1,2
      do lh=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do icol=1,2
      a70h(icol,hq,h2,h3,h4,lh)=czip
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c-----fill the helicity amplitude
      do hq=2,2
      do lh=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
         h(j2)=h2
         h(j3)=h3
         h(j4)=h4
         do j=1,6
            m(j)=amp_qqggg(j1,hq,i2(j),h(i2(j)),i3(j),h(i3(j)),
     &                     i4(j),h(i4(j)),j5,lh,j6,j7)
         enddo
         a70h(1,hq,h2,h3,h4,lh)=m(1)+m(5)+m(4)
         a70h(2,hq,h2,h3,h4,lh)=m(2)+m(3)+m(6)
c--------obtain hq=1 from complex conjugation
         hqc=mod(hq+2,2)+1
         lhc=mod(lh+2,2)+1
         h2c=mod(h2+2,2)+1
         h3c=mod(h3+2,2)+1
         h4c=mod(h4+2,2)+1
         a70h(1,hqc,h2c,h3c,h4c,lhc)=conjg(a70h(1,hq,h2,h3,h4,lh))
         a70h(2,hqc,h2c,h3c,h4c,lhc)=conjg(a70h(2,hq,h2,h3,h4,lh))
      enddo
      enddo
      enddo
      enddo
      enddo
c-----done     
      return
      end

      subroutine xzqqagg_ql(j1,j2,j3,j4,j5,j6,j7,a70h3)
      implicit none
      include 'types.f'
************************************************************************
*     return the helicity amplitudes for                               *
*     0 ---> q(p1)+ph(p2)+g(p3)+g(p4)+qbar(p5)+lb(p6)+l(p7)            *
*     photon is coming from lepton line                                *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: i2(4),i3(4),i4(4),j,lh,h2,h3,h4,hq,h(7)
      integer:: j1,j2,j3,j4,j5,j6,j7,icol
      integer:: hqc,lhc,h2c,h3c,h4c
      complex(dp):: a70h3(2,2,2,2,2,2)
      complex(dp):: amp_qqaag_ql,amp_qqagg_ql
      complex(dp):: m3(2)
      real(dp):: cf1,cf2
      integer:: ii
c-----(3,4) permutation
      i2(1)=j2
      i3(1)=j3
      i4(1)=j4
      i2(2)=j2
      i3(2)=j4
      i4(2)=j3
c initialize helicity amplitudes to zero
c a70h(hq,hg2,hg3,hg4,hl)
      do hq=1,2
      do lh=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do icol=1,2
      a70h3(icol,hq,h2,h3,h4,lh)=czip
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c fill the helicity amplitude
      do hq=2,2
      do lh=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
         h(j2)=h2
         h(j3)=h3
         h(j4)=h4
c--------symmetrize over 3 and 4, 2 is coming from lepton line
         do j=1,2
            m3(j)=amp_qqagg_ql(j1,hq,i2(j),h(i2(j)),i3(j),h(i3(j)),
     &                         i4(j),h(i4(j)),j5,lh,j6,j7)
         enddo
         a70h3(1,hq,h2,h3,h4,lh)=m3(1)
         a70h3(2,hq,h2,h3,h4,lh)=m3(2)
c--------obtain hq=1 from complex conjugation
         hqc=mod(hq+2,2)+1
         lhc=mod(lh+2,2)+1
         h2c=mod(h2+2,2)+1
         h3c=mod(h3+2,2)+1
         h4c=mod(h4+2,2)+1
         a70h3(1,hqc,h2c,h3c,h4c,lhc)=conjg(a70h3(1,hq,h2,h3,h4,lh))
         a70h3(2,hqc,h2c,h3c,h4c,lhc)=conjg(a70h3(2,hq,h2,h3,h4,lh))
      enddo
      enddo
      enddo
      enddo
      enddo
c-----done
      return
      end
 
