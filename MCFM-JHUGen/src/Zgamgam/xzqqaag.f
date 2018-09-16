      subroutine xzqqaag_qq(j1,j2,j3,j4,j5,j6,j7,a70h)
      implicit none
      include 'types.f'
      
************************************************************************
*     Returns the helicity amplitudes for the process                  *
*     0 ---> q(p1)+ph(p2)+ph(p3)+g(p4)+qbar(p5)+lb(p6)+l(p7)           *
*     both photons coming from quark line                              *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: i2(6),i3(6),i4(6),j,lh,h2,h3,h4,hq,h(7)
      integer:: j1,j2,j3,j4,j5,j6,j7
      integer:: hqc,lhc,h2c,h3c,h4c
      complex(dp):: a70h(2,2,2,2,2)
      complex(dp):: atemp,m(6),amp_qqggg
c-----possible permutation of (2,3,4)
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
      a70h(hq,h2,h3,h4,lh)=czip
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
         atemp=czip
         do j=1,6
            m(j)=amp_qqggg(j1,hq,i2(j),h(i2(j)),i3(j),h(i3(j)),
     &                     i4(j),h(i4(j)),j5,lh,j6,j7)
            atemp=atemp+m(j)
         enddo
         a70h(hq,h2,h3,h4,lh)=atemp
c--------obtain hq=1 from complex conjugation
         hqc=mod(hq+2,2)+1
         h2c=mod(h2+2,2)+1
         h3c=mod(h3+2,2)+1
         h4c=mod(h4+2,2)+1
         lhc=mod(lh+2,2)+1
         a70h(hqc,h2c,h3c,h4c,lhc)=(-one)**j4*
     &                             conjg(a70h(hq,h2,h3,h4,lh))
c--------
      enddo
      enddo
      enddo
      enddo
      enddo
c-----done     
      return
      end

      subroutine xzqqaag_ll(j1,j2,j3,j4,j5,j6,j7,a70h)
      implicit none
      include 'types.f'
************************************************************************
*     Returns the helicity amplitudes for the process                  *
*     0 ---> q(p1)+ph(p2)+ph(p3)+g(p4)+qbar(p5)+lb(p6)+l(p7)           *
*     both photons coming from lepton line                             *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: i2(2),i3(2),i4(2),j,lh,h2,h3,h4,hq,h(7)
      integer:: j1,j2,j3,j4,j5,j6,j7
      integer:: hqc,lhc,h2c,h3c,h4c
      complex(dp):: a70hx(2,2,2,2,2),a70h(2,2,2,2,2)
      complex(dp):: amp_qqagg_ql
      complex(dp):: atemp,m(2)
      integer:: ii
c-----permutation of (3,4)
      i2(1)=j2
      i3(1)=j3
      i4(1)=j4
      i2(2)=j2
      i3(2)=j4
      i4(2)=j3
c-----initialize helicity amplitudes to zero
c-----a70h(hq,hg2,hg3,hg4,hl)
      do hq=1,2
      do lh=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      a70hx(hq,h2,h3,h4,lh)=czip
      a70h(hq,h2,h3,h4,lh)=czip
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
         atemp=czip
c--------symmetrize over 3 and 4, 2 is coming from lepton line
         do j=1,2
            m(j)=amp_qqagg_ql(j1,hq,i2(j),h(i2(j)),i3(j),h(i3(j)),
     &                        i4(j),h(i4(j)),j5,lh,j6,j7)
            atemp=atemp+m(j)
         enddo
         a70hx(hq,h2,h3,h4,lh)=atemp
c--------obtain hq=1 from complex conjugation
         hqc=mod(hq+2,2)+1
         h2c=mod(h2+2,2)+1
         h3c=mod(h3+2,2)+1
         h4c=mod(h4+2,2)+1
         lhc=mod(lh+2,2)+1
         a70hx(hqc,h2c,h3c,h4c,lhc)=(-one)**j2*
     &                              conjg(a70hx(hq,h2,h3,h4,lh))
c-----
      enddo
      enddo
      enddo
      enddo
      enddo
c-----flip the helicity of 2(ph) and 4(glu)
c-----as well as 1(q) and 7(l)
c-----to compensate the momentum rotation from ql to ll
      do hq=1,2
      do lh=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      a70h(lh,h4,h3,h2,hq)=a70hx(hq,h2,h3,h4,lh)
      enddo
      enddo
      enddo
      enddo
      enddo
c-----done      
      return
      end
      
      subroutine xzqqaag_ql(j1,j2,j3,j4,j5,j6,j7,a70h3,a70h4)
      implicit none
      include 'types.f'
************************************************************************
*     Returns the helicity amplitudes for the process                  *
*     0 ---> q(p1)+ph(p2)+ph(p3)+g(p4)+qbar(p5)+lb(p6)+l(p7)           *
*     1 photon coming from quark line                                  *
*     1 photon coming from lepton line                                 *
************************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ipsgen.f'
      integer:: i2(4),i3(4),i4(4),j,lh,h2,h3,h4,hq,h(7)
      integer:: j1,j2,j3,j4,j5,j6,j7
      integer:: hqc,lhc,h2c,h3c,h4c
      complex(dp):: a70h3(2,2,2,2,2)
      complex(dp):: a70h4(2,2,2,2,2)
      complex(dp):: amp_qqagg_ql
      complex(dp):: atemp3,m3(2)
      complex(dp):: atemp4,m4(2)
      integer:: ii
c-----permutation over (2,3,4) for but 4 cant be in i2
      i2(1)=j2
      i3(1)=j3
      i4(1)=j4
      i2(2)=j2
      i3(2)=j4
      i4(2)=j3
      i2(3)=j3
      i3(3)=j2
      i4(3)=j4
      i2(4)=j3
      i3(4)=j4
      i4(4)=j2
c-----initialize helicity amplitudes to zero
c-----a70h(hq,hg2,hg3,hg4,hl)
      do hq=1,2
      do lh=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      a70h3(hq,h2,h3,h4,lh)=czip
      a70h4(hq,h2,h3,h4,lh)=czip
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
         atemp3=czip
         atemp4=czip
c--------symmetrize over 3 and 4, 2 is coming from lepton line
         do j=1,2
            if ((ipsgen == 2) .or. (ipsgen == 3)) then
            m3(j)=amp_qqagg_ql(j1,hq,i2(j),h(i2(j)),i3(j),h(i3(j)),
     &                         i4(j),h(i4(j)),j5,lh,j6,j7)
            atemp3=atemp3+m3(j)
            endif
         enddo
         do j=3,4
            if ((ipsgen == 3) .or. (ipsgen == 4)) then
            m4(j-2)=amp_qqagg_ql(j1,hq,i2(j),h(i2(j)),i3(j),h(i3(j)),
     &                         i4(j),h(i4(j)),j5,lh,j6,j7)
            atemp4=atemp4+m4(j-2)
            endif
         enddo
         a70h3(hq,h2,h3,h4,lh)=atemp3
         a70h4(hq,h2,h3,h4,lh)=atemp4
c--------obtain hq=1 from complex conjugation
         hqc=mod(hq+2,2)+1
         h2c=mod(h2+2,2)+1
         h3c=mod(h3+2,2)+1
         h4c=mod(h4+2,2)+1
         lhc=mod(lh+2,2)+1
         a70h3(hqc,h2c,h3c,h4c,lhc)=(-one)**j4*
     &                              conjg(a70h3(hq,h2,h3,h4,lh))
         a70h4(hqc,h2c,h3c,h4c,lhc)=(-one)**j4*
     &                              conjg(a70h4(hq,h2,h3,h4,lh))
c-----
      enddo
      enddo
      enddo
      enddo
      enddo
c-----done
      return
      end
 
