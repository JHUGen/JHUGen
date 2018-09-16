      subroutine qqb_zgam_v(p,msqv)
      implicit none
      include 'types.f'
      
C----Author R.K.Ellis October 2002
C====Virtual corrections to
c     q(-p1)+qbar(-p2)-->e^-(p3)+e-(p4)+gamma(p5)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'scheme.f'
      include 'zprods_decl.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer:: j,h12,h34,h5
      real(dp):: msqv(-nf:nf,-nf:nf),p(mxpart,4),
     & qbq(2),qqb(2)
      real(dp):: fac
      complex(dp):: qbqlo(2,2,2,2),qqblo(2,2,2,2),
     &              qbqnlo(2,2,2,2),qqbnlo(2,2,2,2)
      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      scheme='dred'

      call spinoru(5,p,za,zb)
      call zvirtamps(1,2,3,4,5,za,zb,qbqnlo,qbqlo)
      call zvirtamps(2,1,3,4,5,za,zb,qqbnlo,qqblo)

      fac=cf*ason2pi*aveqq*8._dp*esq**3*xn
      do j=1,2
      qbq(j)=0._dp
      qqb(j)=0._dp
      do h12=1,2
      do h34=1,2
      do h5= 1,2
      qbq(j)=qbq(j)
     & +2._dp*Dble(conjg(qbqlo(j,h12,h34,h5))*qbqnlo(j,h12,h34,h5))
      qqb(j)=qqb(j)
     & +2._dp*Dble(conjg(qqblo(j,h12,h34,h5))*qqbnlo(j,h12,h34,h5))
      enddo
      enddo
      enddo
      qbq(j)=fac*qbq(j)
      qqb(j)=fac*qqb(j)
      enddo

      msqv(:,:)=0._dp
      do j=-nf,nf
        if     (j > 0) then
          msqv(j,-j)=qqb(jj(j))
        elseif (j < 0) then
          msqv(j,-j)=qbq(-jj(j))
        endif
      enddo

      return
      end


      subroutine zvirtamps(p1,p2,p3,p4,p5,za,zb,qbqnlo,qbqlo)
      implicit none
      include 'types.f'
      
C  Amplitudes for virtual corrections taken from Dixon,Kunszt,Signer
C  hep-ph/9803250 Eqns 4.6 and 4.7
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      integer:: p1,p2,p3,p4,p5,j
      real(dp):: s34
      complex(dp):: qbqlo(2,2,2,2),qbqnlo(2,2,2,2),fagamma,fbgamma,
     & prp34,vpole,vpl

      call zamps(p1,p2,p3,p4,p5,za,zb,qbqlo)

c--   calculate propagator factors
      s34=s(p3,p4)      
      prp34=s34/cplx2((s34-zmass**2),zmass*zwidth)
      vpl=vpole(s(p1,p2))

      do j=1,2

      qbqnlo(j,1,1,2)=
     &  +(Q(j)*(Q(j)*q1+L(j)*l1*prp34)
     & *(fagamma(p1,p2,p3,p4,p5,za,zb)
     &  +fbgamma(p1,p2,p3,p4,p5,za,zb)))
     &          +vpl*qbqlo(j,1,1,2)
      qbqnlo(j,1,2,2)=
     &  +(Q(j)*(Q(j)*q1+L(j)*r1*prp34)
     & *(fagamma(p1,p2,p4,p3,p5,za,zb)
     &  +fbgamma(p1,p2,p4,p3,p5,za,zb)))
     &          +vpl*qbqlo(j,1,2,2)
      qbqnlo(j,1,1,1)=
     &  +(Q(j)*(Q(j)*q1+L(j)*l1*prp34)
     .*(fagamma(p2,p1,p4,p3,p5,zb,za)
     & +fbgamma(p2,p1,p4,p3,p5,zb,za)))
     &          +vpl*qbqlo(j,1,1,1)
      qbqnlo(j,1,2,1)=
     &  +(Q(j)*(Q(j)*q1+L(j)*r1*prp34)
     .*(fagamma(p2,p1,p3,p4,p5,zb,za)
     & +fbgamma(p2,p1,p3,p4,p5,zb,za)))
     &          +vpl*qbqlo(j,1,2,1)

      qbqnlo(j,2,2,1)=      
     &  -(Q(j)*(Q(j)*q1+R(j)*r1*prp34)
     & *(fagamma(p1,p2,p3,p4,p5,zb,za)
     &  +fbgamma(p1,p2,p3,p4,p5,zb,za)))
     &          +vpl*qbqlo(j,2,2,1)
      qbqnlo(j,2,1,1)=      
     &  -(Q(j)*(Q(j)*q1+R(j)*l1*prp34)
     & *(fagamma(p1,p2,p4,p3,p5,zb,za)
     &  +fbgamma(p1,p2,p4,p3,p5,zb,za)))
     &          +vpl*qbqlo(j,2,1,1)
      qbqnlo(j,2,2,2)=      
     &  -(Q(j)*(Q(j)*q1+R(j)*r1*prp34)
     & *(fagamma(p2,p1,p4,p3,p5,za,zb)
     &  +fbgamma(p2,p1,p4,p3,p5,za,zb)))
     &          +vpl*qbqlo(j,2,2,2)
      qbqnlo(j,2,1,2)=      
     &  -(Q(j)*(Q(j)*q1+R(j)*l1*prp34)
     & *(fagamma(p2,p1,p3,p4,p5,za,zb)
     &  +fbgamma(p2,p1,p3,p4,p5,za,zb)))
     &          +vpl*qbqlo(j,2,1,2)
      enddo
      return
      end
