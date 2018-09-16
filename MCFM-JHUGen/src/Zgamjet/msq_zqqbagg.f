      subroutine zagg_a70h(j1,j2,j3,j4,j5,j6,j7,a70h1,a70h3)
      implicit none
      include 'types.f'
*******************************************************************
* 0 -> q(-p1) + qb(-p5) + a(p2) + g(p3) + g(p4) + lb(p6) + l(p7)
* return helicity amplitudes for each channel
********************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer:: j1,j2,j3,j4,j5,j6,j7
      integer:: hq,h2,h3,h4,lh,ic
      complex(dp):: a70h1(2,2,2,2,2,2),a70h3(2,2,2,2,2,2)
c-----initialize
      do hq=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do lh=1,2
      do ic=1,2
      a70h1(ic,hq,h2,h3,h4,lh)=czip
      a70h3(ic,hq,h2,h3,h4,lh)=czip
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
c----compute amplitude
      call xzqqagg_qq(j1,j2,j3,j4,j5,j6,j7,a70h1)
      call xzqqagg_ql(j1,j2,j3,j4,j5,j6,j7,a70h3)
c-----done
      return 
      end

      subroutine zagg_m70sq(qi,a70h1,a70h3,msq)
      implicit none
      include 'types.f'
********************************************************************
* 0 -> q(-p1) + qb(-p5) + a(p2) + g(p3) + g(p4) + lb(p6) + l(p7)
* return matrix element squared, for initial flavor qi
* given the helicity amplitudes from each channel
********************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'new_pspace.f'
      include 'ipsgen.f'
      real(dp):: msq,t,qq,gg,ee,m70hsq(32),CF1,CF2
      real(dp):: m70hsqAA(32),m70hsqBB(32),m70hsqAB(32)
      integer:: i,j,k,ihel,ic
      integer:: h2,h3,h4,qi
      complex(dp):: a70h1(2,2,2,2,2,2),a70h3(2,2,2,2,2,2)
      complex(dp):: m70hA(2,32),m70hB(2,32)
      complex(dp):: propQQ,propLL,propQL3,propQL4
c-----
      ee=sqrt(esq)
      gg=sqrt(gsq)
      qq=abs(q1)
c-----
      propQQ =s(6,7)/cplx2(s(6,7)-zmass**2,zwidth*zmass)
      propQL3=t(2,6,7)/cplx2(t(2,6,7)-zmass**2,zwidth*zmass)
c-----      
      ihel=1
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do ic=1,2
         m70hA(ic,ihel)=four*ee**3*gg**2*(
     &      (q1*Q(qi)+l(qi)*l1*propQQ)*Q(qi)*a70h1(ic,1,h2,h3,h4,1))
         m70hB(ic,ihel)=four*ee**3*gg**2*(
     & +qq*(q1*Q(qi)+l(qi)*l1*propQL3)*(-one)*a70h3(ic,1,h2,h3,h4,1))
      enddo
      ihel=ihel+1
      enddo
      enddo
      enddo
c-----
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do ic=1,2
         m70hA(ic,ihel)=four*ee**3*gg**2*(
     &      (q1*Q(qi)+r(qi)*l1*propQQ)*Q(qi)*a70h1(ic,2,h2,h3,h4,1))
         m70hB(ic,ihel)=four*ee**3*gg**2*(
     & -qq*(q1*Q(qi)+r(qi)*l1*propQL3)*(-one)*a70h3(ic,2,h2,h3,h4,1))
      enddo
      ihel=ihel+1
      enddo
      enddo
      enddo
c-----
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do ic=1,2
         m70hA(ic,ihel)=four*ee**3*gg**2*(
     &      (q1*Q(qi)+r(qi)*r1*propQQ)*Q(qi)*a70h1(ic,2,h2,h3,h4,2))
         m70hB(ic,ihel)=four*ee**3*gg**2*(
     & +qq*(q1*Q(qi)+r(qi)*r1*propQL3)*(-one)*a70h3(ic,2,h2,h3,h4,2))
      enddo
      ihel=ihel+1
      enddo
      enddo
      enddo
c-----
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do ic=1,2
         m70hA(ic,ihel)=four*ee**3*gg**2*(
     &     (q1*Q(qi)+l(qi)*r1*propQQ)*Q(qi)*a70h1(ic,1,h2,h3,h4,2))
         m70hB(ic,ihel)=four*ee**3*gg**2*(
     & -qq*(q1*Q(qi)+l(qi)*r1*propQL3)*(-one)*a70h3(ic,1,h2,h3,h4,2))
      enddo
      ihel=ihel+1
      enddo
      enddo
      enddo
c-----overcount ihel
      ihel=ihel-1
c-----square them up
      CF1=64/three
      CF2=-8/three
      do i=1,32
         m70hsqAA(i)=CF1*(abs(m70hA(1,i))**2+abs(m70hA(2,i))**2)
     &              +CF2*two*real(m70hA(1,i)*conjg(m70hA(2,i)))
         m70hsqBB(i)=CF1*(abs(m70hB(1,i))**2+abs(m70hB(2,i))**2)
     &              +CF2*two*real(m70hB(1,i)*conjg(m70hB(2,i)))
         m70hsqAB(i)=CF1*( 2._dp*real(conjg(m70hA(1,i))*m70hB(1,i))
     &                    +2._dp*real(conjg(m70hA(2,i))*m70hB(2,i)) )
     &              +CF2*( two*real(conjg(m70hA(1,i))*m70hB(2,i))
     &                    +two*real(conjg(m70hB(1,i))*m70hA(2,i)) )
      enddo
      do i=1,32
c         m70hsq(i)=m70hsqAA(i)+m70hsqBB(i)+m70hsqAB(i)
         if     (ipsgen == 1) then
             m70hsq(i)=m70hsqAA(i)
         elseif (ipsgen == 2) then
             m70hsq(i)=m70hsqBB(i)+m70hsqAB(i)
         else
            write(6,*) 'Parameter ipsgen should be 1 or 2'
            write(6,*) 'ipsgen = ',ipsgen
            stop
         endif
c         if (new_pspace) then
c             m70hsq(i)=m70hsqBB(i)+m70hsqAB(i)
c         else
c             m70hsq(i)=m70hsqAA(i)
c         endif
      enddo
c-----and sum them up
      msq=zip
      do i=1,32
         msq=msq+m70hsq(i)
      enddo
c-----include color factors
c-----no averaging (2 from removing one Ta)
c-----no identical factor
      msq=msq/2._dp
c-----done
      return
      end

