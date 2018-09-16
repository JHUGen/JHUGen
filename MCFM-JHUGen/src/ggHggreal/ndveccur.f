      subroutine ndveccur(i1,i2,n,p,vcm)
      implicit none
      include 'types.f'

C-----Author R.K. Ellis
C-----October 2005
C     Calculate the basic vector current current
C         <i1-|Gamma(mu)|i2->*n(mu)
C     and <i2-|Gamma(mu)|i1->*n(mu)
C
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      integer::i1,i2
      real(dp)::p(mxpart,4),n(4),rtbp,rtpp
      complex(dp)::fac,c23b,c23p,fb,fp,vc(4),vcm(mxpart,mxpart)
         if (p(i1,4) > zip) then
C-----positive energy case
            rtbp=sqrt(p(i1,4)+p(i1,1))
            c23b=cplx2(p(i1,3),-p(i1,2))
            fb=cone
         else
C-----negative energy case
            rtbp=sqrt(-p(i1,4)-p(i1,1))
            c23b=cplx2(-p(i1,3),p(i1,2))
            fb=im
         endif

         if (p(i2,4) > zip) then
C-----positive energy case
            rtpp=sqrt(p(i2,4)+p(i2,1))
            c23p=cplx2(p(i2,3),-p(i2,2))
            fp=cone
         else
C-----negative energy case
            rtpp=sqrt(-p(i2,4)-p(i2,1))
            c23p=cplx2(-p(i2,3),p(i2,2))
            fp=im
         endif

      fac=fb*fp
      vc(3)=+fac
     & *((rtbp/rtpp)*conjg(c23p)+(rtpp/rtbp)*c23b)
      vc(2)=fac
     & *(-im*(rtbp/rtpp)*conjg(c23p)+im*(rtpp/rtbp)*c23b)
      vc(1)=fac
     & *(cplx2(rtbp*rtpp,zip)-c23b*conjg(c23p)/(rtbp*rtpp))
      vc(4)=fac
     & *(cplx2(rtbp*rtpp,zip)+c23b*conjg(c23p)/(rtbp*rtpp))

      vcm(i1,i2)=n(4)*vc(4)-n(1)*vc(1)-n(2)*vc(2)-n(3)*vc(3)

      vc(3)=fac
     & *((rtbp/rtpp)*c23p+(rtpp/rtbp)*conjg(c23b))
      vc(2)=fac
     & *(im*((rtbp/rtpp))*c23p-im*(rtpp/rtbp)*conjg(c23b))
      vc(1)=fac
     & *(cplx2(rtbp*rtpp,zip)-conjg(c23b)*c23p/(rtbp*rtpp))
      vc(4)=fac
     & *(cplx2(rtbp*rtpp,zip)+conjg(c23b)*c23p/(rtbp*rtpp))

      vcm(i2,i1)=n(4)*vc(4)-n(1)*vc(1)-n(2)*vc(2)-n(3)*vc(3)
c      write(6,*) 'vcm(i1,i2)',vcm(i1,i2)
c      write(6,*) 'vcm(i2,i1)',vcm(i2,i1)
      return
      end

