      subroutine ndveccur(i1,i2,n,p,vcm)
      implicit none
C-----Author R.K. Ellis
C-----October 2005
C     Calculate the basic vector current current 
C         <i1-|Gamma(mu)|i2->*n(mu)
C     and <i2-|Gamma(mu)|i1->*n(mu)
C      
      include 'constants.f'
      integer i1,i2
      double precision p(mxpart,4),n(4),rtbp,rtpp
      double complex fac,c23b,c23p,fb,fp,vc(4),vcm(mxpart,mxpart) 
         if (p(i1,4) .gt. 0d0) then
C-----positive energy case
            rtbp=dsqrt(p(i1,4)+p(i1,1))
            c23b=dcmplx(p(i1,3),-p(i1,2))
            fb=cone
         else
C-----negative energy case
            rtbp=dsqrt(-p(i1,4)-p(i1,1))
            c23b=dcmplx(-p(i1,3),p(i1,2))
            fb=im
         endif

         if (p(i2,4) .gt. 0d0) then
C-----positive energy case
            rtpp=dsqrt(p(i2,4)+p(i2,1))
            c23p=dcmplx(p(i2,3),-p(i2,2))
            fp=cone
         else
C-----negative energy case
            rtpp=dsqrt(-p(i2,4)-p(i2,1))
            c23p=dcmplx(-p(i2,3),p(i2,2))
            fp=im
         endif

      fac=fb*fp
      vc(3)=+fac
     . *(dcmplx(rtbp/rtpp)*dconjg(c23p)+dcmplx(rtpp/rtbp)*c23b)
      vc(2)=fac
     . *(-im*dcmplx(rtbp/rtpp)*dconjg(c23p)+im*dcmplx(rtpp/rtbp)*c23b)
      vc(1)=fac
     . *(dcmplx(rtbp*rtpp)-c23b*dconjg(c23p)/dcmplx(rtbp*rtpp))
      vc(4)=fac
     . *(dcmplx(rtbp*rtpp)+c23b*dconjg(c23p)/dcmplx(rtbp*rtpp))

      vcm(i1,i2)=n(4)*vc(4)-n(1)*vc(1)-n(2)*vc(2)-n(3)*vc(3)

      vc(3)=fac
     . *(dcmplx(rtbp/rtpp)*c23p+dcmplx(rtpp/rtbp)*dconjg(c23b))
      vc(2)=fac
     . *(im*(dcmplx(rtbp/rtpp))*c23p-im*dcmplx(rtpp/rtbp)*dconjg(c23b))
      vc(1)=fac
     . *(dcmplx(rtbp*rtpp)-dconjg(c23b)*c23p/dcmplx(rtbp*rtpp))
      vc(4)=fac
     . *(dcmplx(rtbp*rtpp)+dconjg(c23b)*c23p/dcmplx(rtbp*rtpp))

      vcm(i2,i1)=n(4)*vc(4)-n(1)*vc(1)-n(2)*vc(2)-n(3)*vc(3)
c      write(6,*) 'vcm(i1,i2)',vcm(i1,i2)
c      write(6,*) 'vcm(i2,i1)',vcm(i2,i1)
      return
      end 

