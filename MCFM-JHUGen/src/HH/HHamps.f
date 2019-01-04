       subroutine HHamps(p1,p2,p3,p4,gauge)
C--    Formula taken from Glover and van der Bij
C--    Nucl. Phys. B309 (1988) 202
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'scale.f'
      include 'first.f'
      double precision p1(4),p2(4),p3(4),p4(4),ss,tt,uu,
     & p3sq,p4sq,mQsq,mhsq
      double complex triangle(1,2),box(1,2),gauge(1,2),
     & D123,D213,D132,C12,C23,C13,C34,qlI4,qlI3

c--- initialize QCDLoop, if necessary
      if (first) then
        call qlinit
        first=.false.
      endif

C----Signs for momenta chosen such that p1+p2+p3+p4=0
      ss=(p1(4)+p2(4))**2
     &  -(p1(1)+p2(1))**2-(p1(2)+p2(2))**2-(p1(3)+p2(3))**2
      tt=(p1(4)+p3(4))**2
     &  -(p1(1)+p3(1))**2-(p1(2)+p3(2))**2-(p1(3)+p3(3))**2
      uu=(p2(4)+p3(4))**2
     &  -(p2(1)+p3(1))**2-(p2(2)+p3(2))**2-(p2(3)+p3(3))**2
      p3sq=p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2
      p4sq=p4(4)**2-p4(1)**2-p4(2)**2-p4(3)**2
      mQsq=mt**2
      mhsq=hmass**2

      C12=qlI3(0d0,0d0 ,ss,mQsq,mQsq,mQsq,musq,0)
      C23=qlI3(0d0,p3sq,uu,mQsq,mQsq,mQsq,musq,0)
      C13=qlI3(0d0,p3sq,tt,mQsq,mQsq,mQsq,musq,0)
      C34=qlI3(p3sq,p4sq,ss,mQsq,mQsq,mQsq,musq,0)
      D123=qlI4(0d0,0d0,p3sq,p4sq,ss,uu,mQsq,mQsq,mQsq,mQsq,musq,0)
      D213=qlI4(0d0,0d0,p3sq,p4sq,ss,tt,mQsq,mQsq,mQsq,mQsq,musq,0)
      D132=qlI4(0d0,p3sq,0d0,p4sq,tt,uu,mQsq,mQsq,mQsq,mQsq,musq,0)

      triangle(1,1)=-12d0*mhsq*mQsq/dcmplx(ss-mhsq,hmass*hwidth)
     & *(2d0+(4d0*mQsq-ss)*C12)
      triangle(1,2)=czip

      box(1,1)=-4d0*mQsq*(
     & mQsq*(8d0*mQsq-ss-2d0*mhsq)*(D123+D213+D132)
     & +(uu*tt-mhsq**2)/ss*(4d0*mQsq-mhsq)*D132+2d0+4d0*mQsq*C12
     & +2d0/ss*(mhsq-4d0*mQsq)*((tt-mhsq)*C13+(uu-mhsq)*C23)) 

      box(1,2)=-2d0*mQsq*(
     & 2d0*(8d0*mQsq+ss-2d0*mhsq)*(mQsq*(D123+D213+D132)-C34)
     & -2d0*(ss*C12+(tt-mhsq)*C13+(uu-mhsq)*C23)
     & +1d0/(uu*tt-mhsq**2)*(ss*uu*(8d0*uu*mQsq-uu**2-mhsq**2)*D123
     & +ss*tt*(8d0*tt*mQsq-tt**2-mhsq**2)*D213
     & +(8d0*mQsq+ss-2d0*mhsq)*(ss*(ss-2d0*mhsq)*C12
     & +ss*(ss-4d0*mhsq)*C34
     & +2d0*tt*(mhsq-tt)*C13+2d0*uu*(mhsq-uu)*C23)))

      gauge(1,1)=triangle(1,1)+box(1,1)
      gauge(1,2)=triangle(1,2)+box(1,2)
      return
      end
