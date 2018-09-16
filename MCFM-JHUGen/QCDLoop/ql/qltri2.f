      subroutine qltri2(p1sq,p2sq,musq,Ival)
      implicit none
      include 'qlconstants.f'
      double precision p1sq,p2sq,musq,ron2
      double complex qllnrat,Ival(-2:0),wlog1,wlog2

      wlog1=qllnrat(musq,-p1sq)
      wlog2=qllnrat(musq,-p2sq)
      ron2=0.5d0*(p2sq-p1sq)/p1sq

      Ival(-2)=czip
      if (abs(ron2) .lt. 1d-6) then
      Ival(-1)=-dcmplx(1d0/p1sq)*dcmplx(1d0-ron2)
      Ival(0)=Ival(-1)*wlog1+dcmplx(ron2/p1sq)
      else
      Ival(-1)=(wlog1-wlog2)/dcmplx(p1sq-p2sq)
      Ival(0)=chalf*Ival(-1)*(wlog1+wlog2)
      endif
      return
      end
