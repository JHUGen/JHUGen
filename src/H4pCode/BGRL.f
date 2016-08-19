c--- Combinations of bubble integrals, 0704.3914v3 Eq. (3.21) 
c--- The hatted versions of these functions are further defined
c--- by Eqs. (3.25), (3.26) and (3.27)
      double complex function BGRL1(s,t)
      implicit none
      double complex lnrat
      double precision s,t
      BGRL1=lnrat(-s,-t)/dcmplx(s-t)
      return
      end

      double complex function BGRL2(s,t)
      implicit none
      double complex lnrat
      double precision s,t
      BGRL2=lnrat(-s,-t)/dcmplx(s-t)**2
      return
      end

      double complex function BGRL2hat(s,t)
      implicit none
      double complex lnrat
      double precision s,t
      BGRL2hat=lnrat(-s,-t)/dcmplx(s-t)**2
     . -dcmplx(0.5d0*(s+t)/((s-t)*s*t))
      return
      end

      double complex function BGRL3(s,t)
      implicit none
      double complex lnrat
      double precision s,t
      BGRL3=lnrat(-s,-t)/dcmplx(s-t)**3
      return
      end

      double complex function BGRL3hat(s,t)
      implicit none
      double complex lnrat
      double precision s,t
      BGRL3hat=lnrat(-s,-t)/dcmplx(s-t)**3
     . -dcmplx(0.5d0*(s+t)/((s-t)**2*s*t))
      return
      end
