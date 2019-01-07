      real(kind=qp) zip,zero,half,one,two,three,four,eight
      parameter(zip=0._qp,zero=0._qp,half=0.5_qp,one=1._qp,two=2._qp)
      parameter(three=3._qp,four=4._qp,eight=8._qp) 
      real(kind=qp) pi,pisq,pisqo6,pi3,pi6
      parameter(pi=3.14159265358979323846264338327948_qp,
     & pisq=pi*pi,pisqo6=pisq/6._qp,pi3=pisq/three,pi6=pisqo6)

      complex(kind=qp) im,impi,czip,chalf,cone,ctwo,c2ipi
      parameter(im=cmplx(zip,one,kind=qp),impi=cmplx(zip,pi,kind=qp),
     & czip=cmplx(zip,zip,kind=qp),chalf=cmplx(half,zip,kind=qp),
     & cone=cmplx(one,zip,kind=qp),ctwo=cmplx(two,zip,kind=qp), 
     & c2ipi=two*impi)

      complex(kind=qp) nan
      parameter (nan = (1.e123_qp, 1.e123_qp))


