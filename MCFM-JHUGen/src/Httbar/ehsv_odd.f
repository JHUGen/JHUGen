c--- These are the functions similar to equations (A.8)-(A.14)
c--- but which describe the case of a pseudoscalar (CP odd) Higgs
c--- These results are adapted from Spira et al., hep-ph/9504378.

      double complex function ehsva4_odd(s,t,u)
C     ehsv:EqnA.8
      implicit none
      double precision s,t,u
      double complex ehsvb4_odd
      ehsva4_odd=ehsvb4_odd(s,t,u)+ehsvb4_odd(u,s,t)+ehsvb4_odd(t,u,s)
      return 
      end

      double complex function ehsva2_odd(s,t,u)
C     ehsv:EqnA.9
      implicit none
      double precision s,t,u
      double complex ehsvb2_odd
      ehsva2_odd=ehsvb2_odd(s,t,u)+ehsvb2_odd(s,u,t)
      return 
      end

      double complex function ehsvb4_odd(s,t,u)
      implicit none
C     ehsv:EqnA.10
      include 'masses.f'
      double precision hmass2,s,t,u
      double complex w2,w3
      hmass2=s+t+u
      ehsvb4_odd=mbsq/hmass2*(w2(hmass2)-w2(s)-w3(s,t,u,hmass2))/6d0
      return 
      end

      double complex function ehsvb2_odd(s,t,u)
C     ehsv:EqnA.11
      implicit none
      include 'masses.f'
      double precision hmass2,s,t,u
      double complex w2,w3
      hmass2=s+t+u
      ehsvb2_odd=mbsq/hmass2**2*(
     .                       -s*w3(s,t,u,hmass2)+s/2d0*w3(t,s,u,hmass2)
     .                       +2d0*s*(0.75d0-u/(s+u))*w2(hmass2)
     .                       -s/2d0*w2(s)
     .                       -s*(1d0-2d0*u/(s+u))*w2(t))/6d0
      return 
      end

      double complex function ehsva5_odd(s,t,u)
C     ehsv:EqnA.14
      implicit none
      include 'masses.f'
      double precision hmass2,s,t,u
      double complex w2
      hmass2=s+t+u
      ehsva5_odd=mbsq/hmass2*2d0/3d0*(w2(s)-w2(hmass2))
      return 
      end

