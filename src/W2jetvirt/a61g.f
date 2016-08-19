      double complex function a61g(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
C---hep-ph/9708239, Eqn 2.13
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6
      character*9 st
      double complex a6g,a6sg,a6fg,a6tg
      a61g=
     . +a6g(st,j1,j2,j3,j4,j5,j6,za,zb)
     . -a6g(st,j1,j4,j3,j2,j5,j6,za,zb)/xnsq
     . +(-dfloat(nf)*(a6sg(st,j1,j2,j3,j4,j5,j6,za,zb)
     .               +a6fg(st,j1,j2,j3,j4,j5,j6,za,zb))
     . +a6tg(st,j1,j2,j3,j4,j5,j6,za,zb))/xn
      return
      end

      double complex function a63g(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
C---hep-ph/9708239, Eqn 2.13
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6
      character*9 st,stamp(6)
      double complex a6g
      
c      write(*,*) 'Called a63g with arguments'
c      write(*,*) 'st',st
c      write(*,*) 'j1,j4,j2,j3,j5,j6',j1,j4,j2,j3,j5,j6
c      pause
c--- although a63 is labelled as (q1,qb4,g2,g3), the individual
c---  amplitudes are labelled by (q ,g ,g ,qb) [for 1 and 2]
c---                         and (q, g ,qb, g) [for 3 and 4]
c---                         and (q,qb , g, g) [for 5 and 6]

      if     (st .eq. 'q+qb-g+g+') then
        stamp(1)='q+g+g+qb-'
        stamp(2)='q+g+g+qb-'
        stamp(5)='q+qb-g+g+'
        stamp(6)='q+qb-g+g+'
      elseif (st .eq. 'q+qb-g+g-') then
        stamp(1)='q+g+g-qb-'
        stamp(2)='q+g-g+qb-'
        stamp(5)='q+qb-g+g-'
        stamp(6)='q+qb-g-g+'
      elseif (st .eq. 'q+qb-g-g+') then
        stamp(1)='q+g-g+qb-'
        stamp(2)='q+g+g-qb-'
        stamp(5)='q+qb-g-g+'
        stamp(6)='q+qb-g+g-'
      endif

      a63g=a6g(stamp(1),j1,j2,j3,j4,j5,j6,za,zb)
     .    +a6g(stamp(2),j1,j3,j2,j4,j5,j6,za,zb)
     .    +a6g(stamp(5),j1,j4,j2,j3,j5,j6,za,zb)
     .    +a6g(stamp(6),j1,j4,j3,j2,j5,j6,za,zb)

      if     (st .eq. 'q+qb-g+g+') then
        stamp(3)='q+g+qb-g+'
        stamp(4)='q+g+qb-g+'
        a63g=a63g
     .    +a6g(stamp(3),j1,j2,j4,j3,j5,j6,za,zb)
     .    +a6g(stamp(4),j1,j3,j4,j2,j5,j6,za,zb)
      elseif (st .eq. 'q+qb-g+g-') then
        stamp(3)='q+g+qb-g-'
        stamp(4)='q+g+qb-g-'
        a63g=a63g
     .    +a6g(stamp(3),j1,j2,j4,j3,j5,j6,za,zb)
     .    +a6g(stamp(4),j4,j3,j1,j2,j6,j5,zb,za)
      elseif (st .eq. 'q+qb-g-g+') then
        stamp(3)='q+g+qb-g-'
        stamp(4)='q+g+qb-g-'
        a63g=a63g
     .    +a6g(stamp(3),j4,j2,j1,j3,j6,j5,zb,za)
     .    +a6g(stamp(4),j1,j3,j4,j2,j5,j6,za,zb)
      else
        write(*,*) 'Unimplemented st in a63g'
        stop
      endif

      return
      end

      double complex function a64vg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
C---hep-ph/9708239, Eqn 2.13
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6
      character*9 st
      double complex a6vsg,a6vfg
      a64vg=-a6vsg(st,j1,j4,j2,j3,j5,j6,za,zb)
     .      -a6vfg(st,j1,j4,j2,j3,j5,j6,za,zb)
      return
      end

      double complex function a64axg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
C---hep-ph/9708239, Eqn 2.13
      include 'constants.f'
      include 'zprods_decl.f'
      integer j1,j2,j3,j4,j5,j6
      character*9 st
      double complex a6axg
      a64axg=a6axg(st,j1,j4,j2,j3,j5,j6,za,zb)
      return
      end

      double complex function a65axg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
C---hep-ph/9708239, Eqn 2.13
      include 'constants.f'
      include 'zprods_decl.f'
      character*9 st
      integer j1,j2,j3,j4,j5,j6
      double complex a6axslg
      a65axg=a6axslg(st,j1,j4,j2,j3,j5,j6,za,zb)
      return
      end

      double complex function a6sg(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      character*9 st
      integer j1,j2,j3,j4,j5,j6
c--- modified 4/18/00, according to BDK hep-ph/9708239, eq. (8.2)
      a6sg=czip
      if (st .eq. 'q+g+g+qb-') 
     .  a6sg=one/three/za(j2,j3)**2/s(j5,j6)*(
     .        -(zb(j6,j1)*za(j1,j3)+zb(j6,j2)*za(j2,j3))*zb(j3,j1)
     .         *za(j4,j5)/(s(j1,j2)+s(j1,j3)+s(j2,j3))
     .        +(za(j5,j4)*zb(j4,j3)+za(j5,j2)*zb(j2,j3))*za(j3,j4)
     .         *zb(j1,j6)/(s(j2,j3)+s(j2,j4)+s(j3,j4)))
      return
      end

      double complex function a6fg(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      character*9 st
      integer j1,j2,j3,j4,j5,j6
c--- checked 4/18/00, this piece vanishes identically, eq. (8.3)
      a6fg=czip
      return
      end

      double complex function a6tg(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      character*9 st
      integer j1,j2,j3,j4,j5,j6
      double complex a6sg
c--- modified 4/18/00, according to BDK hep-ph/9708239, eq. (8.3)
      a6tg=czip
      if (st .eq. 'q+g+g+qb-') 
     .  a6tg=one/20d0*s(j2,j3)/mt**2*a6sg(st,j1,j2,j3,j4,j5,j6,za,zb)
      return
      end

      double complex function a6vfg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      character*9 st
      integer j1,j2,j3,j4,j5,j6
      a6vfg=czip
      return
      end

      double complex function a6vsg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      character*9 st
      integer j1,j2,j3,j4,j5,j6
      a6vsg=czip
      return
      end

      double complex function a6axg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      character*9 st
      integer j1,j2,j3,j4,j5,j6
      a6axg=czip
      return
      end

      double complex function a6axslg(st,j1,j4,j2,j3,j5,j6,za,zb)
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      character*9 st
      integer j1,j2,j3,j4,j5,j6
      a6axslg=czip
      return
      end
