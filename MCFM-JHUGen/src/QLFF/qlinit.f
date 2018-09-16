      subroutine qlinit
      implicit none
      include 'mpicommon.f'

      if (rank == 0) then
      write(*,*) '===================================================='
      write(*,*) '  This is QCDLoop - version 1.96                    '
      write(*,*) '  Authors: Keith Ellis and Giulia Zanderighi        '
      write(*,*) '  (ellis@fnal.gov, g.zanderighi1@physics.ox.ac.uk)  '
      write(*,*) '  For details see FERMILAB-PUB-07-633-T,OUTP-07/16P '
      write(*,*) '  arXiv:0712.1851 [hep-ph], published in            '
      write(*,*) '  JHEP 0802:002,2008.                               '
      write(*,*) '===================================================='
      endif

      call ffini

      end
