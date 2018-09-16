      program mcfm
      implicit none
      include 'types.f'
      include 'mpif.h'
      include "omp_lib.h"
      include 'mpicommon.f'
      integer ierr,mylen,support
      character (len=mpi_max_processor_name)::procname
      real(dp):: r0,er0,t0,t1,dumran

      call srand(5699)
*
      call mpi_init_thread(mpi_thread_multiple,support,ierr)
      call mpi_comm_rank(mpi_comm_world,rank,ierr)
      call mpi_comm_size(mpi_comm_world,size,ierr)
      call mpi_get_processor_name(procname,mylen,ierr)

      t0=mpi_wtime()
      call mcfmsub(r0,er0)
      t1=mpi_wtime()
      if (rank.eq.0) write(*,*) 'time to finish',t1-t0,rand()
      call mpi_finalize(ierr)
      end
