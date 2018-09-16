      subroutine read_dm_params
      implicit none
      include 'types.f'
       
      include 'dm_params.f' 
      character*200 line 
      integer:: j
      open(unit=35,file='dm_parameters.DAT',status='old',err=440) 

      read(35,*) xmass  
      read(35,*) dm_lam
      read(35,*) effective_th
      read(35,*) yukawa_scal 
      read(35,*) line
!-----
      do j=1,5 
         read(35,*) dmL(j) 
      enddo
      read(35,*) line 
!-----
      do j=1,5 
         read(35,*) dmR(j) 
      enddo
      read(35,*) line 
      
      read(35,*) medmass
      read(35,*) medwidth
      read(35,*) g_dmx
      read(35,*) g_dmq
      close(35)

      write(6,*) '************** DM PARAMETERS ****************' 
      write(6,43) '* DM Mass ',xmass,'          *' 
!      write(6,43) '* DM coupling ',gx,'           *' 
      write(6,43) '* Lambda  ',dm_lam,'              *'      
      write(6,*) '* Effective theory : ',effective_th,'  * '
      write(6,44)  '* dm L (d,u,s,c,b),',dmL,' *'
      write(6,44)  '* dm R (d,u,s,c,b),',dmR,' *'
      if(effective_th.eqv..false.) then 
         write(6,*) 'Mediator Mass,',medmass
         write(6,*) 'Mediator Width,',medwidth
         write(6,*) 'g_x ',g_dmx
         write(6,*) 'g_q ',g_dmq
      endif
 !     write(6,43) '* Mediator Mass ',medmass,'        *' 
      write(6,*) '***************************************' 

c--- for scalar or pseudoscalar mediator, set scalar couplings
      if (  (dm_mediator == 'scalar')
     & .or. (dm_mediator == 'pseudo') ) then
        call set_scalar_coups()
      endif
      
c--- for axial-vector or pseudoscalar mediator, check couplings
      if (  (dm_mediator == 'axvect')
     & .or. (dm_mediator == 'pseudo') ) then
        call check_dmAxC()
      endif

 43   format(1x,a10,f8.3,a20) 

 44   format(1x,a20,5(f8.4,' '),a4)
      return 
 440  write(6,*) 'Could not open dm_parameters.DAT' 
      stop 
      
      return 
      end 




      subroutine check_dmAxC
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'dm_params.f' 
!----- check that dmL = -dmR for axial coupling 
      integer:: j 
      logical:: reset 

      reset=.false. 

      do j=1,nf
         if(abs(dmL(j)+dmR(j))>1.e-8_dp) then 
          dmR(j)=-dmL(j)
          reset=.true. 
          endif
      enddo

      
      if(reset) then 
         write(6,*) 'Found that dmL is not equal to -dm R as required'
         write(6,*) 'Fixing dmR = -dmL (input) '
         write(6,*) 'to ensure correct axial behaviour.'
         write(6,*) 'Check manual for details.'
      endif

      return 
      end 
      

      
      subroutine set_scalar_coups
      implicit none
      include 'types.f'
       
      include 'dm_params.f' 
      include 'masses.f' 
      include 'ewcouple.f' 
      include 'first.f' 
      integer:: j 
      real(dp):: inv_v
!----- inv v = 1/v = sqrt(gwsq/4m_W^2) 
      inv_v=sqrt(gwsq/(4._dp*wmass**2)) 

      if(yukawa_scal) then 
      if(effective_th) then 
         dmL(1)=0._dp 
         dmL(2)=0._dp 
         dmL(3)=0._dp 
         dmL(4)=mc/dm_lam 
         dmL(5)=mb/dm_lam 
         do j=1,5 
            dmR(j)=dmL(j) 
         enddo
c         if(first) then 
            write(6,*) 'Setting up Scalar Yukawa Couplings to DM' 
            write(6,*) 'Couplings are as follows ' 
            write(6,44) 'dm L :',dmL
            write(6,44) 'dm R :',dmR 
c            first=.false.
c         endif
      else
         dmL(1)=0._dp 
         dmL(2)=0._dp 
         dmL(3)=0._dp 
         dmL(4)=mc*inv_v/g_dmq
         dmL(5)=mb*inv_v/g_dmq
         do j=1,5 
            dmR(j)=dmL(j) 
         enddo
c         if(first) then 
            write(6,*) 'Setting up Scalar Yukawa Couplings to DM'
            write(6,*) 'We are in full theory now, have removed factor'
            write(6,*) 'g_dmq from input (since fixed by UV)'
            write(6,*) 'Output below is m_q/v/g_dmq(input) ' 
            write(6,44) 'dm L :',dmL(1:5)
            write(6,44) 'dm R :',dmR(1:5) 
c            first=.false.
c         endif

      endif
      endif
      
 44   format(1x,a10,5(f8.4,' '))
      return 
      end 
      
