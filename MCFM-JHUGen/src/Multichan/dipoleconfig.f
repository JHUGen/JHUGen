      subroutine dipoleconfig(maxdip,dipconfig)
************************************************************************
*                                                                      *
*     Given a process, return an array that specifies all the          *
*     dipole configurations that are present in the real subtractions  *
*                                                                      *
*          maxdip: maximum no. of dipoles for this process             *
*       dipconfig: array of dipole configs (i,j,k) in standard MCFM    *
*                  notation, i.e. (emitter, emitted, spectator)        *
*                                                                      *
*     Author: J.M.Campbell                                             *
*       Date: 19th March 2009                                          *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'process.f'
      include 'ptilde.f'
      include 'frag.f'
      integer maxdip,dipconfig(maxd,3)
      
      if    ((case .eq. 'W_only') .or. (case .eq. 'Z_only')) then
        maxdip=2
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
      elseif ((case .eq. 'WWqqbr') .or. (case .eq. 'WWnpol')
     .   .or. (case .eq. 'WZbbar') .or. (case .eq. 'ZZlept')
     .   .or. (case .eq. 'HZZ_4l') .or. (case .eq. 'HmZZ4l')) then
        maxdip=2
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
      elseif (case .eq. 'qq_Hqq') then
        maxdip=8
        dipconfig(1,:)= (/ 1,7,5 /)
        dipconfig(2,:)= (/ 5,7,1 /)
        dipconfig(3,:)= (/ 2,7,6 /)
        dipconfig(4,:)= (/ 6,7,2 /)
        dipconfig(5,:)= (/ 1,5,2 /)
        dipconfig(6,:)= (/ 2,6,1 /)
        dipconfig(7,:)= (/ 1,6,2 /)
        dipconfig(8,:)= (/ 2,7,1 /)
      elseif ((case .eq. 'W_1jet') .or. (case .eq. 'Z_1jet')
     &    .or.(case .eq. 'gmgmjt')) then
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
        dipconfig(3,:)= (/ 1,6,2 /)
        dipconfig(4,:)= (/ 2,6,1 /)
        dipconfig(5,:)= (/ 1,5,6 /)
        dipconfig(6,:)= (/ 5,6,1 /)
        dipconfig(7,:)= (/ 1,6,5 /)
        dipconfig(8,:)= (/ 2,6,5 /)
        dipconfig(9,:)= (/ 5,6,2 /)
        dipconfig(10,:)=(/ 2,5,6 /)
        if ((frag) .and. (case .eq. 'gmgmjt')) then
        dipconfig(11,:)=(/ 3,5,1 /)
        dipconfig(12,:)=(/ 3,6,1 /)
        dipconfig(13,:)=(/ 4,5,1 /)
        dipconfig(14,:)=(/ 4,6,1 /)
        maxdip=14
        else
        maxdip=10
        endif
      elseif ((case .eq. 'H_tjet') .or. (case .eq. 'Z_tjet')) then
        maxdip=8
        dipconfig(1,:)= (/ 1,7,6 /)
        dipconfig(2,:)= (/ 6,7,1 /)
        dipconfig(3,:)= (/ 2,7,6 /)
        dipconfig(4,:)= (/ 6,7,2 /)
        dipconfig(5,:)= (/ 1,7,2 /)
        dipconfig(6,:)= (/ 2,7,1 /)
        dipconfig(7,:)= (/ 1,6,2 /)
        dipconfig(8,:)= (/ 2,6,1 /)
      elseif ((case .eq. 'H_tdkj') .or. (case .eq. 'Z_tdkj')) then
        maxdip=8
        dipconfig(1,:)= (/ 1,9,8 /)
        dipconfig(2,:)= (/ 8,9,1 /)
        dipconfig(3,:)= (/ 2,9,8 /)
        dipconfig(4,:)= (/ 8,9,2 /)
        dipconfig(5,:)= (/ 1,9,2 /)
        dipconfig(6,:)= (/ 2,9,1 /)
        dipconfig(7,:)= (/ 1,8,2 /)
        dipconfig(8,:)= (/ 2,8,1 /)
      elseif (case .eq. 'hflgam') then
        maxdip=6
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
        dipconfig(3,:)= (/ 4,5,1 /)
        dipconfig(4,:)= (/ 1,5,4 /)
        dipconfig(5,:)= (/ 4,5,2 /)
        dipconfig(6,:)= (/ 2,5,4 /)
      elseif (case .eq.'gamgam') then 
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
        if (frag) then
        dipconfig(3,:)= (/ 3,5,0 /)
        dipconfig(4,:)= (/ 4,5,0 /)
        maxdip=4
        else
        maxdip=2
        endif
      elseif (case .eq.'dirgam') then 
        dipconfig(1,:)= (/ 1,4,2 /)
        dipconfig(2,:)= (/ 2,4,1 /)
        dipconfig(3,:)= (/ 1,5,2 /)
        dipconfig(4,:)= (/ 2,5,1 /)
        dipconfig(5,:)= (/ 4,5,2 /)
        dipconfig(6,:)= (/ 4,5,1 /)
        if (frag) then
        dipconfig(7,:)= (/ 3,4,1 /)
        dipconfig(8,:)= (/ 3,5,1 /)
        maxdip=8
        else
        maxdip=6
        endif
      elseif (case .eq. 'Wgamma') then 
         dipconfig(1,:)= (/ 1,6,2 /)
         dipconfig(2,:)= (/ 2,6,1 /)
        if (frag) then
           dipconfig(3,:)= (/ 5,6,2 /)
           maxdip=3
        else
           maxdip=2
        endif
      elseif (case .eq. 'Zgamma') then 
         dipconfig(1,:)= (/ 1,6,2 /)
         dipconfig(2,:)= (/ 2,6,1 /)
        if (frag) then
           dipconfig(3,:)= (/ 5,6,2 /)
           maxdip=3
        else
           maxdip=2
        endif
      elseif (case .eq. 'Zgajet') then 
         dipconfig(1,:)= (/ 1,6,2 /)
         dipconfig(2,:)= (/ 2,6,1 /)
         dipconfig(3,:)= (/ 1,7,2 /)
         dipconfig(4,:)= (/ 2,7,1 /)
         dipconfig(5,:)= (/ 6,7,1 /)
         dipconfig(6,:)= (/ 6,7,2 /)
        if (frag) then
           dipconfig(7,:)= (/ 5,6,2 /)
           dipconfig(8,:)= (/ 5,7,2 /)
           maxdip=8
        else
           maxdip=6
        endif
      elseif ((case .eq. 'W_2gam') .or. (case .eq. 'Z_2gam')) then 
         dipconfig(1,:)= (/ 1,7,2 /)
         dipconfig(2,:)= (/ 2,7,1 /)
        if (frag) then
           dipconfig(3,:)= (/ 5,7,2 /)
           dipconfig(4,:)= (/ 6,7,2 /)
           maxdip=4
        else
           maxdip=2
        endif
      elseif (case .eq. 'dm_gam') then 
         dipconfig(1,:)= (/ 1,6,2 /)
         dipconfig(2,:)= (/ 2,6,1 /)
        if (frag) then
           dipconfig(3,:)= (/ 5,6,2 /)
           maxdip=3
        else
           maxdip=2
        endif
      elseif (case .eq. 'trigam') then
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        if (frag) then
        dipconfig(3,:)= (/ 3,6,2 /)
        dipconfig(4,:)= (/ 4,6,2 /)
        dipconfig(5,:)= (/ 5,6,2 /)
        maxdip=5
        else
        maxdip=2
        endif
      elseif (case .eq. 'tottth') then
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        maxdip=2
      elseif (case .eq. 'W_cjet') then
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        maxdip=2
      elseif (case .eq. 'WW_jet') then
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
        dipconfig(3,:)= (/ 1,8,2 /)
        dipconfig(4,:)= (/ 2,8,1 /)
        dipconfig(5,:)= (/ 1,7,8 /)
        dipconfig(6,:)= (/ 7,8,1 /)
        dipconfig(7,:)= (/ 1,8,7 /)
        dipconfig(8,:)= (/ 2,8,7 /)
        dipconfig(9,:)= (/ 7,8,2 /)
        dipconfig(10,:)=(/ 2,7,8 /)
        maxdip=10
      elseif(case.eq.'fourga') then 
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
        if (frag) then
           dipconfig(3,:)= (/ 3,7,1 /)
           dipconfig(4,:)= (/ 4,7,1 /)
           dipconfig(5,:)= (/ 5,7,1 /)
           dipconfig(6,:)= (/ 6,7,1 /)
           maxdip=6
        else
           maxdip=2
        endif
      else
        write(6,*) 'Dipole configurations for this process not'
        write(6,*) 'properly specified in dipoleconfig.f'
        stop
      endif
      
      return
      end
      
