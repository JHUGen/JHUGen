      subroutine dipoleconfig(maxdip,dipconfig)
      implicit none
      include 'types.f'
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
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      include 'ptilde.f'
      include 'frag.f'
      include 'taucut.f'
      include 'kpart.f'
      include 'npart.f'
      include 'hdecaymode.f'
      integer:: maxdip,dipconfig(maxd,3),i7,i8
      
c      if (usescet) then
c        if ((case .eq. 'W_1jet') .or. (case .eq. 'Z_1jet')
c     &  .or.(case .eq. 'ggfus1')) then
c           if ((part .eq. klord) .or. (part .eq. kvirt)) then
c             maxdip=2
c             dipconfig(1,:)= (/ 1,5,2 /)
c             dipconfig(2,:)= (/ 2,5,1 /)
c             return
c           endif
c!---------- real branch
c!---------- npart = 3, doing 2->3 splitting as above
c           if (npart .eq. 3) then
c             maxdip=2
c             dipconfig(1,:)= (/ 1,5,2 /)
c             dipconfig(2,:)= (/ 2,5,1 /)
c             return
c!---------- npart = 4, doing 3->4 splitting, c.f. W+1 jet
c           elseif (npart .eq. 4) then
c             maxdip=10
c             dipconfig(1,:)= (/ 1,5,2 /)
c             dipconfig(2,:)= (/ 2,5,1 /)
c             dipconfig(3,:)= (/ 1,6,2 /)
c             dipconfig(4,:)= (/ 2,6,1 /)
c             dipconfig(5,:)= (/ 1,5,6 /)
c             dipconfig(6,:)= (/ 5,6,1 /)
c             dipconfig(7,:)= (/ 1,6,5 /)
c             dipconfig(8,:)= (/ 2,6,5 /)
c             dipconfig(9,:)= (/ 5,6,2 /)
c             dipconfig(10,:)=(/ 2,5,6 /)
c             return
c           else
c             write(6,*) 'Error in dipoleconfig: npart=',npart
c             stop
c           endif
c           write(6,*) 'npart=',npart
cc           write(6,*) 'dipoleconfig not written yet for taucut > 0'
c           stop
c        endif
c      endif
      
      if    ((kcase==kW_only) .or. (kcase==kZ_only)
     &  .or. (kcase==kggfus0)) then
        maxdip=2
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
      elseif ((kcase==kWWqqbr) .or. (kcase==kWWnpol)
     &   .or. (kcase==kWZbbar) .or. (kcase==kZZlept)
     &   .or. (kcase==kWHbbar) .or. (kcase==kZHbbar)
     &   .or. (kcase==kWHgaga) .or. (kcase==kZHgaga)
     &   .or. (kcase==kHZZ_4l) .or. (kcase==kHmZZ4l)
     &   .or. (kcase==kZHlept)) then
        maxdip=2
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
      elseif (kcase==kqq_Hqq) then
        maxdip=8
        dipconfig(1,:)= (/ 1,7,5 /)
        dipconfig(2,:)= (/ 5,7,1 /)
        dipconfig(3,:)= (/ 2,7,6 /)
        dipconfig(4,:)= (/ 6,7,2 /)
        dipconfig(5,:)= (/ 1,5,2 /)
        dipconfig(6,:)= (/ 2,6,1 /)
        dipconfig(7,:)= (/ 1,6,2 /)
        dipconfig(8,:)= (/ 2,7,1 /)
      elseif ((kcase==kW_1jet) .or. (kcase==kZ_1jet)
     &    .or.(kcase==kggfus1) .or. (kcase==kgmgmjt)
     &    .or.(kcase==kHgagaj)) then
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
        if ((frag) .and. (kcase==kgmgmjt)) then
        dipconfig(11,:)=(/ 3,5,1 /)
        dipconfig(12,:)=(/ 3,6,1 /)
        dipconfig(13,:)=(/ 4,5,1 /)
        dipconfig(14,:)=(/ 4,6,1 /)
        maxdip=14
        else
        maxdip=10
        endif
      elseif ((kcase==kW_2jet) .or. (kcase==kZ_2jet)
     &    .or.(kcase==kggfus2)) then
        maxdip=27
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
        dipconfig(3,:)= (/ 1,7,5 /)
        dipconfig(4,:)= (/ 2,7,6 /)
        dipconfig(5,:)= (/ 1,7,6 /)
        dipconfig(6,:)= (/ 2,7,5 /)
        dipconfig(7,:)= (/ 5,7,1 /)
        dipconfig(8,:)= (/ 6,7,2 /)
        dipconfig(9,:)= (/ 6,7,1 /)
        dipconfig(10,:)=(/ 5,7,2 /)
        dipconfig(11,:)=(/ 5,7,6 /)
        dipconfig(12,:)=(/ 6,7,5 /)

        dipconfig(13,:)=(/ 1,6,2 /)
        dipconfig(14,:)=(/ 2,6,1 /)
        dipconfig(15,:)=(/ 1,6,5 /)
        dipconfig(16,:)=(/ 2,6,7 /)
        dipconfig(17,:)=(/ 1,6,7 /)
        dipconfig(18,:)=(/ 2,6,5 /)
        dipconfig(19,:)=(/ 5,6,1 /)
c        dipconfig(20,:)=(/ 7,6,2 /)
c        dipconfig(21,:)=(/ 7,6,1 /)
        dipconfig(20,:)=(/ 5,6,2 /)
        dipconfig(21,:)=(/ 5,6,7 /)
c        dipconfig(24,:)=(/ 7,6,5 /)

        dipconfig(22,:)=(/ 1,5,2 /)
        dipconfig(23,:)=(/ 2,5,1 /)
        dipconfig(24,:)=(/ 1,5,7 /)
        dipconfig(25,:)=(/ 2,5,6 /)
        dipconfig(26,:)=(/ 1,5,6 /)
        dipconfig(27,:)=(/ 2,5,7 /)
c        dipconfig(28,:)=(/ 7,5,1 /)
c        dipconfig(32,:)=(/ 6,5,2 /)
c        dipconfig(33,:)=(/ 6,5,1 /)
c        dipconfig(34,:)=(/ 7,5,2 /)
c        dipconfig(35,:)=(/ 7,5,6 /)
c        dipconfig(36,:)=(/ 6,5,7 /)
      elseif ((kcase==kH_tjet) .or. (kcase==kZ_tjet)) then
        maxdip=8
        dipconfig(1,:)= (/ 1,7,6 /)
        dipconfig(2,:)= (/ 6,7,1 /)
        dipconfig(3,:)= (/ 2,7,6 /)
        dipconfig(4,:)= (/ 6,7,2 /)
        dipconfig(5,:)= (/ 1,7,2 /)
        dipconfig(6,:)= (/ 2,7,1 /)
        dipconfig(7,:)= (/ 1,6,2 /)
        dipconfig(8,:)= (/ 2,6,1 /)
      elseif ((kcase==kH_tdkj) .or. (kcase==kZ_tdkj)) then
        maxdip=8
        dipconfig(1,:)= (/ 1,9,8 /)
        dipconfig(2,:)= (/ 8,9,1 /)
        dipconfig(3,:)= (/ 2,9,8 /)
        dipconfig(4,:)= (/ 8,9,2 /)
        dipconfig(5,:)= (/ 1,9,2 /)
        dipconfig(6,:)= (/ 2,9,1 /)
        dipconfig(7,:)= (/ 1,8,2 /)
        dipconfig(8,:)= (/ 2,8,1 /)
      elseif (kcase==khflgam) then
        maxdip=6
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
        dipconfig(3,:)= (/ 4,5,1 /)
        dipconfig(4,:)= (/ 1,5,4 /)
        dipconfig(5,:)= (/ 4,5,2 /)
        dipconfig(6,:)= (/ 2,5,4 /)
      elseif ((kcase==kgamgam).or.(kcase==kgg2gam)) then 
        dipconfig(1,:)= (/ 1,5,2 /)
        dipconfig(2,:)= (/ 2,5,1 /)
        if (frag) then
        dipconfig(3,:)= (/ 3,5,0 /)
        dipconfig(4,:)= (/ 4,5,0 /)
        maxdip=4
        else
        maxdip=2
        endif
      elseif (kcase==kdirgam) then 
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
      elseif (kcase==kWgamma) then 
         dipconfig(1,:)= (/ 1,6,2 /)
         dipconfig(2,:)= (/ 2,6,1 /)
        if (frag) then
           dipconfig(3,:)= (/ 5,6,2 /)
           maxdip=3
        else
           maxdip=2
        endif
      elseif (kcase==kZgamma) then 
         dipconfig(1,:)= (/ 1,6,2 /)
         dipconfig(2,:)= (/ 2,6,1 /)
        if (frag) then
           dipconfig(3,:)= (/ 5,6,2 /)
           maxdip=3
        else
           maxdip=2
        endif
      elseif (kcase==kZgajet) then 
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
      elseif ((kcase==kW_2gam) .or. (kcase==kZ_2gam)) then 
         dipconfig(1,:)= (/ 1,7,2 /)
         dipconfig(2,:)= (/ 2,7,1 /)
        if (frag) then
           dipconfig(3,:)= (/ 5,7,2 /)
           dipconfig(4,:)= (/ 6,7,2 /)
           maxdip=4
        else
           maxdip=2
        endif
      elseif (kcase==kdm_gam) then 
         dipconfig(1,:)= (/ 1,6,2 /)
         dipconfig(2,:)= (/ 2,6,1 /)
        if (frag) then
           dipconfig(3,:)= (/ 5,6,2 /)
           maxdip=3
        else
           maxdip=2
        endif
      elseif (kcase==ktrigam) then
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
      elseif (kcase==ktottth) then
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        maxdip=2
      elseif (kcase==kW_cjet) then
        dipconfig(1,:)= (/ 1,6,2 /)
        dipconfig(2,:)= (/ 2,6,1 /)
        maxdip=2
      elseif ((kcase==kWW_jet) .or. (kcase==kWH1jet)
     &    .or.(kcase==kZH1jet)) then
        if (hdecaymode == 'wpwm') then
          i7=9
          i8=10
        else
          i7=7
          i8=8
        endif
        dipconfig(1,:)= (/ 1,i7,2 /)
        dipconfig(2,:)= (/ 2,i7,1 /)
        dipconfig(3,:)= (/ 1,i8,2 /)
        dipconfig(4,:)= (/ 2,i8,1 /)
        dipconfig(5,:)= (/ 1,i7,i8 /)
        dipconfig(6,:)= (/ i7,i8,1 /)
        dipconfig(7,:)= (/ 1,i8,i7 /)
        dipconfig(8,:)= (/ 2,i8,i7 /)
        dipconfig(9,:)= (/ i7,i8,2 /)
        dipconfig(10,:)=(/ 2,i7,i8 /)
        maxdip=10
      elseif(kcase==kfourga) then 
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
      
