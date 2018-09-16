      subroutine checkEGZres
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'masses.f'
      include 'scale.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'nflav.f'
      include 'scheme.f'
      include 'deltar.f'
      real(dp):: p(mxpart,4),Hqarbvsqanal,Hqaqavsqanal,
     & HAQggvsqanal,Hggggvsqanal,check,check_a,check_b,check_i

c--- Set up momenta, Higgs mass and scale for special point

c--- These are per the paper
c      p(1,4)=+0.30674037867_dp
c      p(1,1)=-0.17738694693_dp
c      p(1,2)=-0.01664472021_dp
c      p(1,3)=-0.24969277974_dp
c
c      p(2,4)=+0.34445032281_dp
c      p(2,1)=+0.14635282800_dp
c      p(2,2)=-0.10707762397_dp
c      p(2,3)=+0.29285022975_dp
c
c      p(3,4)=0.5_dp
c      p(3,1)=0._dp
c      p(3,2)=0._dp
c      p(3,3)=0._dp
c
c      p(4,4)=0.5_dp
c      p(4,1)=0._dp
c      p(4,2)=0._dp
c      p(4,3)=0._dp
c
c      p(5,4)=+0.22091667641_dp
c      p(5,1)=+0.08911915938_dp
c      p(5,2)=+0.19733901856_dp
c      p(5,3)=+0.04380941793_dp
c
c      p(6,4)=+0.12789262211_dp
c      p(6,1)=-0.05808504045_dp
c      p(6,2)=-0.07361667438_dp
c      p(6,3)=-0.08696686795_dp

c--- These are direct from Giulia's code (more s.f.)
      p(1,4)=+0.306740378669096_dp
      p(1,1)=-0.1773869469327211_dp
      p(1,2)=-0.01664472021061108_dp
      p(1,3)=-0.2496927797375384_dp

      p(2,4)=+0.3444503228094516_dp
      p(2,1)=+0.1463528280021222_dp
      p(2,2)=-0.1070776239705892_dp
      p(2,3)=+0.2928502297491496_dp

      p(3,4)=0.5_dp
      p(3,1)=0._dp
      p(3,2)=0._dp
      p(3,3)=0._dp

      p(4,4)=0.5_dp
      p(4,1)=0._dp
      p(4,2)=0._dp
      p(4,3)=0._dp

      p(5,4)=+0.2209166764078266_dp
      p(5,1)=+0.08911915938215637_dp
      p(5,2)=+0.1973390185605316_dp
      p(5,3)=+0.04380941793341205_dp

      p(6,4)=+0.1278926221136257_dp
      p(6,1)=-0.05808504045155738_dp
      p(6,2)=-0.0736166743793313_dp
      p(6,3)=-0.0869668679450233_dp

      hmass=1._dp
      scale=1._dp
      musq=1._dp
      nflav=5
      scheme='tH-V'
      deltar=1._dp

      call dotem(6,p,s)
      call spinoru(6,p,za,zb)

c--- only comparing finite piece
      epinv=0._dp

      call HqarbLO(6,2,5,1,check)
      write(6,99) 'Proc A: EGZ check B',check
      check=Hqarbvsqanal(5,6,1,2)
      write(6,99) 'Proc A: EGZ check V',check

      call HqaqaLO(6,2,5,1,check,check_a,check_b,check_i)
      write(6,99) 'Proc B: EGZ check B',check
      check=Hqarbvsqanal(5,6,1,2)+Hqarbvsqanal(1,6,5,2)
     &     +Hqaqavsqanal(5,6,1,2)
      write(6,99) 'Proc B: EGZ check V',check

      call HQAggLO(1,2,5,6,check,check_a,check_b,check_i)
      write(6,99) 'Proc C: EGZ check B',check
      check=HAQggvsqanal(1,2,5,6)
      write(6,99) 'Proc C: EGZ check V',check

      call HggggLO(1,2,5,6,check,check_a,check_b,check_i)
      write(6,99) 'Proc D: EGZ check B',check
      check=Hggggvsqanal(1,2,5,6)
      write(6,99) 'Proc D: EGZ check V',check

c      pause

      return

   99 format(a20,e24.15)

      end

