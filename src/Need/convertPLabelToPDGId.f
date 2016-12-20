      subroutine convertPLabelsToPDGIds()
      implicit none
      include 'mxpart.f'
      include 'plabel.f'
      include 'pid_pdg.f'
      integer convertPLabelToPDGId
      integer j
      do j=1,mxpart
         pid_pdg(j)=convertPLabelToPDGId(j)
      enddo
      return
      end

      function convertPLabelToPDGId(j)
      implicit none
      include 'mxpart.f'
      include 'plabel.f'
      integer convertPLabelToPDGId
      integer j

      convertPLabelToPDGId=-9000
      if (j.le.mxpart) then
         if (plabel(j) .eq. 'pp') then
         convertPLabelToPDGId=0
         elseif (plabel(j) .eq. 'qj') then
         convertPLabelToPDGId=0
         elseif (plabel(j) .eq. 'ig') then
         convertPLabelToPDGId=21
         elseif (plabel(j) .eq. 'gj') then
         convertPLabelToPDGId=21

         elseif (plabel(j) .eq. 'uq') then
         convertPLabelToPDGId=2
         elseif (plabel(j) .eq. 'dq') then
         convertPLabelToPDGId=1
         elseif (plabel(j) .eq. 'cq') then
         convertPLabelToPDGId=4
         elseif (plabel(j) .eq. 'sq') then
         convertPLabelToPDGId=3
         elseif (plabel(j) .eq. 'bq') then
         convertPLabelToPDGId=5
         elseif (plabel(j) .eq. 'nl') then
         convertPLabelToPDGId=12
         elseif (plabel(j) .eq. 'nm') then
         convertPLabelToPDGId=14
         elseif (plabel(j) .eq. 'nt') then
         convertPLabelToPDGId=16
         elseif (plabel(j) .eq. 'el') then
         convertPLabelToPDGId=11
         elseif (plabel(j) .eq. 'ml') then
         convertPLabelToPDGId=13
         elseif (plabel(j) .eq. 'tl') then
         convertPLabelToPDGId=15

         elseif (plabel(j) .eq. 'ua') then
         convertPLabelToPDGId=-2
         elseif (plabel(j) .eq. 'da') then
         convertPLabelToPDGId=-1
         elseif (plabel(j) .eq. 'ca') then
         convertPLabelToPDGId=-4
         elseif (plabel(j) .eq. 'sa') then
         convertPLabelToPDGId=-3
         elseif (plabel(j) .eq. 'ba') then
         convertPLabelToPDGId=-5
         elseif (plabel(j) .eq. 'na') then
         convertPLabelToPDGId=-12
         elseif (plabel(j) .eq. 'ea') then
         convertPLabelToPDGId=-11
         elseif (plabel(j) .eq. 'ma') then
         convertPLabelToPDGId=-13
         elseif (plabel(j) .eq. 'ta') then
         convertPLabelToPDGId=-15
         endif
      endif

      return
      end
