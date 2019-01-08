      subroutine Carraysetup
      implicit none
      include 'pvCnames.f'
      include 'Carraydef.f'
      include 'Carrays.f'  


      ci(1)=cc1
      ci(2)=cc2


      cii(z2(1,1))=cc11
      cii(z2(1,2))=cc12
      cii(z2(2,2))=cc22
       
      ciii(z3(1,1,1))=cc111
      ciii(z3(1,1,2))=cc112
      ciii(z3(1,2,2))=cc122
      ciii(z3(2,2,2))=cc222

      ciiii(z4(1,1,1,1))=cc1111
      ciiii(z4(1,1,1,2))=cc1112
      ciiii(z4(1,1,2,2))=cc1122
      ciiii(z4(1,2,2,2))=cc1222
      ciiii(z4(2,2,2,2))=cc2222

      ciiiii(z5(1,1,1,1,1))=cc11111
      ciiiii(z5(1,1,1,1,2))=cc11112
      ciiiii(z5(1,1,1,2,2))=cc11122
      ciiiii(z5(1,1,2,2,2))=cc11222
      ciiiii(z5(1,2,2,2,2))=cc12222
      ciiiii(z5(2,2,2,2,2))=cc22222

      ciiiiii(z6(1,1,1,1,1,1))=cc111111
      ciiiiii(z6(1,1,1,1,1,2))=cc111112
      ciiiiii(z6(1,1,1,1,2,2))=cc111122
      ciiiiii(z6(1,1,1,2,2,2))=cc111222
      ciiiiii(z6(1,1,2,2,2,2))=cc112222
      ciiiiii(z6(1,2,2,2,2,2))=cc122222
      ciiiiii(z6(2,2,2,2,2,2))=cc222222

      ciiiiiii(z7(1,1,1,1,1,1,1))=cc1111111
      ciiiiiii(z7(1,1,1,1,1,1,2))=cc1111112
      ciiiiiii(z7(1,1,1,1,1,2,2))=cc1111122
      ciiiiiii(z7(1,1,1,1,2,2,2))=cc1111222
      ciiiiiii(z7(1,1,1,2,2,2,2))=cc1112222
      ciiiiiii(z7(1,1,2,2,2,2,2))=cc1122222
      ciiiiiii(z7(1,2,2,2,2,2,2))=cc1222222
      ciiiiiii(z7(2,2,2,2,2,2,2))=cc2222222


      czzi(1)=cc001
      czzi(2)=cc002

      czzii(z2(1,1))=cc0011
      czzii(z2(1,2))=cc0012
      czzii(z2(2,2))=cc0022

      czziii(z3(1,1,1))=cc00111
      czziii(z3(1,1,2))=cc00112
      czziii(z3(1,2,2))=cc00122
      czziii(z3(2,2,2))=cc00222

      czziiii(z4(1,1,1,1))=cc001111
      czziiii(z4(1,1,1,2))=cc001112
      czziiii(z4(1,1,2,2))=cc001122
      czziiii(z4(1,2,2,2))=cc001222
      czziiii(z4(2,2,2,2))=cc002222

      czziiiii(z5(1,1,1,1,1))=cc0011111
      czziiiii(z5(1,1,1,1,2))=cc0011112
      czziiiii(z5(1,1,1,2,2))=cc0011122
      czziiiii(z5(1,1,2,2,2))=cc0011222
      czziiiii(z5(1,2,2,2,2))=cc0012222
      czziiiii(z5(2,2,2,2,2))=cc0022222

      czzzzi(1)=cc00001
      czzzzi(2)=cc00002

      czzzzii(z2(1,1))=cc000011
      czzzzii(z2(1,2))=cc000012
      czzzzii(z2(2,2))=cc000022

      czzzziii(z3(1,1,1))=cc0000111
      czzzziii(z3(1,1,2))=cc0000112
      czzzziii(z3(1,2,2))=cc0000122
      czzzziii(z3(2,2,2))=cc0000222

      czzzzzzi(1)=cc0000001
      czzzzzzi(2)=cc0000002


      return
      end
