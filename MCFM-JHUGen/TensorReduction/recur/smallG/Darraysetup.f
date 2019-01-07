      subroutine Darraysetup
      implicit none
      include 'pvDnames.f'
      include 'Darraydef.f'
      include 'Darrays.f'  


      di(1)=dd1
      di(2)=dd2
      di(3)=dd3


      dii(z2(1,1))=dd11
      dii(z2(1,2))=dd12
      dii(z2(1,3))=dd13
      dii(z2(2,2))=dd22
      dii(z2(2,3))=dd23
      dii(z2(3,3))=dd33
       
      diii(z3(1,1,1))=dd111
      diii(z3(1,1,2))=dd112
      diii(z3(1,1,3))=dd113

      diii(z3(1,2,2))=dd122
      diii(z3(1,2,3))=dd123
      diii(z3(1,3,3))=dd133

      diii(z3(2,2,2))=dd222
      diii(z3(2,2,3))=dd223
      diii(z3(2,3,3))=dd233
      diii(z3(3,3,3))=dd333

      diiii(z4(1,1,1,1))=dd1111
      diiii(z4(1,1,1,2))=dd1112
      diiii(z4(1,1,1,3))=dd1113

      diiii(z4(1,1,2,2))=dd1122
      diiii(z4(1,1,2,3))=dd1123
      diiii(z4(1,1,3,3))=dd1133

      diiii(z4(1,2,2,2))=dd1222
      diiii(z4(1,2,2,3))=dd1223
      diiii(z4(1,2,3,3))=dd1233

      diiii(z4(1,3,3,3))=dd1333
      diiii(z4(2,2,2,2))=dd2222
      diiii(z4(2,2,2,3))=dd2223

      diiii(z4(2,2,3,3))=dd2233
      diiii(z4(2,3,3,3))=dd2333
      diiii(z4(3,3,3,3))=dd3333


      diiiii(z5(1,1,1,1,1))=dd11111
      diiiii(z5(1,1,1,1,2))=dd11112
      diiiii(z5(1,1,1,1,3))=dd11113

      diiiii(z5(1,1,1,2,2))=dd11122
      diiiii(z5(1,1,1,2,3))=dd11123
      diiiii(z5(1,1,1,3,3))=dd11133

      diiiii(z5(1,1,2,2,2))=dd11222
      diiiii(z5(1,1,2,2,3))=dd11223
      diiiii(z5(1,1,2,3,3))=dd11233

      diiiii(z5(1,1,3,3,3))=dd11333
      diiiii(z5(1,2,2,2,2))=dd12222
      diiiii(z5(1,2,2,2,3))=dd12223

      diiiii(z5(1,2,2,3,3))=dd12233
      diiiii(z5(1,2,3,3,3))=dd12333
      diiiii(z5(1,3,3,3,3))=dd13333

      diiiii(z5(2,2,2,2,2))=dd22222
      diiiii(z5(2,2,2,2,3))=dd22223
      diiiii(z5(2,2,2,3,3))=dd22233

      diiiii(z5(2,2,3,3,3))=dd22333
      diiiii(z5(2,3,3,3,3))=dd23333
      diiiii(z5(3,3,3,3,3))=dd33333


      diiiiii(z6(1,1,1,1,1,1))=dd111111
      diiiiii(z6(1,1,1,1,1,2))=dd111112

      diiiiii(z6(1,1,1,1,1,3))=dd111113
      diiiiii(z6(1,1,1,1,2,2))=dd111122

      diiiiii(z6(1,1,1,1,2,3))=dd111123
      diiiiii(z6(1,1,1,1,3,3))=dd111133

      diiiiii(z6(1,1,1,2,2,2))=dd111222
      diiiiii(z6(1,1,1,2,2,3))=dd111223

      diiiiii(z6(1,1,1,2,3,3))=dd111233
      diiiiii(z6(1,1,1,3,3,3))=dd111333

      diiiiii(z6(1,1,2,2,2,2))=dd112222
      diiiiii(z6(1,1,2,2,2,3))=dd112223

      diiiiii(z6(1,1,2,2,3,3))=dd112233
      diiiiii(z6(1,1,2,3,3,3))=dd112333

      diiiiii(z6(1,1,3,3,3,3))=dd113333
      diiiiii(z6(1,2,2,2,2,2))=dd122222

      diiiiii(z6(1,2,2,2,2,3))=dd122223
      diiiiii(z6(1,2,2,2,3,3))=dd122233

      diiiiii(z6(1,2,2,3,3,3))=dd122333
      diiiiii(z6(1,2,3,3,3,3))=dd123333

      diiiiii(z6(1,3,3,3,3,3))=dd133333
      diiiiii(z6(2,2,2,2,2,2))=dd222222

      diiiiii(z6(2,2,2,2,2,3))=dd222223
      diiiiii(z6(2,2,2,2,3,3))=dd222233

      diiiiii(z6(2,2,2,3,3,3))=dd222333
      diiiiii(z6(2,2,3,3,3,3))=dd223333

      diiiiii(z6(2,3,3,3,3,3))=dd233333
      diiiiii(z6(3,3,3,3,3,3))=dd333333


      diiiiiii(z7(1,1,1,1,1,1,1))=dd1111111      

      diiiiiii(z7(1,1,1,1,1,1,2))=dd1111112     
      diiiiiii(z7(1,1,1,1,1,1,3))=dd1111113  
          
      diiiiiii(z7(1,1,1,1,1,2,2))=dd1111122     
      diiiiiii(z7(1,1,1,1,1,2,3))=dd1111123      
      diiiiiii(z7(1,1,1,1,1,3,3))=dd1111133      

      diiiiiii(z7(1,1,1,1,2,2,2))=dd1111222     
      diiiiiii(z7(1,1,1,1,2,2,3))=dd1111223     
      diiiiiii(z7(1,1,1,1,2,3,3))=dd1111233     
      diiiiiii(z7(1,1,1,1,3,3,3))=dd1111333     

      diiiiiii(z7(1,1,1,2,2,2,2))=dd1112222     
      diiiiiii(z7(1,1,1,2,2,2,3))=dd1112223     
      diiiiiii(z7(1,1,1,2,2,3,3))=dd1112233     
      diiiiiii(z7(1,1,1,2,3,3,3))=dd1112333     
      diiiiiii(z7(1,1,1,3,3,3,3))=dd1113333     

      diiiiiii(z7(1,1,2,2,2,2,2))=dd1122222     
      diiiiiii(z7(1,1,2,2,2,2,3))=dd1122223     
      diiiiiii(z7(1,1,2,2,2,3,3))=dd1122233     
      diiiiiii(z7(1,1,2,2,3,3,3))=dd1122333     
      diiiiiii(z7(1,1,2,3,3,3,3))=dd1123333     
      diiiiiii(z7(1,1,3,3,3,3,3))=dd1133333     

      diiiiiii(z7(1,2,2,2,2,2,2))=dd1222222     
      diiiiiii(z7(1,2,2,2,2,2,3))=dd1222223     
      diiiiiii(z7(1,2,2,2,2,3,3))=dd1222233     
      diiiiiii(z7(1,2,2,2,3,3,3))=dd1222333     
      diiiiiii(z7(1,2,2,3,3,3,3))=dd1223333     
      diiiiiii(z7(1,2,3,3,3,3,3))=dd1233333     
      diiiiiii(z7(1,3,3,3,3,3,3))=dd1333333     

      diiiiiii(z7(2,2,2,2,2,2,2))=dd2222222     
      diiiiiii(z7(2,2,2,2,2,2,3))=dd2222223     
      diiiiiii(z7(2,2,2,2,2,3,3))=dd2222233     
      diiiiiii(z7(2,2,2,2,3,3,3))=dd2222333     
      diiiiiii(z7(2,2,2,3,3,3,3))=dd2223333     
      diiiiiii(z7(2,2,3,3,3,3,3))=dd2233333     
      diiiiiii(z7(2,3,3,3,3,3,3))=dd2333333     
      diiiiiii(z7(3,3,3,3,3,3,3))=dd3333333     


      dzzi(1)=dd001
      dzzi(2)=dd002
      dzzi(3)=dd003

      dzzii(z2(1,1))=dd0011
      dzzii(z2(1,2))=dd0012
      dzzii(z2(1,3))=dd0013
      dzzii(z2(2,2))=dd0022
      dzzii(z2(2,3))=dd0023
      dzzii(z2(3,3))=dd0033

      dzziii(z3(1,1,1))=dd00111
      dzziii(z3(1,1,2))=dd00112
      dzziii(z3(1,1,3))=dd00113
      dzziii(z3(1,2,2))=dd00122
      dzziii(z3(1,2,3))=dd00123
      dzziii(z3(1,3,3))=dd00133
      dzziii(z3(2,2,2))=dd00222
      dzziii(z3(2,2,3))=dd00223
      dzziii(z3(2,3,3))=dd00233
      dzziii(z3(3,3,3))=dd00333

      dzziiii(z4(1,1,1,1))=dd001111
      dzziiii(z4(1,1,1,2))=dd001112
      dzziiii(z4(1,1,1,3))=dd001113

      dzziiii(z4(1,1,2,2))=dd001122
      dzziiii(z4(1,1,2,3))=dd001123
      dzziiii(z4(1,1,3,3))=dd001133

      dzziiii(z4(1,2,2,2))=dd001222
      dzziiii(z4(1,2,2,3))=dd001223
      dzziiii(z4(1,2,3,3))=dd001233

      dzziiii(z4(1,3,3,3))=dd001333
      dzziiii(z4(2,2,2,2))=dd002222
      dzziiii(z4(2,2,2,3))=dd002223

      dzziiii(z4(2,2,3,3))=dd002233
      dzziiii(z4(2,3,3,3))=dd002333
      dzziiii(z4(3,3,3,3))=dd003333

      dzziiiii(z5(1,1,1,1,1))=dd0011111
      
      dzziiiii(z5(1,1,1,1,2))=dd0011112
      dzziiiii(z5(1,1,1,1,3))=dd0011113

      dzziiiii(z5(1,1,1,2,2))=dd0011122
      dzziiiii(z5(1,1,1,2,3))=dd0011123
      dzziiiii(z5(1,1,1,3,3))=dd0011133

      dzziiiii(z5(1,1,2,2,2))=dd0011222
      dzziiiii(z5(1,1,2,2,3))=dd0011223
      dzziiiii(z5(1,1,2,3,3))=dd0011233
      dzziiiii(z5(1,1,3,3,3))=dd0011333

      dzziiiii(z5(1,2,2,2,2))=dd0012222
      dzziiiii(z5(1,2,2,2,3))=dd0012223
      dzziiiii(z5(1,2,2,3,3))=dd0012233
      dzziiiii(z5(1,2,3,3,3))=dd0012333
      dzziiiii(z5(1,3,3,3,3))=dd0013333

      dzziiiii(z5(2,2,2,2,2))=dd0022222
      dzziiiii(z5(2,2,2,2,3))=dd0022223
      dzziiiii(z5(2,2,2,3,3))=dd0022233
      dzziiiii(z5(2,2,3,3,3))=dd0022333
      dzziiiii(z5(2,3,3,3,3))=dd0023333
      dzziiiii(z5(3,3,3,3,3))=dd0033333

      dzzzzi(1)=dd00001
      dzzzzi(2)=dd00002
      dzzzzi(3)=dd00003

      dzzzzii(z2(1,1))=dd000011
      dzzzzii(z2(1,2))=dd000012
      dzzzzii(z2(1,3))=dd000013
      dzzzzii(z2(2,2))=dd000022
      dzzzzii(z2(2,3))=dd000023
      dzzzzii(z2(3,3))=dd000033

      dzzzziii(z3(1,1,1))=dd0000111

      dzzzziii(z3(1,1,2))=dd0000112
      dzzzziii(z3(1,1,3))=dd0000113

      dzzzziii(z3(1,2,2))=dd0000122
      dzzzziii(z3(1,2,3))=dd0000123
      dzzzziii(z3(1,3,3))=dd0000133

      dzzzziii(z3(2,2,2))=dd0000222
      dzzzziii(z3(2,2,3))=dd0000223
      dzzzziii(z3(2,3,3))=dd0000233
      dzzzziii(z3(3,3,3))=dd0000333

      dzzzzzzi(1)=dd0000001
      dzzzzzzi(2)=dd0000002
      dzzzzzzi(3)=dd0000003

      return
      end
