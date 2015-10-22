#define WorkPath "/home/schulze/1git_JHUGen/VH_NLO/"
#define setP1P1 "0"
#define setP2P2 "0"
#define setP3P3 "MZ^2"
#define setP4P4 "MH^2"
#define setDST   "4"
#define setDSTm4 "0"
#define InterfereDiags1 "1,1"
#define InterfereDiags2 "1,1"



#include /home/schulze/lib/FeynArtsToForm/header2.frm


#include `WorkPath'Virt_qqbVH_input.frm

* external up quark mass
id MU=0;
argument;
   id MU=0;
endargument;

* internal up quark mass
id MQUGen5=MT;
argument;
   id MQUGen5=MT;
endargument;

* internal dn quark mass
id MQDGen5=0;
argument;
   id MQDGen5=0;
endargument;



id IndexDelta(Col1?,Col2?) = 1;
id DID(Col1?)=1;
id SumOver(?args) = 1;



#include `WorkPath'Virt_qqbVH_calc.frm




Format Mathematica;
Print;
Bracket IZ,SME;
.sort;


#write <`WorkPath'Virt_qqbVH_output.dat> "NumAmps = `NumAmps';\n";
#write <`WorkPath'Virt_qqbVH_output.dat> "AmpList = \"  `AmpList'\";\n";

#do TheAmp = `AmpList'
  #write <`WorkPath'Virt_qqbVH_output.dat> "M`TheAmp' = (%E);\n ", `TheAmp';
#enddo




.end
