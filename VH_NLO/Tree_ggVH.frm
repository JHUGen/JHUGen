#define WorkPath "/afs/cern.ch/user/m/maschulz/1git_JHUGen/VH_NLO/"
#define setP1P1 "0"
#define setP2P2 "0"
#define setP3P3 "MZ^2"
#define setP4P4 "MH^2"
#define setDST   "4"
#define setDSTm4 "0"
#define InterfereDiags1 "1,1"
#define InterfereDiags2 "1,1"



#include /afs/cern.ch/user/m/maschulz/lib/FeynArtsToForm/header2.frm


#include `WorkPath'Tree_ggVH_input.frm

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



#include `WorkPath'Tree_ggVH_calc.frm




Format Mathematica;
Print;
Bracket SME,IZ,ICol,EL,GS,SW,MW,Pi;
.sort;


#write <`WorkPath'Tree_ggVH_output.dat> "NumAmps = `NumAmps';\n";
#write <`WorkPath'Tree_ggVH_output.dat> "AmpList = \"  `AmpList'\";\n";

#do TheAmp = `AmpList'
  #write <`WorkPath'Tree_ggVH_output.dat> "M`TheAmp' = (%E);\n ", `TheAmp';
#enddo




.end
