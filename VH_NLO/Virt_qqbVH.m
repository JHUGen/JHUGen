InitPath =  $HomeDirectory<>"/lib/MathLib/";
FeynArtsPath =  $HomeDirectory<>"/lib/FeynArts-3.9/";
FeynArtsToFormPath =  $HomeDirectory<>"/lib/FeynArtsToForm/";
ProjectPath =  $HomeDirectory<>"/1git_JHUGen/VH_NLO/";
Unprotect[Diagonal];
Clear[Diagonal];

Get[ FeynArtsPath<>"FeynArts.m" ];
Get[ FeynArtsToFormPath<>"ExtractAmp.m" ];


TheProcess = {F[3, {1}],-F[3, {1}]} -> {V[2],S[1]};

TheLevel = {Classes};
TheGenericModel = "Lorentz";
TheModel = SMQCD;
PaintIt  = True;
WriteAmp = True;

SetOptions[CreateTopologies, ExcludeTopologies -> {Tadpoles,WFCorrections}];


SetOptions[InsertFields,
    	InsertionLevel -> TheLevel,
    	GenericModel -> TheGenericModel ,
    	Model -> TheModel,
        ExcludeParticles -> { SV[2],SV[3],S[1],S[2],S[3],V[1],V[3],U[1],U[2],U[3],U[4]  },
        LastSelections-> {},
        Restrictions -> { NoLightFHCoupling }
];

SetOptions[Paint, Numbering->Simple, ColumnsXRows->{5,7} ];
CKM = IndexDelta;



(* generate Virt level diagrams *)
TheVirtTop  = CreateTopologies[ 1, 2 -> 2 ];



TheVirtDiag = InsertFields[ TheVirtTop, TheProcess ];

If[PaintIt,
   Paint[TheVirtDiag, DisplayFunction->(Display[ProjectPath<>"Virt_qqbVH.ps",#]&) ];
];

TheVirtAmp   = CreateFeynAmp[ TheVirtDiag,   PreFactor-> 1 ];

Do[  TheVirtAmp[[i,3]]  *=DID[i] ,{i,1,Length[ TheVirtAmp]} ];

(* removing the el.weak correction diagrams *)
TheVirtAmp2 = Delete[TheVirtAmp,{{1},{2},{5},{6},{7},{9},{10},{11},{12},{13}}];


MyVirtAmp   = ExtractAmp[ TheVirtAmp2 ];
If[WriteAmp,
   AmpToFORM[ Join[MyVirtAmp] ,ProjectPath<>"Virt_qqbVH_input.frm"];
];

