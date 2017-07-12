{
  TString LIBMCFMPATH = "../data/slc6_amd64_gcc530/";
  TString LIBMCFM="libmcfm_703.so";
  TString LIBJHUGENMELA="libjhugenmela.so";
  TString LIBMELA="libMELA.so";

  gInterpreter->AddIncludePath("$ROOFITSYS/include/");
  gInterpreter->AddIncludePath("../interface/");
  //////////////////////////////////////
  //these explicit loads are required on
  //some machines but not others
  //not entirely sure why
  //either way, they shouldn't hurt
  gSystem->Load("libRooFit");
  gSystem->Load("libPhysics");
  gSystem->Load("libgfortran");
  //////////////////////////////////////
  gSystem->Load(LIBMCFMPATH + LIBMCFM);
  gSystem->Load(LIBMCFMPATH + LIBJHUGENMELA);
  gSystem->Load(LIBMCFMPATH + LIBMELA);
}
