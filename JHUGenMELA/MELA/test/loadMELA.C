{
  TString loadMELA = __FILE__;
  TString testdir = loadMELA(0, loadMELA.Last('/'));
  TString LIBMCFMPATH = testdir+"/../data/slc6_amd64_gcc530/";
  TString LIBMCFM="libmcfm_705.so";
  TString LIBJHUGENMELA="libjhugenmela.so";
  TString LIBMELA="libMELA.so";
  TString LIBCOLLIER="libcollier.so";

  gInterpreter->AddIncludePath("$ROOFITSYS/include/");
  gInterpreter->AddIncludePath(testdir+"/../interface/");
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
  gSystem->Load(LIBMCFMPATH + LIBCOLLIER);
}
