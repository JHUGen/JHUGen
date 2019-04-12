#ifndef _TMCFM_HH_
#define _TMCFM_HH_

//----------------------------------------------------------
//
// http://www.chiralcomp.com/support/mixing_f77_c_cpp/
//  (almost all) Fortran compilers add, during compilation, an underscore (_) at the end of the Fortran routine names. Our experience is that the f77 compiler in HP Unix environments does not do this.
//    * We noticed that if a fortran subroutine has an underscore anywhere in its name, the GNU g77 compiler adds two (2) underscores at the end of the name.
//
//   In the compiled Fortran code all arguments to functions are passed by their address.
//   If a Fortran function takes a character string as an argument, the string length must be passed as the last argument (i.e. after the "ordinary" argument list). 
//
//----------------------------------------------------------
// MCFM parameters
enum {nSupportedHiggses=2}; // This is not actually an integer parameter in MCFM. As of migration to MCFM 7.0, this is the total number of Higgses supported, at most.
// constants.f 
enum {nf=5}; // nflavors of "pdf" quarks
// constants.f 
enum {mxpart=14}; // maxpart, max. number of particles. MCFM 6.8+: CAUTION!!! IMPORTANT TO CHECK WHEN UPDATING MCFM!!!
//enum {mxpart=12}; // maxpart, max. number of particles. MCFM 6.7
//mxdim.f
enum {ndims=22};

enum {nmsq=11};


extern "C" {
  //---------------------------------
  // Structure
  //---------------------------------
  extern struct{
    int nproc;
  } nproc_;

  extern struct{
    bool qlfirst;
  } qlfirst_;

  extern struct{
    bool Qflag, Gflag, QandGflag;
  } flags_;


#define bveg1_mcfm_ bveg1_
  extern struct{
    double xl[ndims], xu[ndims], acc;
    int ndim, ncall, itmx, nprn;
  } bveg1_mcfm_;

  extern struct {
    int ih1, ih2;
  } density_;

  extern struct{
    double scale, musq;
  } scale_;

  extern struct{
    double facscale;
  } facscale_;

  extern struct {
    int n2; int n3; double mass2; double width2; double mass3; double width3;
  } breit_;

  extern struct{
    int nqcdjets, nqcdstart;
  } nqcdjets_;


  extern struct{
    double xmin;
  } xmin_;


  extern struct{
    int npart;
  } npart_;

  extern struct{
    int nuflav;
  } nuflav_;

  extern struct{
    double vsymfact;
  } vsymfact_;

  extern struct{
    bool interference;
  } interference_;

  extern struct{
    double cutoff;
  } cutoff_;

  extern struct{
    double amz;
  } couple_;

  extern struct{
    double Gf_inp, aemmz_inp, xw_inp, wmass_inp, zmass_inp;
  } ewinput_;

  extern struct{
    int ewscheme;
  } ewscheme_;

  extern struct {
    double gsq, as, ason2pi, ason4pi;
  } qcdcouple_;

  extern struct {
    double Gf, gw, xw, gwsq, esq, vevsq;
  } ewcouple_;

  extern struct {
    double aemmz;
  } em_;

  extern struct {
    double Q[11], tau[11];
  } ewcharge_;

  extern struct {
    double delg1_z, delg1_g, lambda_g, lambda_z, delk_g, delk_z, tevscale;
  } anomcoup_;


  //mcfm/src/Inc/masses.F
#define masses_mcfm_ masses_

  extern struct{
    double md, mu, ms, mc, mb, mt,
    mel, mmu, mtau,
    hmass, hwidth,
    wmass, wwidth,
    zmass, zwidth,
    twidth,
    tauwidth,
    mtausq, mcsq, mbsq;
  } masses_mcfm_;


  extern struct{

    int AllowAnomalousCouplings;
    int distinguish_HWWcouplings;
    int AnomalCouplPR;
    int AnomalCouplDK;
    int channeltoggle_stu;
    int vvhvvtoggle_vbfvh;

    int cz_q1sq, cz_q2sq, cz_q12sq;
    int cw_q1sq, cw_q2sq, cw_q12sq;
    int c2z_q1sq, c2z_q2sq, c2z_q12sq;
    int c2w_q1sq, c2w_q2sq, c2w_q12sq;

    double mb_4gen, mt_4gen;

    double LambdaBSM, Lambda_Q;
    double Lambda_zgs1;
    double Lambda_z1, Lambda_z2, Lambda_z3, Lambda_z4;
    double Lambda_z11, Lambda_z21, Lambda_z31, Lambda_z41;
    double Lambda_z12, Lambda_z22, Lambda_z32, Lambda_z42;
    double Lambda_z10, Lambda_z20, Lambda_z30, Lambda_z40;
    double Lambda_w1, Lambda_w2, Lambda_w3, Lambda_w4;
    double Lambda_w11, Lambda_w21, Lambda_w31, Lambda_w41;
    double Lambda_w12, Lambda_w22, Lambda_w32, Lambda_w42;
    double Lambda_w10, Lambda_w20, Lambda_w30, Lambda_w40;

    double h2mass, h2width;

    double Lambda2BSM, Lambda2_Q;
    double Lambda2_zgs1;
    double Lambda2_z1, Lambda2_z2, Lambda2_z3, Lambda2_z4;
    double Lambda2_z11, Lambda2_z21, Lambda2_z31, Lambda2_z41;
    double Lambda2_z12, Lambda2_z22, Lambda2_z32, Lambda2_z42;
    double Lambda2_z10, Lambda2_z20, Lambda2_z30, Lambda2_z40;
    double Lambda2_w1, Lambda2_w2, Lambda2_w3, Lambda2_w4;
    double Lambda2_w11, Lambda2_w21, Lambda2_w31, Lambda2_w41;
    double Lambda2_w12, Lambda2_w22, Lambda2_w32, Lambda2_w42;
    double Lambda2_w10, Lambda2_w20, Lambda2_w30, Lambda2_w40;

    double kappa_top[2]; double kappa_tilde_top[2];
    double kappa_bot[2]; double kappa_tilde_bot[2];
    double ghg2[2]; double ghg3[2]; double ghg4[2];
    double kappa_4gen_top[2]; double kappa_tilde_4gen_top[2];
    double kappa_4gen_bot[2]; double kappa_tilde_4gen_bot[2];
    double ghg2_4gen[2]; double ghg3_4gen[2]; double ghg4_4gen[2];

    double ghz1[2]; double ghz2[2]; double ghz3[2]; double ghz4[2]; // No additional q2 dependence
    double ghz1_prime[2]; double ghz2_prime[2]; double ghz3_prime[2]; double ghz4_prime[2]; // Dipole ansatz
    double ghz1_prime2[2]; double ghz2_prime2[2]; double ghz3_prime2[2]; double ghz4_prime2[2]; // |q1**2| + |q2**2|
    double ghz1_prime3[2]; double ghz2_prime3[2]; double ghz3_prime3[2]; double ghz4_prime3[2]; // |q1**2| - |q2**2|
    double ghz1_prime4[2]; double ghz2_prime4[2]; double ghz3_prime4[2]; double ghz4_prime4[2]; // (q1 + q2)**2
    double ghz1_prime5[2]; double ghz2_prime5[2]; double ghz3_prime5[2]; double ghz4_prime5[2]; // q1**4 + q2**4
    double ghz1_prime6[2]; double ghz2_prime6[2]; double ghz3_prime6[2]; double ghz4_prime6[2]; // q1**4 - q2**4
    double ghz1_prime7[2]; double ghz2_prime7[2]; double ghz3_prime7[2]; double ghz4_prime7[2]; // |q1**2| * |q2**2|

    double ghzgs1_prime2[2]; double ghzgs2[2]; double ghzgs3[2]; double ghzgs4[2];
    double ghgsgs2[2]; double ghgsgs3[2]; double ghgsgs4[2];

    double ghw1[2]; double ghw2[2]; double ghw3[2]; double ghw4[2]; // No additional q2 dependence
    double ghw1_prime[2]; double ghw2_prime[2]; double ghw3_prime[2]; double ghw4_prime[2]; // Dipole ansatz
    double ghw1_prime2[2]; double ghw2_prime2[2]; double ghw3_prime2[2]; double ghw4_prime2[2]; // |q1**2| + |q2**2|
    double ghw1_prime3[2]; double ghw2_prime3[2]; double ghw3_prime3[2]; double ghw4_prime3[2]; // |q1**2| - |q2**2|
    double ghw1_prime4[2]; double ghw2_prime4[2]; double ghw3_prime4[2]; double ghw4_prime4[2]; // (q1 + q2)**2
    double ghw1_prime5[2]; double ghw2_prime5[2]; double ghw3_prime5[2]; double ghw4_prime5[2]; // q1**4 + q2**4
    double ghw1_prime6[2]; double ghw2_prime6[2]; double ghw3_prime6[2]; double ghw4_prime6[2]; // q1**4 - q2**4
    double ghw1_prime7[2]; double ghw2_prime7[2]; double ghw3_prime7[2]; double ghw4_prime7[2]; // |q1**2| * |q2**2|


    double kappa2_top[2]; double kappa2_tilde_top[2];
    double kappa2_bot[2]; double kappa2_tilde_bot[2];
    double gh2g2[2]; double gh2g3[2]; double gh2g4[2];
    double kappa2_4gen_top[2]; double kappa2_tilde_4gen_top[2];
    double kappa2_4gen_bot[2]; double kappa2_tilde_4gen_bot[2];
    double gh2g2_4gen[2]; double gh2g3_4gen[2]; double gh2g4_4gen[2];

    double gh2z1[2]; double gh2z2[2]; double gh2z3[2]; double gh2z4[2]; // No additional q2 dependence
    double gh2z1_prime[2]; double gh2z2_prime[2]; double gh2z3_prime[2]; double gh2z4_prime[2]; // Dipole ansatz
    double gh2z1_prime2[2]; double gh2z2_prime2[2]; double gh2z3_prime2[2]; double gh2z4_prime2[2]; // |q1**2| + |q2**2|
    double gh2z1_prime3[2]; double gh2z2_prime3[2]; double gh2z3_prime3[2]; double gh2z4_prime3[2]; // |q1**2| - |q2**2|
    double gh2z1_prime4[2]; double gh2z2_prime4[2]; double gh2z3_prime4[2]; double gh2z4_prime4[2]; // (q1 + q2)**2
    double gh2z1_prime5[2]; double gh2z2_prime5[2]; double gh2z3_prime5[2]; double gh2z4_prime5[2]; // q1**4 + q2**4
    double gh2z1_prime6[2]; double gh2z2_prime6[2]; double gh2z3_prime6[2]; double gh2z4_prime6[2]; // q1**4 - q2**4
    double gh2z1_prime7[2]; double gh2z2_prime7[2]; double gh2z3_prime7[2]; double gh2z4_prime7[2]; // |q1**2| * |q2**2|

    double gh2zgs1_prime2[2]; double gh2zgs2[2]; double gh2zgs3[2]; double gh2zgs4[2];
    double gh2gsgs2[2]; double gh2gsgs3[2]; double gh2gsgs4[2];

    double gh2w1[2]; double gh2w2[2]; double gh2w3[2]; double gh2w4[2]; // No additional q2 dependence
    double gh2w1_prime[2]; double gh2w2_prime[2]; double gh2w3_prime[2]; double gh2w4_prime[2]; // Dipole ansatz
    double gh2w1_prime2[2]; double gh2w2_prime2[2]; double gh2w3_prime2[2]; double gh2w4_prime2[2]; // |q1**2| + |q2**2|
    double gh2w1_prime3[2]; double gh2w2_prime3[2]; double gh2w3_prime3[2]; double gh2w4_prime3[2]; // |q1**2| - |q2**2|
    double gh2w1_prime4[2]; double gh2w2_prime4[2]; double gh2w3_prime4[2]; double gh2w4_prime4[2]; // (q1 + q2)**2
    double gh2w1_prime5[2]; double gh2w2_prime5[2]; double gh2w3_prime5[2]; double gh2w4_prime5[2]; // q1**4 + q2**4
    double gh2w1_prime6[2]; double gh2w2_prime6[2]; double gh2w3_prime6[2]; double gh2w4_prime6[2]; // q1**4 - q2**4
    double gh2w1_prime7[2]; double gh2w2_prime7[2]; double gh2w3_prime7[2]; double gh2w4_prime7[2]; // |q1**2| * |q2**2|

    double dV_A[2]; double dP_A[2]; double dM_A[2]; double dFour_A[2];
    double dV_Z[2]; double dP_Z[2]; double dM_Z[2]; double dFour_Z[2];
    double dZZWpWm[2]; double dZAWpWm[2]; double dAAWpWm[2];

  } spinzerohiggs_anomcoupl_;

  extern struct {
    bool srdiags;
  } srdiags_;

  //mcfm/src/Inc/noglue.f
  extern struct{
    bool noglue, ggonly, gqonly, omitgg;
  } noglue_;

  //mcfm/src/Inc/zcouple.F
  extern struct{
    double l[nf], r[nf], q1, l1, r1, q2, l2, r2, le, ln, re, rn, sin2w;
  } zcouple_;

  extern struct{
    int nwz;
  } nwz_;

  extern struct {
    double Vud, Vus, Vub, Vcd, Vcs, Vcb;
  } cabib_;
  
  extern struct {
    double taumin;
  } taumin_;

  extern struct {
    double sqrts;
  } energy_;

  extern struct {
    int nlooprun;
  } nlooprun_;

  extern struct {
    int nflav;
  } nflav_;

  extern struct {
    int lastphot;
  } lastphot_;

  extern struct {
    char runstring[30];
  } runstring_;

  extern struct {
    char pdlabel[7];
  } pdlabel_;

  extern struct {
    char plabel[mxpart][2];
  } plabel_;


  //---------------------------------
  // function
  //---------------------------------


  //##############
  // Initialization
  //##############
  void mcfm_init_(char* inputfile, char* workdir);
  void chooser_();
  void coupling_();
  void coupling2_();
  void qlinit_();
  void fdist_(int* ih1, double* xx, double* pdfscale, double* fx1);
  double alphas_(double* q, double* amz, int* nloop);
  void ckmfill_(int* nwz);


  //###############

  // WW
#define qqb_ww_ qqb_ww_
  void qqb_ww_(double* p, double* msq);

  // WZ
#define qqb_wz_ qqb_wz_
  void qqb_wz_(double* p, double* msq);

  // ZZ->4l
#define qqb_zz_ qqb_zz_
  void qqb_zz_(double* p, double* msq);

#define qqb_zz_stu_ qqb_zz_stu_ // Custom qqb->ZZ for different s, t, u channels
  void qqb_zz_stu_(double* p, double* msq, int* channeltoggle);

#define qq_zzqq_ qq_zzqq_ // WBF-ZZ
  void qq_zzqq_(double* p, double* msq);
//#define qq_zzqq_bkg_ qq_zzqq_bkg_ // WBF-ZZ
//  void qq_zzqq_bkg_(double* p, double* msq);
#define qq_wwqq_ qq_wwqq_ // WBF-WW
  void qq_wwqq_(double* p, double* msq);
#define qq_vvqq_ qq_vvqq_ // WBF-VV
  void qq_vvqq_(double* p, double* msq);

#define qq_zzqqstrong_ qq_zzqqstrong_ // JJQCD-ZZ
  void qq_zzqqstrong_(double* p, double* msq);
#define qq_wwqqstrong_ qq_wwqqstrong_ // JJQCD-WW
  void qq_wwqqstrong_(double* p, double* msq);

#define gg_hzz_tb_ gg_hzz_tb_
  void gg_hzz_tb_(double* p, double* msq);
#define gg_zz_hpi_ gg_zz_hpi_ // Only H+interf, no |gg->ZZ|**2
  void gg_zz_hpi_(double* p, double* msq);
#define gg_zz_all_ gg_zz_all_
  void gg_zz_all_(double* p, double* msq);
#define gg_zz_ gg_zz_
  void gg_zz_(double* p, double* msq);
#define gg_zz_int_ gg_zz_int_
  void gg_zz_int_(double* p, double* msq);


  // WW+ZZ, p as in ZZ
#define gg_hvv_tb_ gg_hvv_tb_
  void gg_hvv_tb_(double* p, double* msq);
#define gg_vv_all_ gg_vv_all_
  void gg_vv_all_(double* p, double* msq);
#define gg_vv_ gg_vv_
  void gg_vv_(double* p, double* msq);

  // ZZ->2l2q
#define qqb_z2jet_ qqb_z2jet_
  void qqb_z2jet_(double* p, double* msq);

  // ZG->2lG
#define qqb_zgam_ qqb_zgam_
  void qqb_zgam_(double* p, double* msq);

#define qqb_hww_ qqb_hww_
  void qqb_hww_(double* p, double* msq);

#define qqb_hzz_ qqb_hzz_
  void qqb_hzz_(double* p, double* msq);

#define qqb_hzz_tb_ qqb_hzz_tb_
  void qqb_hzz_tb_(double* p, double* msq);

#define qqb_z_ qqb_z_
  void qqb_z_(double* p, double *msq);

#define qqb_wgam_ qqb_wgam_
  void qqb_wgam_(double* p, double* msq);

#define qqb_wgam_ qqb_wgam_
  void qqb_wgam_(double* p, double* msq);

  // qqb/gg->GG
#define qqb_gamgam_ qqb_gamgam_
  void qqb_gamgam_(double* p, double* msq);


#define qqb_w_g_ qqb_w_g_
  void qqb_w_g_(double* p, double* msq);

  //###############

  double lowint_(double* r, double* wgt);
  void dotem_(int N, double* p, double* s);
  void boost_mcfm_(double* mass, double* p1, double* p_in, double* p_out);
  void breitw_(double* x1, double* mminsq, double* mmaxsq, double* rmass, double* rwidth, double* msq, double* wt);
  void phi3m0_(double* xth, double* xphi, double* p0, double* p1, double* p2, double* wt);
  void gen2_(double* r, double* p, double* pswt);
  void gen3_(double* r, double* p, double* pswt);
  void gen3jet_(double* r, double* p, double* pswt);
  void gen4_(double* r, double* p, double* wt4);
  void gen4h_(double* r, double* p, double* wt4);
  void gen_njets_(double* r, int*, double* p, double* msq);

}

#endif
