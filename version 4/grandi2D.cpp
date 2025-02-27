// Original model:
// A novel computational model of the human ventricular action potential and Ca transient.
// Grandi E, Pasqualini FS, Bers DM. J Mol Cell Cardiol. 2010 Jan;48(1):112-21. Epub 2009 Oct 14. PMID:19835882
//
// Heart failure model:
// Simulation and Mechanistic Investigation of the Arrhythmogenic Role of the Late Sodium Current in Human Heart Failure.
// Trenor, B., Cardona, K., Gomez, J. F., Rajamani, S., Ferrero, J. M., Belardinelli, L., & Saiz, J. (2012). (M. Rota, Ed.)PLoS ONE, 7(3), e32659. doi:10.1371/journal.pone.0032659.s015
//
// This version was converted from original matlab code to C++ by:
// Author: Tashalee Brown (July 2013)
// v1.0: Non-variable time step implementation
// v2.0: Added APD and CV calculations
// v3.0: Adding vogel gap junctions
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include <time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
//must add flags -I /opt/local/include
using namespace std;


//Model characteristics
const double epi = 0;//Cell type of model is originally endo
const double HF = 0; //To simulate the heart failure model use HF=1
const bool DynamicGapON = true;

//Global stimulation parameters
int num_stim=5;
int stim_equil = 1;
int file_filter = 1000;
double stim_mag=19.5;
double stim_dur=3;
int stim_col_max=5, stim_row_max=5;
double bcl=500;
double DT = 0.005;
char output_file_name[30] = "testGrandi3.dat";
char output_file_name2[30] = "testGrandi4.dat";
char output_file_name3[30] = "params.dat";
//char output_file_name4[30] = "voltage.dat";

const char *neighbors_myo_file = "neighbors_myo_file.dat";
const char *neighbors_fib_file = "neighbors_fib_file.dat";
const char *grid_file = "grid_file.dat";

//Declaration of scaling parameters
double scaleNaL;
double scalehL;
double scaleto;
double scaleK1;
double scaleNaK;
double scaleNab;
double scaleCab;
double scaleNCX;
double scaleSERCA;
double scaleleak;
double ec50SR;


//Function to calculate PDE
void PDE(vector<vector<double> > &, vector<vector<double> > &, vector <vector<int> >, bool);



//Function to read files and retune 2D vector of content
vector<vector<int> > fileToVector(const char *name);

enum cx_model{
  CX43CX43,
  CX43CX45,
  CX45CX45,
};

//Function for dynamic gap junctions
double dyngap(double, cx_model);
//Function for dynamic gap junctions
double dyngap1(double, cx_model);
//Function for dynamic gap junctions
double dyngap2(double, cx_model);

  struct voltage_params{
    double G1H, G2H, G1L, G2L, V1H, V2H, V1L, V2L, VJ;
  };
//Header for root solving
const gsl_root_fdfsolver_type *TypePtr;
gsl_root_fdfsolver *SolverPtr;

gsl_function_fdf FDF_hh;//defines gsl functions
gsl_function_fdf FDF_ll;
gsl_function_fdf FDF_hl;
gsl_function_fdf FDF_lh;

double rootSolver(gsl_function_fdf, double);
//Functions to be solved for roots, each with fcn, deriv, and fdf (which is solving fcn and deriv simultaneously)
double voltagehh (double vhh, void *params);
double voltagehh_deriv (double vhh, void *params);
void voltagehh_fdf (double vhh, void *params, double *_vhh, double *_dvhh);
double voltagell (double vll, void *params);
double voltagell_deriv (double vll, void *params);
void voltagell_fdf (double vll, void *params, double *_vll, double *_dvll);
double voltagehl (double vhl, void *params);
double voltagehl_deriv (double vhl, void *params);
void voltagehl_fdf (double vhl, void *params, double *_vhl, double *_dvhl);
double voltagelh (double vlh, void *params);
double voltagelh_deriv (double vlh, void *params);
void voltagelh_fdf (double vlh, void *params, double *_vlh, double *_dvlh);

//Function prototypes
//void calcAPD(double, double, int);

/////////////////////////////////////////////////////////////////////////////////////////////
int main(){

  //Model Initialization
  //Initialize Root solver
  TypePtr = gsl_root_fdfsolver_newton;//Creates an instance of newton solver
  SolverPtr = gsl_root_fdfsolver_alloc (TypePtr);

  cout <<"\n Initialized GSL Root Solver: "<<gsl_root_fdfsolver_name(SolverPtr)<<endl;

  clock_t t1, t2;
  t1=clock();
  //Scaling parameters for HF
  if (HF == 1){
    scaleNaL=2.;
    scalehL=2.;
    scaleto=.4;
    scaleK1=.68;
    scaleNaK=.9;
    scaleNab=0;
    scaleCab=1.53;
    scaleNCX=2.;
    scaleSERCA=.5;
    scaleleak=3.;
  }
  else {
    scaleNaL=1.;
    scalehL=1.;
    scaleto=1.;
    scaleK1=1.;
    scaleNaK=1.;
    scaleNab=1.;
    scaleCab=1.;
    scaleNCX=1.;
    scaleSERCA=1.;
    scaleleak=1.;
  }

  FILE *output, *output2, *output3;//, *output4;
  output = fopen(output_file_name, "w");
  output2 = fopen(output_file_name2, "w");
  output3 = fopen(output_file_name3, "w");
  // output4 = fopen(output_file_name4, "w");

  vector <vector<int> > neighbors_myo = fileToVector(neighbors_myo_file);
  vector <vector<int> > neighbors_fib = fileToVector(neighbors_fib_file);
  //////////////////////////////////////////////////////////////////////////////////////////
  vector <vector<int> > thegrid = fileToVector(grid_file);
  const int NCOLS = thegrid.size();
  const int NROWS = thegrid[0].size();
  const int stim_length = thegrid.size()*thegrid[0].size();
  int stim_myo[stim_length];
  int stim_fib[stim_length];//split into stim myo and stim fib array  
  for (int k=0; k < stim_length; k++){
    stim_myo[k]=0;
    stim_fib[k]=0;
  }

  int prev=0, curr=0;
  int stim_fib_index = 0;
  int stim_myo_index = 0;

  for (int r=0; r < stim_row_max; r++){
    for (int c=0; c < stim_col_max; c++){
      curr = thegrid[r][c];
      if (curr > 0){//if myocyte
        if (curr != prev){
          cout << "Stimulate myocyte #: " << curr <<endl;
          stim_myo[stim_myo_index] = 1;
          stim_myo_index++;
        }
      }
      else{//if fibroblast
        if (curr != prev){
          cout << "Stimulate non-myocyte #: " << curr <<endl;
          stim_fib[stim_fib_index] = 1;
          stim_fib_index++;
        }
      }
      prev = curr;
    }
  }

  const int CVL_cell1 = floor(NCOLS/4);
  const int CVL_cell2 = floor((NCOLS/4)*3);
  const int CVT_cell1 = floor(NROWS/4);
  const int CVT_cell2 = floor((NROWS/4)*3);

  const int CVL_1 = thegrid[floor(NROWS/2)][CVL_cell1];
  const int CVL_2 = thegrid[floor(NROWS/2)][CVL_cell2];
  const int CVT_1 = thegrid[CVT_cell1][floor(NROWS/2)];
  const int CVT_2 = thegrid[CVT_cell2][floor(NROWS/2)];
  double CVstart;

  cout <<"Longitudinal cells: "<< CVL_1 << " & " << CVL_2<<endl;
  cout <<"Transverse cells: "<< CVT_1 << " & " << CVT_2 <<endl;

  ///////////////////////////////////////////////////////////////
  //Number of myocyte cells and fibroblast cells
  const int CELLS_myo = neighbors_myo.size();
  int CELLS_fib ;
  if (!neighbors_fib.empty()){
    CELLS_fib = neighbors_fib.size();
  }
  else{
    CELLS_fib = 0;
  }
  cout << "Fibroblasts: " << CELLS_fib << endl;
  cout << "Myocytes: " << CELLS_myo << endl;

  //Variables to measure AP parameters
  double V_myo[CELLS_myo];
  double V_fib[CELLS_fib];
  vector<double> Vrest;
  vector<double> apd90;
  vector<int> apd90_start;
  vector<double> apd90_t1;
  vector<double> apd90_t2;
  vector<double> upstroke90;
  vector<double> downstroke90;
  Vrest.resize(CELLS_myo);
  apd90.resize(CELLS_myo);
  apd90_start.resize(CELLS_myo);
  apd90_t1.resize(CELLS_myo);
  apd90_t2.resize(CELLS_myo);
  upstroke90.resize(CELLS_myo);
  downstroke90.resize(CELLS_myo);

  //Constants
  const double R = 8314;       // [J/kmol*K]
  const double Frdy = 96485;   // [C/mol]
  const double Temp = 310;     // [K]
  const double FoRT = Frdy/R/Temp;
  const double Cmem = 1.3810e-10;   // [F] membrane capacitance
  const double Qpow = (Temp-310)/10;
  const double pi = 3.141592653589793;

  //Cell geometry (consider simplifying calculations
  double cellLength = 100;     // cell length [um]
  double cellRadius = 10.25;   // cell radius [um]
  double junctionLength = 160e-3;  // junc length [um]
  double junctionRadius = 15e-3;   // junc radius [um]
  double distSLcyto = 0.45;    // dist. SL to cytosol [um]
  double distJuncSL = 0.5;  // dist. junc to SL [um]
  double DcaJuncSL = 1.64e-6;  // Dca junc to SL [cm^2/sec]
  double DcaSLcyto = 1.22e-6; // Dca SL to cyto [cm^2/sec]
  double DnaJuncSL = 1.09e-5;  // Dna junc to SL [cm^2/sec]
  double DnaSLcyto = 1.79e-5;  // Dna SL to cyto [cm^2/sec]
  double Vcell = pi*cellRadius*cellRadius*cellLength*1e-15;    // [L]
  double Vmyo = 0.65*Vcell;
  double Vsr = 0.035*Vcell;
  double Vsl = 0.02*Vcell;
  double Vjunc = 0.0539*.01*Vcell;
  double SAjunc = 20150*pi*2*junctionLength*junctionRadius;  // [um^2]
  double SAsl = pi*2*cellRadius*cellLength;          // [um^2]
  //Omitted some commented lines from original code
  double J_ca_juncsl =1/1.2134e12; // [L/msec] = 8.2413e-13
  double J_ca_slmyo = 1/2.68510e11; // [L/msec] = 3.2743e-12
  double J_na_juncsl = 1/(1.6382e12/3*100); // [L/msec] = 6.1043e-13
  double J_na_slmyo = 1/(1.8308e10/3*100);  // [L/msec] = 5.4621e-11

  //Fractional currents in compartments
  double Fjunc = 0.11;
  double Fsl = 1-Fjunc;
  double Fjunc_CaL = 0.9;
  double Fsl_CaL = 1-Fjunc_CaL;

  //Fixed ion concentrations
  double Cli = 15;   // Intracellular Cl  [mM]
  double Clo = 150;  // Extracellular Cl  [mM]
  double Ko = 5.4;   // Extracellular K   [mM]
  double Nao = 140;  // Extracellular Na  [mM]
  double Cao = 1.8;  // Extracellular Ca  [mM]
  double Mgi = 1;    // Intracellular Mg  [mM]  <------------------??

  //Nerst potential (non changing)
  double ecl;            // [mV]
  double ena_junc;     // [mV]
  double ena_sl;       // [mV]
  //Some commented code was removed from original code
  double ek;	        // [mV]
  double eca_junc;   // [mV]
  double eca_sl;     // [mV]

  //Na transport parameters
  double GNa=23;
  double GNaL=0.015;// ---------------------------CONDUCTANCE OF INAL FROM HUMAN
  double GNaB = 0.597e-3;    // [mS/uF] 0.897e-3
  double IbarNaK = 1.0*1.8;//1.90719;     // [uA/uF]
  double KmNaip = 11;         // [mM]11
  double KmKo =1.5;         // [mM]1.5
  double Q10NaK = 1.63;
  double Q10KmNai = 1.39;

  //// K current parameters
  double pNaK = 0.01833;
  double gkp = 2*0.001;

  // Cl current parameters
  double GClCa =0.5* 0.109625;   // [mS/uF]
  double GClB = 1*9e-3;        // [mS/uF]
  double KdClCa = 100e-3;    // [mM]

  // I_Ca parameters
  double pNa = 0.50*1.5e-8;       // [cm/sec]
  double pCa = 0.50*5.4e-4;       // [cm/sec]
  double pK = 0.50*2.7e-7;        // [cm/sec]
  double Q10CaL = 1.8;

  //// Ca transport parameters
  double IbarNCX = 1.0*4.5;      // [uA/uF]5.5 before - 9 in rabbit
  double KmCai = 3.59e-3;    // [mM]
  double KmCao = 1.3;        // [mM]
  double KmNai = 12.29;      // [mM]
  double KmNao = 87.5;       // [mM]
  double ksat = 0.32;        // [none]
  double nu = 0.27;          // [none]
  double Kdact = 0.150e-3;   // [mM]
  double Q10NCX = 1.57;      // [none]
  double IbarSLCaP = 0.0673; // IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
  double KmPCa = 0.5e-3;     // [mM]
  double GCaB = 5.513e-4;    // [uA/uF] 3
  double Q10SLCaP = 2.35;    // [none]

  // SR flux parameters
  double Q10SRCaP = 2.6;          // [none]
  double Vmax_SRCaP = 1.0*5.3114e-3;  // [mM/msec] (286 umol/L cytosol/sec)
  double Kmf = 0.246e-3;          // [mM] default
  //double //Kmf = 0.175e-3;          // [mM]
  double Kmr = 1.7;               // [mM]L cytosol
  double hillSRCaP = 1.787;       // [mM]
  double ks = 25;                 // [1/ms]
  double koCa = 10;               // [mM^-2 1/ms]   //default 10   modified 20
  double kom = 0.06;              // [1/ms]
  double kiCa = 0.5;              // [1/mM/ms]
  double kim = 0.005;             // [1/ms]

  if (HF == 1){
    ec50SR = 0.4;           // [mM]
  }
  else{
    ec50SR = 0.45;           // [mM]  in HF
  }

  // Buffering parameters
  // koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
  double Bmax_Naj = 7.561;       // [mM] // Bmax_Naj = 3.7; (c-code difference?)  // Na buffering
  double Bmax_Nasl = 1.65;       // [mM]
  double koff_na = 1e-3;         // [1/ms]
  double kon_na = 0.1e-3;        // [1/mM/ms]
  double Bmax_TnClow = 70e-3;    // [mM]                      // TnC low affinity
  double koff_tncl = 19.6e-3;    // [1/ms]
  double kon_tncl = 32.7;        // [1/mM/ms]
  double Bmax_TnChigh = 140e-3;  // [mM]                      // TnC high affinity
  double koff_tnchca = 0.032e-3; // [1/ms]
  double kon_tnchca = 2.37;      // [1/mM/ms]
  double koff_tnchmg = 3.33e-3;  // [1/ms]
  double kon_tnchmg = 3e-3;      // [1/mM/ms]
  double Bmax_CaM = 24e-3;       // [mM] **? about setting to 0 in c-code**   // CaM buffering
  double koff_cam = 238e-3;      // [1/ms]
  double kon_cam = 34;           // [1/mM/ms]
  double Bmax_myosin = 140e-3;   // [mM]                      // Myosin buffering
  double koff_myoca = 0.46e-3;   // [1/ms]
  double kon_myoca = 13.8;       // [1/mM/ms]
  double koff_myomg = 0.057e-3;  // [1/ms]
  double kon_myomg = 0.0157;     // [1/mM/ms]
  double Bmax_SR = 19*.9e-3;     // [mM] (Bers text says 47e-3) 19e-3
  double koff_sr = 60e-3;        // [1/ms]
  double kon_sr = 100;           // [1/mM/ms]
  double Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl;        // [mM]    // SL buffering
  double Bmax_SLlowj = 4.6e-3*Vmyo/Vjunc*0.1;    // [mM]    //Fei *0.1!!! junction reduction factor
  double koff_sll = 1300e-3;     // [1/ms]
  double kon_sll = 100;          // [1/mM/ms]
  double Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl;       // [mM]
  double Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1;  // [mM] //Fei *0.1!!! junction reduction factor
  double koff_slh = 30e-3;       // [1/ms]
  double kon_slh = 100;          // [1/mM/ms]
  double Bmax_Csqn = 140e-3*Vmyo/Vsr;            // [mM] // Bmax_Csqn = 2.6;      // Csqn buffering
  double koff_csqn = 65;         // [1/ms]
  double kon_csqn = 100;         // [1/mM/ms]

  //Some commented code on membrane currents emitted

  //Gating variables
  //I_Na: Fast Na Current
  double mss, taum, ah, bh, tauh, hss, aj, bj, tauj, jss;
  double aml, bml, mlss, tauml, hlss, tauhl;
  double sigma, fnak;
  double xrss, tauxr, rkr, xsss, tauxs, xtoss, tauxtof, ytoss, tauytof;
  double xkurss, tauxkur, ykurss, tauykur, aki, bki;
  double dss, taud, fss, tauf;
  double I_Na_junc, I_Na_sl, I_Na, I_NaL_junc, I_NaL_sl, I_NaL, I_nabk_junc, I_nabk_sl, I_nabk;
  double I_nak_junc, I_nak_sl, I_nak;
  double I_kr, I_ks, I_kp, I_kur, I_ki;
  double I_ClCa_junc, I_ClCa_sl, I_ClCa, I_Clbk;
  double ibarca_j, ibarca_sl, ibark, ibarna_j, ibarna_sl, I_Ca_junc, I_Ca_sl, I_Ca, CaSR;
  double I_CaK, I_CaNa_junc, I_CaNa_sl, I_CaNa, I_Catot;
  //Added variables
  double kiSRCa, koSRCa, MinSR, MaxSR, fcaCaj, fcaCaMSL, kiss, GtoFast, GtoSlow, I_tos;
  double I_tof, tauytos, tauxtos, I_kp_sl, I_kp_junc, kp_kp, I_ks_sl,eks, I_ks_junc, gks_sl, gks_junc, gkr;

  double Ka_junc, Ka_sl, s1_junc, s1_sl, s2_junc, s2_sl, s3_junc, s3_sl, I_ncx_junc, I_ncx_sl, I_ncx;
  double I_pca_junc, I_pca_sl, I_pca;
  double I_cabk_junc, I_cabk_sl, I_cabk;
  double I_Na_tot, I_Cl_tot, I_Ca_tot, I_tot, I_stim;
  double kCaSR, RI, J_SRCarel, J_serca, J_SRleak, J_CaB_cytosol;
  double J_CaB_junction, J_CaB_sl;
  double I_Na_tot_junc, I_Na_tot_sl, I_Na_tot_sl2, I_Na_tot_junc2;
  double I_K_tot, I_Ca_tot_junc, I_Ca_tot_sl, I_to;
  float r;

  //Fibroblast parameters
  //---------------------------------------------------------------------------------------------------

  double C_fib  = 6.3; //cell capacitance in pF
  double GKv = 0.25; //The maximal conductance of IKv_fib in nS/pF or uS/nF
  double GK1_fib = 0.4822; //The maximal conductance of this inward rectifying K+ current in nS/pF of uS/nF
  double B     = -200; //The empirically determined constant
  double V_rev = -150; //mV The reversal potential of the electrogenic pump

  //K_mK and KmNa are binding constants
  double KmK_fib  = 1.0; //mM
  double KmNa_fib = 11.0; //mM
  double Max_INaK = 2.002; //The maximum current generated by INaK_fib in pA/pF or nA/nF
  double GbNa_fib = 0.0095; //The leak conductance in nS/pF of uS/nF

  // Initial parameters
  double Nai_fib  = 8.0; //CellML from Dr. Fink
  double ENa_fib = 70.0; //millivolts
  double EK_fib = -87.0; //millivolts
  double s_inf_fib, tau_s_fib;
  double r_inf_fib, tau_r_fib;
  double alpha_K1_fib, beta_K1_fib;
  double IK1_fib, IKv_fib, INaK_fib, IbNa_fib, Itotal_fib, I_stim_fib;

  //----------------------------------------------------------------------------------------------------
  vector< vector< double > > y(CELLS_myo, vector< double >( 42 ) );
  vector< vector< double > > ydot(CELLS_myo, vector< double >( 42 ) );
  vector< vector< double > > yf(CELLS_fib, vector< double >( 3 ) );
  vector< vector< double > > yfdot(CELLS_fib, vector< double >( 3 ) );
  //cout << "vectors loaded"<<endl;

  //Initialize myocyte variables
  for (int node=0; node<CELLS_myo; node++){
    V_myo[node]=-8.1386313e+1;//previous myo voltage
    Vrest[node]=0;
    upstroke90[node]=0;
    downstroke90[node]=0;
    apd90_start[node]=-1;
    apd90_t1[node]=0;
    apd90_t2[node]=0;
    apd90[node]=0;
    y[node][0]=-8.1386313e+1;//Duplicate of membrane potential
    y[node][1]=3.8467219e-3;//m
    y[node][2]=6.2384201e-1;//h
    y[node][3]=6.2163026e-1;//j
    y[node][4]=2.9571680e-6;//d
    y[node][5]=9.9509556e-1;
    y[node][6]=2.5423478e-2;
    y[node][7]=1.5493027e-2;
    y[node][8]=4.4273905e-4;
    y[node][9]=7.9406731e-1;
    y[node][10]=4.4273111e-4;
    y[node][11]=9.9999579e-1;
    y[node][12]=2.8804776e-2;
    y[node][13]=4.3004119e-3;
    y[node][14]=8.9309469e-1;
    y[node][15]=8.4202753e-7;
    y[node][16]=1.0079157e-7;
    y[node][17]=3.4239645e+0;
    y[node][18]=7.4707318e-1;
    y[node][19]=9.3020592e-3;
    y[node][20]=1.1949450e-1;
    y[node][21]=9.6153876e-3;
    y[node][22]=3.0933180e-4;
    y[node][23]=2.1350328e-3;
    y[node][24]=1.3734533e-1;
    y[node][25]=2.2626469e-3;
    y[node][26]=7.6307679e-3;
    y[node][27]=1.0040086e-2;
    y[node][28]=7.5004475e-2;
    y[node][29]=1.1576060e-1;
    y[node][30]=1.2287127e+0;
    y[node][31]=5.8241866e-1;
    y[node][32]=8.2712566e+0;
    y[node][33]=8.2702314e+0;
    y[node][34]=8.2704125e+0;
    y[node][35]=120;
    y[node][36]=1.8131347e-4;
    y[node][37]=1.0826539e-4;
    y[node][38]=9.1363224e-5;
    y[node][39]=-8.1386313e+1;
    y[node][40]=2.8110845e-3;
    y[node][41]=1.6274168e-1;
  }
  //Initialize fib variables
  for (int node = 0; node<CELLS_fib; node++){
    yf[node][0] = -49.21; //CellML from Dr. Fink
    yf[node][1] = 0.0; //r gate
    yf[node][2] = 1.0; //s gate
  }
  //Local pacing & stimulus parameters
  int count, bcl_int, stim_dur_int;
  double time;
  int num = 0;
  stim_dur_int   = stim_dur/DT;
  bcl_int = bcl/DT;



  while ( num < num_stim )	{
    num = num + 1;//Number is the num of beats.
    count = 0;//Integer representation of time

    for (int node = 0; node<CELLS_myo; node++){
      // Values for APD Calculation
      fprintf(output3,"%.12f\t", apd90[node]);
      apd90_start[node] = -1;
      Vrest[node] = y[node][0];
      upstroke90[node] = Vrest[node] - (Vrest[node]*0.5);
      downstroke90[node] = Vrest[node] - (Vrest[node]*0.1);
    }
    fprintf(output3,"\n");

    while ( count < bcl_int ){
      count = count + 1;
      time = count * DT;//Time since start of cycle in ms


      //Start of the myocyte dimensions
      for (int x = 0; x<CELLS_myo; x++){

        //Nerst potential without time dependence
        ecl = (1/FoRT)*log(Cli/Clo);            // [mV]
        ek = (1/FoRT)*log(Ko/y[x][35]);	        // [mV]

        //Nerst potentials
        ena_junc = (1/FoRT)*log(Nao/y[x][32]);     // [mV]
        ena_sl = (1/FoRT)*log(Nao/y[x][33]);       // [mV]
        //Some commented code was removed from original code
        eca_junc = (1/FoRT/2)*log(Cao/y[x][36]);   // [mV]
        eca_sl = (1/FoRT/2)*log(Cao/y[x][37]);     // [mV]
        eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y[x][35]+pNaK*y[x][34]));

        //I_Na: Fast Na Current
        mss = 1 / pow((1 + exp( -(56.86 + y[x][0]) / 9.03 )),2);
        taum = 0.1292 * exp(-pow(((y[x][0]+45.79)/15.54),2)) + 0.06487 * exp(-pow(((y[x][0]-4.823)/51.12),2));

        if (y[x][0] >= -40.) {
          ah = 0;
          bh = (0.77 / (0.13*(1 + exp( -(y[x][0] + 10.66) / 11.1 ))));
          aj = 0;
          bj = ((0.6 * exp( 0.057 * y[x][0])) / (1 + exp( -0.1 * (y[x][0] + 32) )));
        }
        else{
          ah = (0.057 * exp( -(y[x][0] + 80) / 6.8 ));
          bh = ((2.7 * exp( 0.079 * y[x][0]) + 3.1e5 * exp(0.3485 * y[x][0])));
          aj = (((-2.5428 * 1.0e4*exp(0.2444*y[x][0]) - 6.948e-6 * exp(-0.04391*y[x][0])) * (y[x][0] + 37.78)) / (1 + exp( 0.311 * (y[x][0] + 79.23) )));
          bj = ((0.02424 * exp( -0.01052 * y[x][0] )) / (1 + exp( -0.1378 * (y[x][0] + 40.14) )));
        }
        tauh = 1 / (ah + bh);
        hss = 1 / pow((1 + exp( (y[x][0] + 71.55)/7.43 )),2);
        tauj = 1 / (aj + bj);
        jss = 1 / ((1 + exp( (y[x][0] + 71.55)/7.43 ))*(1 + exp( (y[x][0] + 71.55)/7.43 )));

        ydot[x][1] = (mss - y[x][1]) / taum;
        ydot[x][2] = (hss - y[x][2]) / tauh;
        ydot[x][3] = (jss - y[x][3]) / tauj;

        I_Na_junc = Fjunc*GNa*y[x][1]*y[x][1]*y[x][1]*y[x][2]*y[x][3]*(y[x][0]-ena_junc);
        I_Na_sl = Fsl*GNa*y[x][1]*y[x][1]*y[x][1]*y[x][2]*y[x][3]*(y[x][0]-ena_sl);
        I_Na = I_Na_junc+I_Na_sl;

        // Late I_na
        aml=0.32*(y[x][0]+47.13)/(1-exp(-0.1*(y[x][0]+47.13)));
        bml=0.08*exp(-y[x][0]/11.0);
        mlss=aml/(aml + bml);
        tauml = 1 / (aml + bml);

        hlss=1/(1+exp((y[x][0]+91)/6.1));
        tauhl=scalehL*233;// modified in HF

        ydot[x][40] = (mlss - y[x][40]) / tauml;//<------------------------??
        ydot[x][41] = (hlss - y[x][41]) / tauhl;

        I_NaL_junc = scaleNaL*Fjunc*GNaL*pow(y[x][40],3)*y[x][41]*(y[x][0]-ena_junc);//modified in HF
        I_NaL_sl =scaleNaL*Fsl*GNaL*pow(y[x][40],3)*y[x][41]*(y[x][0]-ena_sl); //modified in HF
        I_NaL = I_NaL_junc+I_NaL_sl;

        //I_nabk: Na Background Current
        I_nabk_junc = scaleNab*Fjunc*GNaB*(y[x][0]-ena_junc); //modified in HF
        I_nabk_sl = scaleNab*Fsl*GNaB*(y[x][0]-ena_sl);//modified in HF
        I_nabk = I_nabk_junc+I_nabk_sl;

        //I_nak: Na/K Pump Current
        sigma = (exp(Nao/67.3)-1)/7;
        fnak = 1/(1+0.1245*exp(-0.1*y[x][0]*FoRT)+0.0365*sigma*exp(-y[x][0]*FoRT));
        I_nak_junc = scaleNaK*1*Fjunc*IbarNaK*fnak*Ko /(1+pow((KmNaip/y[x][32]),4)) /(Ko+KmKo);//modified in HF
        I_nak_sl = scaleNaK*1*Fsl*IbarNaK*fnak*Ko /(1+pow((KmNaip/y[x][33]),4)) /(Ko+KmKo);//modified in HF
        I_nak = I_nak_junc+I_nak_sl;

        //I_kr: Rapidly Activating K Current
        gkr =1.0*0.035*sqrt(Ko/5.4);
        xrss = 1/(1+exp(-(y[x][0]+10)/5));
        tauxr = 550/(1+exp((-22-y[x][0])/9))*6/(1+exp((y[x][0]-(-11))/9))+230/(1+exp((y[x][0]-(-40))/20));
        ydot[x][12] = (xrss-y[x][12])/tauxr;
        rkr = 1/(1+exp((y[x][0]+74)/24));

        //Omitted commented lines
        I_kr =gkr*y[x][12]*rkr*(y[x][0]-ek);
        //Omitted commented lines
        //Using none markov IKS
        gks_junc=1*0.0035;
        gks_sl=1*0.0035; //FRA
        xsss = 1 / (1+exp(-(y[x][0] + 3.8)/14.25)); // fitting Fra
        tauxs=990.1/(1+exp(-(y[x][0]+2.436)/14.12));
        ydot[x][13] = (xsss-y[x][13])/tauxs;
        I_ks_junc = Fjunc*gks_junc*y[x][13]*y[x][13]*(y[x][0]-eks);
        I_ks_sl = Fsl*gks_sl*y[x][13]*y[x][13]*(y[x][0]-eks);
        I_ks = I_ks_junc+I_ks_sl;

        //I_kp: Plateau K current
        kp_kp = 1/(1+exp(7.488-y[x][0]/5.98));
        I_kp_junc = Fjunc*gkp*kp_kp*(y[x][0]-ek);
        I_kp_sl = Fsl*gkp*kp_kp*(y[x][0]-ek);
        I_kp = I_kp_junc+I_kp_sl;

        //I_to: Transient Outward K Current (slow and fast components)
        // modified for human myocytes
        if (epi==1){
          GtoSlow=1.0*0.13*0.12; //epi
          GtoFast=1.0*0.13*0.88; //epi
        }
        else{
          GtoSlow=0.13*0.3*0.964; //endo
          GtoFast=0.13*0.3*0.036; //endo
        }

        xtoss = 1/(1+exp(-(y[x][0]-19.0)/13));
        ytoss = 1/(1+exp((y[x][0]+19.5)/5));
        tauxtos = 9/(1+exp((y[x][0]+3.0)/15))+0.5;
        tauytos = 800/(1+exp((y[x][0]+60.0)/10))+30;
        ydot[x][8] = (xtoss-y[x][8])/tauxtos;
        ydot[x][9] = (ytoss-y[x][9])/tauytos;

        I_tos = scaleto*GtoSlow*y[x][8]*y[x][9]*(y[x][0]-ek);    // [uA/uF]   modified in HF

        tauxtof = 8.5*exp(-pow(((y[x][0]+45)/50),2))+0.5;
        tauytof = 85*exp((-pow((y[x][0]+40),2)/220))+7;
        ydot[x][10] = (xtoss-y[x][10])/tauxtof;
        ydot[x][11] = (ytoss-y[x][11])/tauytof;

        I_tof = scaleto*GtoFast*y[x][10]*y[x][11]*(y[x][0]-ek);
        I_to = I_tos + I_tof;

        //I_ki: Time-Independent K Current
        aki = 1.02/(1+exp(0.2385*(y[x][0]-ek-59.215)));
        bki =(0.49124*exp(0.08032*(y[x][0]+5.476-ek)) + exp(0.06175*(y[x][0]-ek-594.31))) /(1 + exp(-0.5143*(y[x][0]-ek+4.753)));
        kiss = aki/(aki+bki);
        I_ki =scaleK1*1* 0.35*sqrt(Ko/5.4)*kiss*(y[x][0]-ek); // modified in HF

        //I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
        I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/y[x][36])*(y[x][0]-ecl);
        I_ClCa_sl = Fsl*GClCa/(1+KdClCa/y[x][37])*(y[x][0]-ecl);
        I_ClCa = I_ClCa_junc+I_ClCa_sl;
        I_Clbk = GClB*(y[x][0]-ecl);

        //Omitted commented lines
        fss = 1/(1+exp((y[x][0]+35)/9))+0.6/(1+exp((50-y[x][0])/20));
        //Omitted commented lines

        dss = 1/(1+exp(-(y[x][0]+5)/6.0));
        taud = dss*(1-exp(-(y[x][0]+5)/6.0))/(0.035*(y[x][0]+5));
        tauf = 1/(0.0197*exp( -pow((0.0337*(y[x][0]+14.5)),2) )+0.02);
        ydot[x][4] = (dss-y[x][4])/taud;
        ydot[x][5] = (fss-y[x][5])/tauf;
        ydot[x][6] = 1.7*y[x][36]*(1-y[x][6])-11.9e-3*y[x][6]; // fCa_junc   koff!!!!!!!!
        ydot[x][7] = 1.7*y[x][37]*(1-y[x][7])-11.9e-3*y[x][7]; // fCa_sl
        fcaCaMSL= 0.1/(1+(0.01/y[x][37]));
        fcaCaj= 0.1/(1+(0.01/y[x][36]));
        fcaCaMSL=0;
        fcaCaj= 0;

        ibarca_j = pCa*4*(y[x][0]*Frdy*FoRT) * (0.341*y[x][36]*exp(2*y[x][0]*FoRT)-0.341*Cao) /(exp(2*y[x][0]*FoRT)-1);
        ibarca_sl = pCa*4*(y[x][0]*Frdy*FoRT) * (0.341*y[x][37]*exp(2*y[x][0]*FoRT)-0.341*Cao) /(exp(2*y[x][0]*FoRT)-1);
        ibark = pK*(y[x][0]*Frdy*FoRT)*(0.75*y[x][35]*exp(y[x][0]*FoRT)-0.75*Ko) /(exp(y[x][0]*FoRT)-1);
        ibarna_j = pNa*(y[x][0]*Frdy*FoRT) *(0.75*y[x][32]*exp(y[x][0]*FoRT)-0.75*Nao)  /(exp(y[x][0]*FoRT)-1);
        ibarna_sl = pNa*(y[x][0]*Frdy*FoRT) *(0.75*y[x][33]*exp(y[x][0]*FoRT)-0.75*Nao)  /(exp(y[x][0]*FoRT)-1);
        I_Ca_junc = (Fjunc_CaL*ibarca_j*y[x][4]*y[x][5]*((1-y[x][6])+fcaCaj)*pow(Q10CaL,Qpow))*0.45*1;
        I_Ca_sl = (Fsl_CaL*ibarca_sl*y[x][4]*y[x][5]*((1-y[x][7])+fcaCaMSL)*pow(Q10CaL,Qpow))*0.45*1;
        I_Ca = I_Ca_junc+I_Ca_sl;
        I_CaK = (ibark*y[x][4]*y[x][5]*(Fjunc_CaL*(fcaCaj+(1-y[x][6]))+Fsl_CaL*(fcaCaMSL+(1-y[x][7])))*pow(Q10CaL, Qpow))*0.45*1;
        I_CaNa_junc = (Fjunc_CaL*ibarna_j*y[x][4]*y[x][5]*((1-y[x][6])+fcaCaj)*pow(Q10CaL, Qpow))*0.45*1;
        I_CaNa_sl = (Fsl_CaL*ibarna_sl*y[x][4]*y[x][5]*((1-y[x][7])+fcaCaMSL)*pow(Q10CaL, Qpow))*0.45*1;
        I_CaNa = I_CaNa_junc+I_CaNa_sl;
        I_Catot = I_Ca+I_CaK+I_CaNa;

        //I_ncx: Na/Ca Exchanger flux
        Ka_junc = 1/(1+pow((Kdact/y[x][36]),2));
        Ka_sl = 1/(1+pow((Kdact/y[x][37]),2));
        s1_junc = exp(nu*y[x][0]*FoRT)*y[x][32]*y[x][32]*y[x][32]*Cao;
        s1_sl = exp(nu*y[x][0]*FoRT)*y[x][33]*y[x][33]*y[x][33]*Cao;
        s2_junc = exp((nu-1)*y[x][0]*FoRT)*Nao*Nao*Nao*y[x][36];
        s3_junc = KmCai*Nao*Nao*Nao*(1+pow((y[x][32]/KmNai),3)) + KmNao*KmNao*KmNao*y[x][36]*(1+y[x][36]/KmCai)+KmCao*y[x][32]*y[x][32]*y[x][32]+y[x][32]*y[x][32]*y[x][32]*Cao+Nao*Nao*Nao*y[x][36];
        s2_sl = exp((nu-1)*y[x][0]*FoRT)*Nao*Nao*Nao*y[x][37];
        s3_sl = KmCai*Nao*Nao*Nao*(1+pow((y[x][33]/KmNai), 3)) + KmNao*KmNao*KmNao*y[x][37]*(1+y[x][37]/KmCai)+KmCao*y[x][33]*y[x][33]*y[x][33]+y[x][33]*y[x][33]*y[x][33]*Cao+Nao*Nao*Nao*y[x][37];

        //Omitted commented lines
        I_ncx_junc = scaleNCX*Fjunc*IbarNCX*pow(Q10NCX,Qpow)*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+ksat*exp((nu-1)*y[x][0]*FoRT));//modified in HF
        I_ncx_sl = scaleNCX*Fsl*IbarNCX*pow(Q10NCX,Qpow)*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+ksat*exp((nu-1)*y[x][0]*FoRT)); //modified in HF
        I_ncx = I_ncx_junc+I_ncx_sl;

        //I_pca: Sarcolemmal Ca Pump Current
        I_pca_junc = Fjunc*pow(Q10SLCaP,Qpow)*IbarSLCaP*pow(y[x][36],1.6)/(pow(KmPCa,1.6)+pow(y[x][36],1.6));
        I_pca_sl = Fsl*pow(Q10SLCaP,Qpow)*IbarSLCaP*pow(y[x][37],1.6)/(pow(KmPCa,1.6)+pow(y[x][37],1.6));
        I_pca = I_pca_junc+I_pca_sl;

        //I_cabk: Ca Background Current
        I_cabk_junc =scaleCab*Fjunc*GCaB*(y[x][0]-eca_junc);//modified in HF
        I_cabk_sl = scaleCab*Fsl*GCaB*(y[x][0]-eca_sl);//modified in HF
        I_cabk = I_cabk_junc+I_cabk_sl;

        //SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
        MaxSR = 15;
        MinSR = 1;
        kCaSR = MaxSR - (MaxSR-MinSR)/pow(1+(ec50SR/y[x][31]),2.5);
        koSRCa = koCa/kCaSR;
        kiSRCa = kiCa*kCaSR;
        RI = 1-y[x][14]-y[x][15]-y[x][16];
        ydot[x][14] = (kim*RI-kiSRCa*y[x][36]*y[x][14])-(koSRCa*y[x][36]*y[x][36]*y[x][14]-kom*y[x][15]);   // R
        ydot[x][15] = (koSRCa*y[x][36]*y[x][36]*y[x][14]-kom*y[x][15])-(kiSRCa*y[x][36]*y[x][15]-kim*y[x][16]);// O
        ydot[x][16] = (kiSRCa*y[x][36]*y[x][15]-kim*y[x][16])-(kom*y[x][16]-koSRCa*y[x][36]*y[x][36]*RI);   // I
        J_SRCarel = ks*y[x][15]*(y[x][31]-y[x][36]);          // [mM/ms]
        J_serca = 1*pow(Q10SRCaP,Qpow)*Vmax_SRCaP*(pow((y[x][38]/Kmf), hillSRCaP)-pow((y[x][31]/Kmr),hillSRCaP))/(1+pow((y[x][38]/Kmf),hillSRCaP)+pow((y[x][31]/Kmr),hillSRCaP));
        J_SRleak = 5.348e-6*(y[x][31]-y[x][36]);           //   [mM/ms]

        //Sodium and Calcium Buffering
        ydot[x][17] = kon_na*y[x][32]*(Bmax_Naj-y[x][17])-koff_na*y[x][17];        // NaBj      [mM/ms]
        ydot[x][18] = kon_na*y[x][33]*(Bmax_Nasl-y[x][18])-koff_na*y[x][18];       // NaBsl     [mM/ms]

        //Cytosolic Ca Buffers
        ydot[x][19] = kon_tncl*y[x][38]*(Bmax_TnClow-y[x][19])-koff_tncl*y[x][19];            // TnCL      [mM/ms]
        ydot[x][20] = kon_tnchca*y[x][38]*(Bmax_TnChigh-y[x][20]-y[x][21])-koff_tnchca*y[x][20]; // TnCHc     [mM/ms]
        ydot[x][21] = kon_tnchmg*Mgi*(Bmax_TnChigh-y[x][20]-y[x][21])-koff_tnchmg*y[x][21];   // TnCHm     [mM/ms]
        ydot[x][22] = kon_cam*y[x][38]*(Bmax_CaM-y[x][22])-koff_cam*y[x][22];                 // CaM       [mM/ms]
        ydot[x][23] = kon_myoca*y[x][38]*(Bmax_myosin-y[x][23]-y[x][24])-koff_myoca*y[x][23];    // Myosin_ca [mM/ms]
        ydot[x][24] = kon_myomg*Mgi*(Bmax_myosin-y[x][23]-y[x][24])-koff_myomg*y[x][24];      // Myosin_mg [mM/ms]
        ydot[x][25] = kon_sr*y[x][38]*(Bmax_SR-y[x][25])-koff_sr*y[x][25];                    // SRB       [mM/ms]
        J_CaB_cytosol = ydot[x][19]+ydot[x][20]+ydot[x][21]+ydot[x][22]+ydot[x][23]+ydot[x][24]+ydot[x][25];

        //Junctional and SL Ca Buffers
        ydot[x][26] = kon_sll*y[x][36]*(Bmax_SLlowj-y[x][26])-koff_sll*y[x][26];       // SLLj      [mM/ms]
        ydot[x][27] = kon_sll*y[x][37]*(Bmax_SLlowsl-y[x][27])-koff_sll*y[x][27];      // SLLsl     [mM/ms]
        ydot[x][28] = kon_slh*y[x][36]*(Bmax_SLhighj-y[x][28])-koff_slh*y[x][28];      // SLHj      [mM/ms]
        ydot[x][29] = kon_slh*y[x][37]*(Bmax_SLhighsl-y[x][29])-koff_slh*y[x][29];     // SLHsl     [mM/ms]
        J_CaB_junction = ydot[x][26]+ydot[x][28];
        J_CaB_sl = ydot[x][27]+ydot[x][29];

        // Ion concentrations
        // SR Ca Concentrations
        ydot[x][30] = kon_csqn*y[x][31]*(Bmax_Csqn-y[x][30])-koff_csqn*y[x][30];       // Csqn      [mM/ms]
        ydot[x][31] = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel)-ydot[x][30];         // Ca_sr     [mM/ms] //Ratio 3 leak current

        I_Na_tot_junc = I_Na_junc+I_NaL_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   // [uA/uF] With INaL---
        I_Na_tot_sl = I_Na_sl+I_NaL_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl; // [uA/uF]
        I_Na_tot_sl2 = 3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   // [uA/uF]
        I_Na_tot_junc2 = 3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   // [uA/uF]

        ydot[x][32] = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(y[x][33]-y[x][32])-ydot[x][17];
        ydot[x][33] = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(y[x][32]-y[x][33])+J_na_slmyo/Vsl*(y[x][34]-y[x][33])-ydot[x][18];
        ydot[x][34] = J_na_slmyo/Vmyo*(y[x][33]-y[x][34]);             // [mM/msec]

        //Potassium Concentration
        I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp;     // [uA/uF]
        ydot[x][35] =0; // -I_K_tot*Cmem/(Vmyo*Frdy);

        //Calcium Concentrations
        I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   // [uA/uF]
        I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            // [uA/uF]
        ydot[x][36] = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(y[x][37]-y[x][36]) -J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  // Ca_j
        ydot[x][37] = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(y[x][36]-y[x][37]) + J_ca_slmyo/Vsl*(y[x][38]-y[x][37])-J_CaB_sl;   // Ca_sl
        ydot[x][38] = -J_serca*Vsr/Vmyo-J_CaB_cytosol +J_ca_slmyo/Vmyo*(y[x][37]-y[x][38]);

        /* Current injection */
        if ( (( count <= stim_dur_int) && (stim_myo[x] == 1))&& (num > stim_equil) ){
          I_stim = stim_mag;

        }
        else {

          I_stim = 0.0;
        }

        // Membrane Potential
        I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          // [uA/uF]
        I_Cl_tot = I_ClCa+I_Clbk;                        // [uA/uF]
        I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
        I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
        ydot[x][0] = -(I_tot-I_stim);

        //Update membrane potential and gates
        for (int cc=0; cc<42; cc++) y[x][cc] = y[x][cc]+DT*ydot[x][cc];

        //Calculate APD//////////////////////////////////////////////
        if ( apd90_start[x] == -1 && y[x][0] > upstroke90[x] ){	
          apd90_t1[x] = time;
          apd90_start[x] = 0;
        }
        if ( apd90_start[x] == 0 && y[x][0] < downstroke90[x] ){
          apd90_t2[x] = time;
          apd90_start[x] = 1;
          apd90[x] = apd90_t2[x] - apd90_t1[x];
        }
        //Calculate APD//////////////////////////////////////////////
        if(count%file_filter==0) 
        {
          fprintf(output,"%.12f\t",y[x][0]);
          //  fprintf(output3,"%.12f\t",ydot[x][0]);
        }
      }//End of myo CELLS
      if(count%file_filter==0){
        fprintf(output,"\n");
        //  fprintf(output3,"\n");
      }
      //Start of the fib dimensions
      ///////////////////////////////////////////////////////////////
      for (int x = 0; x<CELLS_fib; x++){
        alpha_K1_fib = 0.1  / (1 + exp(0.06  * (yf[x][0] - EK_fib - 200.)));
        beta_K1_fib = (3.  * exp(0.0002 * (yf[x][0] - EK_fib + 100. )) + 1  * exp(0.1 * (yf[x][0] - EK_fib - 10. ))) / (1  + exp(-(0.5 ) * (yf[x][0] - EK_fib)));
        IK1_fib= GK1_fib * alpha_K1_fib / (alpha_K1_fib + beta_K1_fib) * (yf[x][0] - EK_fib);

        IKv_fib = GKv * yf[x][1] * yf[x][2] * (yf[x][0] - EK_fib);

        INaK_fib = Max_INaK * Ko / (KmK_fib + Ko) * pow(Nai_fib, 1.5) / (pow(KmNa_fib, 1.5 ) + pow(Nai_fib, 1.5 )) * (yf[x][0] + 150. ) / (yf[x][0] + 200. );

        IbNa_fib = GbNa_fib * (yf[x][0] - ENa_fib);

        // Current injection
        if ( (( count <= stim_dur_int)&& (stim_fib[x] == 1))&& (num> stim_equil)) {
          I_stim_fib = -stim_mag;
        }
        else {
          I_stim_fib = 0.0;

        }


        Itotal_fib = IKv_fib + IK1_fib + INaK_fib + IbNa_fib + I_stim_fib;


        //update gates
        r_inf_fib = 1. / (1. + exp(-((yf[x][0]  + 20. )) / 11. ));
        tau_r_fib = 20.3  + 138.  * exp(-(pow((yf[x][0]  + 20.) / 25.9, 2 )));
        yfdot[x][1] = (r_inf_fib - yf[x][1])/tau_r_fib;
        //r_gate_fib = r_inf_fib - (r_inf_fib - r) * exp(-DT / tau_r_fib);

        s_inf_fib = 1./(1. + exp ((yf[x][0]  + 23.0) / 7.0));
        tau_s_fib = 1574. + 5268. * exp(-(pow((yf[x][0]  + 23. ) / 22.7, 2 )));
        //s_gate_fib = s_inf_fib - (s_inf_fib - s) * exp(-DT / tau_s_fib);
        yfdot[x][2] = (s_inf_fib - yf[x][2])/tau_s_fib;


        yfdot[x][0] = (-Itotal_fib);

        //Update membrane potential and gates
        for (int cc=0; cc<3; cc++) yf[x][cc] = yf[x][cc] + DT * yfdot[x][cc];

        if(count%file_filter==0) fprintf(output2,"%.12f\t",yf[x][0]);

      } //end of cell loop fibroblast
      // for (int bb=0; bb<18; bb++) printf("yf[node][%d] = %.12f;\n",bb,yf[0][bb]);
      ////////////////////////////////////////////////////////////////////////////////////////////
      ////END OF CELL LOOP
      // for (int cc=35; cc<=41; cc++) fprintf(output,"%.12f\t",y[x][cc]);
      if(count%file_filter==0) fprintf(output2,"\n");
      //cout <<"Done with cell loops"<<endl;

      PDE(y, yf, neighbors_myo, DynamicGapON);
      //cout << "I out function: " << ydot[0][0] << endl;
      //cout<<"first pde"<<endl;

      PDE(y, yf, neighbors_fib, DynamicGapON);
      //cout<<"second pde"<<endl;

      //Calculate Conduction Velocity/////////////////////////////////////////////////////////////

      if ((CVL_1 > 0) && (CVL_2 > 0)){
        if (y[abs(CVL_1)][0] >= -40 &&  V_myo[abs(CVL_1)]<-40) CVstart = time;
        if (y[abs(CVL_2)][0] >= -40 &&  V_myo[abs(CVL_2)]<-40){
          cout <<"Longitudinal CV = " <<((CVL_cell2 - CVL_cell1)* 25e-4)/(time - CVstart)*1000 << endl;

        }
      }
      if ((CVL_1 > 0) && (CVL_2 < 0)){
        if (y[abs(CVL_1)][0] >= -40 &&  V_myo[abs(CVL_1)]<-40) CVstart = time;
        if (yf[abs(CVL_2)][0] >= -40 &&  V_fib[abs(CVL_2)]<-40){
          cout <<"Longitudinal CV = " <<((CVL_cell2 - CVL_cell1)* 25e-4)/(time - CVstart)*1000 << endl;

        }
      }
      if ((CVL_1 < 0) && (CVL_2 > 0)){
        if (yf[abs(CVL_1)][0] >= -40 &&  V_fib[abs(CVL_1)]<-40) CVstart = time;
        if (y[abs(CVL_2)][0] >= -40 &&  V_myo[abs(CVL_2)]<-40){
          cout <<"Longitudinal CV = " <<((CVL_cell2 - CVL_cell1)* 25e-4)/(time - CVstart)*1000 << endl;
        }
      }
      if ((CVL_1 < 0) && (CVL_2 < 0)){
        if (yf[abs(CVL_1)][0] >= -40 &&  V_fib[abs(CVL_1)]<-40) CVstart = time;
        if (yf[abs(CVL_2)][0] >= -40 &&  V_fib[abs(CVL_2)]<-40){
          cout <<"Longitudinal CV = " <<((CVL_cell2 - CVL_cell1)* 25e-4)/(time - CVstart)*1000 << endl;
        }
      }

      /*


         if (y[CVT_1][0] >= -40 &&  V_myo[CVT_1]<-40) CVstart = time;
         if (y[CVT_2][0] >= -40 &&  V_myo[CVT_2]<-40){
         cout <<"Transverse CV = " <<((CVT_cell2 - CVT_cell1)* 25e-4)/(time - CVstart)*1000 << endl;
         }
         */


      //Stores previous voltage values for myocyte and fibroblast
      for (int aa=0; aa<CELLS_myo; aa++) {V_myo[aa] = y[aa][0];}
      for (int dd=0; dd<CELLS_fib; dd++) {V_fib[dd] = yf[dd][0];}

      // cout<< "Calculate conduction velocity" << endl;

      //Calculate Conduction Velocity/////////////////////////////////////////////////////////////


    } //end of bcl loop


  }// end of stimulation loop
  for (int bb=0; bb<42; bb++) printf("y[node][%d] = %.12f;\n",bb,y[0][bb]);

  //for (int gg=0; gg<3; gg++) printf("yf[node][%d] = %.12f;\n",gg,yf[0][gg]);

  fclose(output);
  fclose(output2);
  fclose(output3);
  // fclose(output4);
  t2=clock();
  cout <<"Runtime: "<< ((float)t2-(float)t1)/(CLOCKS_PER_SEC*60) << " min(s)" <<endl;
  return 0;
}
void PDE( vector< vector<double> >  &Vmyo, vector< vector<double> > &Vfib, vector< vector<int> >  neighbor, bool dynamicGapOn ){//& pass by reference changes the original parameter
  double G, Vj;
  vector< vector<double> > temp;
  vector< vector<double> > temp2;
  temp = Vmyo;
  temp2 = Vfib;
  for (int p = 0; p < neighbor.size(); p++){
    double sum_myo=0,sum_fib=0;
    for (int N = 1; N < neighbor[p].size(); N++){//starts at 1 because first entry (0) is self
      if(neighbor[p][0] > 0){//If self is a myocyte, determine neighbor
        if(neighbor[p][N] > 0){ //If neighbor is a myocyte
          if((neighbor[p][N] == p+2)||(neighbor[p][N] == p-0)){//lazy way of indicating end to end connections
            G = .02;//.06 ; //If the neighbor is end to end then conductance should be 600
          }
          else{
            G = .02; //200 nS converted to uS for the currents
          }
          Vj = (Vmyo[ abs(neighbor[p][N])-1][0] - Vmyo[p][0]);//Neighbor minus 1 to get the right vector index?
          if(dynamicGapOn){G = G * dyngap1(Vj, CX43CX43);}
        }
        else{ //If neighbor is a fibroblast
          //G = 8e-3;
          Vj = (Vfib[ abs(neighbor[p][N])-1][0] - Vmyo[p][0]);//Neighbor minus 1 to get the right vector index?
          G = 0.0008;
          if(dynamicGapOn){ G = G * dyngap1(Vj, CX43CX45);}//This is a Cx43Cx45 channel
        }
        sum_myo += G * Vj;//current in units of uA.
      }
      else {//If self is a fibroblast
        if(neighbor[p][N] > 0){ //If neighbor is a myocyte
          Vj = (Vmyo[ abs(neighbor[p][N])-1][0] - Vfib[p][0]);//Neighbor minus 1 to get the right vector index?
          G = 0.0008; 
          if(dynamicGapOn){G = G * dyngap1(-Vj, CX43CX45);}//This is a Cx45Cx43 channel so just switch sign of Vj in the Cx43Cx45 model
        }
        else{//If neighbor is fibroblast
          Vj = (Vfib[ abs(neighbor[p][N])-1][0] - Vfib[p][0]);//Neighbor minus 1 to get the right vector index?
          G = 0.0008;
          if(dynamicGapOn){G = G * dyngap1(Vj, CX45CX45);}//This is a Cx45Cx45 channel
        }
        sum_fib += G * Vj;//current in units of pA
      }
    }
    if(neighbor[p][0] > 0){
      temp[p][0] +=  sum_myo * (1./125);//125 pF myocyte converted to uF //Still unsure if I should be multiplying by DT
    }
    else{
      temp2[p][0] += sum_fib * (1./25);//25 pF fibroblast 
    }
  }
  Vmyo = temp;
  Vfib = temp2;
}
vector<vector<int> > fileToVector(const char *name){
  vector<vector<int> > result;
  ifstream input (name);
  string lineData;
  while(getline(input, lineData))
  {
    if(lineData.empty()){
      cout<<"File empty"<<endl;
    }
    else{
      double d;
      vector<int> row;
      stringstream lineStream(lineData);

      while (lineStream >> d){
        row.push_back(d);
      }
      result.push_back(row);
    }
  }
  return result;
}
double dyngap(double VJ, cx_model cx_choice){
  double A1, A2;
  double V1, V2;
  double g_res, g_max;
  // double VJ;
  //cout << VJ << endl;
  //  VJ = Vright - Vleft;
  switch (cx_choice){//Opposing gates model of gap junctions by Chen-Izu
    case 0:
      //Parameters for homotypic Cx43 gap junction model
      A1=0.058;
      V1=61.3;
      A2=0.058;
      V2=61.3;
      g_res=0.29;
      g_max=1.00;
      break;
    case 1://Cell on the left is Cx43, cell on the right is Cx45
      A1=0.088;
      V1 = 22.3;
      A2=0.027;
      V2=125.3;
      g_res=0.04;
      g_max=1.03;
      break;
    case 2://Parameters for homotypic Cx45 gap junciton model
      A1=0.110;
      V1 = 10.2;
      A2=0.110;
      V2=10.2;
      g_res=0.06;
      g_max=1.56;
      break;
    default: 
      break;
  }
  return (g_res + ((g_max - g_res)/(1 + exp(A1 * (-VJ - V1)) + exp(A2 * (VJ - V2)))));//Contingent gating model

}
double dyngap1(double VJ, cx_model cx_choice){
  double alpha_1, alpha_2, alpha_3, alpha_4, beta_1, beta_2, beta_3, beta_4;
  double V_LL1, V_LL2, V_LH1, V_LH2, V_HL1, V_HL2, V_HH1, V_HH2, Vj;
  double n_HH = .5;
  double n_LH = .25;
  double n_HL = .25;
  double n_LL = 0;
  double g_1, g_2, delta_g_1, delta_g_2, R_1, R_2, A, B, C, D, determinant, V_1, V_2, gamma_1, gamma_2;
  double g_HH, g_LH, g_HL, g_LL;
  double g_junction, R_junction, I_junction;
  double tolerance = 0.001;
  double t = 0.0;
  int N_channels;

  double V_H, V_L, gamma_H, gamma_L, alpha_coef, beta_coef, V_alpha, V_beta;
  double V_H2, V_L2, gamma_H2, gamma_L2, alpha_coef2, beta_coef2, V_alpha2, V_beta2;

  // double VJ;
  //cout << VJ << endl;
  //  VJ = Vright - Vleft;
  switch (cx_choice){//Opposing gates model of gap junctions by Chen-Izu
    case 0:
      //Gj_inst = 0.0733
      N_channels=13;//13.643 floored
      //Parameters for homotypic Cx43 gap junction model
      V_H= 145.9;
      V_L= 299;
      gamma_H= 146.6;
      gamma_L= 13.1;
      alpha_coef= 181.5e-3;
      beta_coef= 0.007e-3;
      V_alpha= 8.437;
      V_beta= 8.675;

      V_H2= 145.9; //units mV*/
      V_L2= 299;
      gamma_H2= 146.6;
      gamma_L2= 13.1;
      alpha_coef2= 181.5e-3; //units ms
        beta_coef2= 0.007e-3;
        V_alpha2= 8.437;
        V_beta2= 8.675;
        break;
    case 1://Cell on the left is Cx43, cell on the right is Cx45
      //Gj_inst = 0.0471
      N_channels=21;//should be 21.23
      V_H= 55.9;
      V_L= 299;
      gamma_H= 126.6;
      gamma_L= 13.1;
      alpha_coef= 181.5e-3;
      beta_coef= 0.007e-3;
      V_alpha= 8.437;
      V_beta= 8.675;

      V_H2= 615; //units mV*/
      V_L2= 299;
      gamma_H2= 75;
      gamma_L2= 4;
      alpha_coef2= 181.5e-3; //units ms
      beta_coef2= 0.007e-3;
      V_alpha2= 2.437;
      V_beta2= 6.675;
      break;
    case 2://Parameters for homotypic Cx45 gap junciton model
      //Gj_inst = 0.0284
      N_channels=35;//should be 35.211
      V_H=  113.86;
      V_L= 345.23;
      gamma_H= 57;
      gamma_L= 4;
      alpha_coef=  0.2468;
      beta_coef= 0.000469;
      V_alpha= 4.6867;
      V_beta= 75.9036;

      V_H2=  113.86;
      V_L2= 345.23;
      gamma_H2= 57;
      gamma_L2= 4;
      alpha_coef2=  0.2468;
      beta_coef2= 0.000469;
      V_alpha2= 4.6867;
      V_beta2= 75.9036;
      break;
    default: 
      break;
  }
  // Procedure for calculating g_LL

  g_1 = 10.0;													// Initial guess for g_1
  g_2 = 10.0;													// Initial guess for g_2

  V_1 = V_L;													// State of hemichannel 1
  V_2 = V_L2;													// State of hemichannel 2
  gamma_1 = gamma_L;
  gamma_2 = gamma_L2;

  R_1 = gamma_1*exp((Vj/V_1)*(g_2/(g_1+g_2)))-g_1;								// Residual equation 1
  R_2 = gamma_2*exp(-(Vj/V_2)*(g_1/(g_1+g_2)))-g_2;								// Residual equation 2

  while ( fabs(R_1) > tolerance || fabs(R_2) > tolerance )
  {
    A = -((gamma_1*Vj*g_2/(V_1*(g_1+g_2)*(g_1+g_2)))*exp((Vj/V_1)*(g_2/(g_1+g_2))))-1.0;			// Partial derivative of equation 1 wrt g_1
    B = (gamma_1*Vj*g_1/(V_1*(g_1+g_2)*(g_1+g_2)))*exp((Vj/V_1)*(g_2/(g_1+g_2)));				// Partial derivative of equation 1 wrt g_2
    C = -(gamma_2*Vj*g_2/(V_2*(g_1+g_2)*(g_1+g_2)))*exp(-(Vj/V_2)*(g_1/(g_1+g_2)));				// Partial derivative of equation 2 wrt g_1
    D = ((gamma_2*Vj*g_1/(V_2*(g_1+g_2)*(g_1+g_2)))*exp(-(Vj/V_2)*(g_1/(g_1+g_2))))-1.0;			// Partial derivative of equation 2 wrt g_2

    determinant = A*D-B*C;											// Determinant
    delta_g_1 = (1.0/determinant)*(D*(-R_1)-B*(-R_2));							// Change in g_1
    delta_g_2 = (1.0/determinant)*(-C*(-R_1)+A*(-R_2));							// Change in g_2

    g_1 = g_1+delta_g_1;											// Updated value of g_1
    g_2 = g_2+delta_g_2;											// Updated value of g_2

    R_1 = gamma_1*exp((Vj/V_1)*(g_2/(g_1+g_2)))-g_1;							// Residual for equation 1
    R_2 = gamma_2*exp(-(Vj/V_2)*(g_1/(g_1+g_2)))-g_2;							// Residual for equation 2
  }

  g_LL = (g_1*g_2)/(g_1+g_2);											// Junction conductance in the LL state
  V_LL1 = -Vj*(g_2/(g_1+g_2));											// Voltage drop across the first hemichannel
  V_LL2 = Vj*(g_1/(g_1+g_2));											// Voltage drop across the second hemichannel


  // Procedure for calculating g_LH

  g_1 = 10.0;
  g_2 = 10.0;

  V_1 = V_L;
  V_2 = V_H2;
  gamma_1 = gamma_L;
  gamma_2 = gamma_H2;

  R_1 = gamma_1*exp((Vj/V_1)*(g_2/(g_1+g_2)))-g_1;
  R_2 = gamma_2*exp(-(Vj/V_2)*(g_1/(g_1+g_2)))-g_2;

  while ( fabs(R_1) > tolerance || fabs(R_2) > tolerance )
  {
    A = -((gamma_1*Vj*g_2/(V_1*(g_1+g_2)*(g_1+g_2)))*exp((Vj/V_1)*(g_2/(g_1+g_2))))-1.0;
    B = (gamma_1*Vj*g_1/(V_1*(g_1+g_2)*(g_1+g_2)))*exp((Vj/V_1)*(g_2/(g_1+g_2)));
    C = -(gamma_2*Vj*g_2/(V_2*(g_1+g_2)*(g_1+g_2)))*exp(-(Vj/V_2)*(g_1/(g_1+g_2)));
    D = ((gamma_2*Vj*g_1/(V_2*(g_1+g_2)*(g_1+g_2)))*exp(-(Vj/V_2)*(g_1/(g_1+g_2))))-1.0;

    determinant = A*D-B*C;
    delta_g_1 = (1.0/determinant)*(D*(-R_1)-B*(-R_2));
    delta_g_2 = (1.0/determinant)*(-C*(-R_1)+A*(-R_2));

    g_1 = g_1+delta_g_1;
    g_2 = g_2+delta_g_2;

    R_1 = gamma_1*exp((Vj/V_1)*(g_2/(g_1+g_2)))-g_1;
    R_2 = gamma_2*exp(-(Vj/V_2)*(g_1/(g_1+g_2)))-g_2;
  }

  g_LH = (g_1*g_2)/(g_1+g_2);
  V_LH1 = -Vj*(g_2/(g_1+g_2));
  V_LH2 = Vj*(g_1/(g_1+g_2));


  // Procedure for calculating g_HL

  g_1 = 10.0;
  g_2 = 50.0;

  V_1 = V_H;
  V_2 = V_L2;
  gamma_1 = gamma_H;
  gamma_2 = gamma_L2;

  R_1 = gamma_1*exp((Vj/V_1)*(g_2/(g_1+g_2)))-g_1;
  R_2 = gamma_2*exp(-(Vj/V_2)*(g_1/(g_1+g_2)))-g_2;

  while ( fabs(R_1) > tolerance || fabs(R_2) > tolerance )
  {
    A = -((gamma_1*Vj*g_2/(V_1*(g_1+g_2)*(g_1+g_2)))*exp((Vj/V_1)*(g_2/(g_1+g_2))))-1.0;
    B = (gamma_1*Vj*g_1/(V_1*(g_1+g_2)*(g_1+g_2)))*exp((Vj/V_1)*(g_2/(g_1+g_2)));
    C = -(gamma_2*Vj*g_2/(V_2*(g_1+g_2)*(g_1+g_2)))*exp(-(Vj/V_2)*(g_1/(g_1+g_2)));
    D = ((gamma_2*Vj*g_1/(V_2*(g_1+g_2)*(g_1+g_2)))*exp(-(Vj/V_2)*(g_1/(g_1+g_2))))-1.0;

    determinant = A*D-B*C;
    delta_g_1 = (1.0/determinant)*(D*(-R_1)-B*(-R_2));
    delta_g_2 = (1.0/determinant)*(-C*(-R_1)+A*(-R_2));

    g_1 = g_1+delta_g_1;
    g_2 = g_2+delta_g_2;

    R_1 = gamma_1*exp((Vj/V_1)*(g_2/(g_1+g_2)))-g_1;
    R_2 = gamma_2*exp(-(Vj/V_2)*(g_1/(g_1+g_2)))-g_2;
  }

  g_HL = (g_1*g_2)/(g_1+g_2);
  V_HL1 = -Vj*(g_2/(g_1+g_2));
  V_HL2 = Vj*(g_1/(g_1+g_2));


  // Procedure for calculating g_HH

  g_1 = 10.0;
  g_2 = 10.0;

  V_1 = V_H;
  V_2 = V_H2;
  gamma_1 = gamma_H;
  gamma_2 = gamma_H2;

  R_1 = gamma_1*exp((Vj/V_1)*(g_2/(g_1+g_2)))-g_1;
  R_2 = gamma_2*exp(-(Vj/V_2)*(g_1/(g_1+g_2)))-g_2;

  while ( fabs(R_1) > tolerance || fabs(R_2) > tolerance )
  {
    A = -((gamma_1*Vj*g_2/(V_1*(g_1+g_2)*(g_1+g_2)))*exp((Vj/V_1)*(g_2/(g_1+g_2))))-1.0;
    B = (gamma_1*Vj*g_1/(V_1*(g_1+g_2)*(g_1+g_2)))*exp((Vj/V_1)*(g_2/(g_1+g_2)));
    C = -(gamma_2*Vj*g_2/(V_2*(g_1+g_2)*(g_1+g_2)))*exp(-(Vj/V_2)*(g_1/(g_1+g_2)));
    D = ((gamma_2*Vj*g_1/(V_2*(g_1+g_2)*(g_1+g_2)))*exp(-(Vj/V_2)*(g_1/(g_1+g_2))))-1.0;

    determinant = A*D-B*C;
    delta_g_1 = (1.0/determinant)*(D*(-R_1)-B*(-R_2));
    delta_g_2 = (1.0/determinant)*(-C*(-R_1)+A*(-R_2));

    g_1 = g_1+delta_g_1;
    g_2 = g_2+delta_g_2;

    R_1 = gamma_1*exp((Vj/V_1)*(g_2/(g_1+g_2)))-g_1;
    R_2 = gamma_2*exp(-(Vj/V_2)*(g_1/(g_1+g_2)))-g_2;
  }

  g_HH = (g_1*g_2)/(g_1+g_2);
  V_HH1 = -Vj*(g_2/(g_1+g_2));
  V_HH2 = Vj*(g_1/(g_1+g_2));


  alpha_1 = 2.0*alpha_coef/(1.0+exp(-V_LH1/V_alpha));
  alpha_2 = 2.0*alpha_coef2/(1.0+exp(-V_HL2/V_alpha2));
  alpha_3 = 2.0*alpha_coef/(1.0+exp(-V_LL1/V_alpha));
  alpha_4 = 2.0*alpha_coef2/(1.0+exp(-V_LL2/V_alpha2));

  beta_1 = beta_coef*exp(-V_HH1/V_beta);
  beta_2 = beta_coef2*exp(-V_HH2/V_beta2);
  beta_3 = beta_coef*exp(-V_HL1/V_beta);
  beta_4 = beta_coef2*exp(-V_LH2/V_beta2);

  n_LL = n_LL+DT*(beta_4*n_LH+beta_3*n_HL-(alpha_3+alpha_4)*n_LL);
  n_LH = n_LH+DT*(beta_1*n_HH-(alpha_1+beta_4)*n_LH+alpha_4*n_LL);
  n_HL = n_HL+DT*(beta_2*n_HH-(alpha_2+beta_3)*n_HL+alpha_3*n_LL);
  n_HH = n_HH+DT*(-(beta_1+beta_2)*n_HH+alpha_1*n_LH+alpha_2*n_HL);

  return N_channels*(n_LL*g_LL+n_LH*g_LH+n_HL*g_HL+n_HH*g_HH)*pow(10.0,-3);	//units nanoSiemens
}
double dyngap2(double VJ, cx_model cx_choice){
  //State variables
  double  n1H2H,n1L2H,n1H2L,n1L2L;
  double n[4];
  //Conductance of each channel state
  double g1H2H,g1L2H,g1H2L,g1L2L;
  //Voltage across each hemichannel
  double V1H,V2H,V1L,V2L;
  double V1H2H,V1L2L,V1H2L,V1L2H;
  //Conductance of high and low hemichannel conductance state
  double G1H,G2H,G1L,G2L;
  //Rate constants (1/ms)
  double alpha_1, beta_1, alpha_2, beta_2;
  double a1, a2, a3, a4, b1, b2, b3, b4;
  //Alpha and beta potential units (mV)
  double V1a, V1b;
  double V2a, V2b;

  // cout << "In dyngap" << endl;
  //cout << VJ <<endl;
  //  VJ = Vright - Vleft;
  switch (cx_choice){
    case 0:
      //Parameters for homotypic Cx43 gap junction model
      //voltage across the channels units: mV
      V1H = 145.9;
      V2H = 145.9;
      V1L = 299.0;
      V2L = 299.0;
      //Conductance at high state and low state (pS)
      G1H = 146.6;
      G2H = 146.6;
      G1L = 13.1;
      G2L = 13.1;
      //Rate constants unit: 1/ms
      alpha_1 = 181.5e-3;
      beta_1 = 0.007e-3;
      alpha_2 = 181.5e-3;
      beta_2 = 0.007e-3;

      //alpha and beta potential units (mV)
      V1a = 8.437;
      V1b = 8.675;
      V2a = 8.437;
      V2b = 8.675;
      break;
    case 1:
      //Cell on the left is Cx43, cell on the right is Cx45
      //Parameters for homotypic Cx43 gap junction model

      //voltage across the channels units: mV
      V1H = 55.9;
      V2H = 615;
      V1L = 299.0;
      V2L = 299.0;
      //Conductance at high state and low state (pS)
      G1H = 126.6;
      G2H = 75;
      G1L = 13.1;
      G2L = 4;
      //Rate constants unit: 1/ms
      alpha_1 = 181.5e-3;
      beta_1 = 0.007e-3;
      alpha_2 = 181.5e-3;
      beta_2 = 0.007e-3;

      //alpha and beta potential units (mV)
      V1a = 8.437;
      V1b = 8.675;
      V2a = 2.437;
      V2b = 6.675;
      break;
    case 2:
      //Parameters for homotypic Cx45 gap junction model
      //voltage across the channels units: mV
      V1H = 113.86;
      V2H = 113.86;
      V1L = 345.23;
      V2L = 345.23;
      //conductance at high state and low state (pS)
      G1H = 57;
      G2H = 57;
      G1L = 4;
      G2L = 4;
      //Rate constants unit: 1/ms
      alpha_1 = 0.2468;
      beta_1 = 0.000469;
      alpha_2 = 0.2468;
      beta_2 = 0.000469;
      //alpha and beta potential units (mV)
      V1a = 4.6867;
      V2a = 4.6867;
      V1b = 75.9036;
      V2b = 75.9036;
      break;
    default: 
      cout << "Default case chosen in ERROR"<< endl;
      break;
  }

  //nXY is the probabilities of the channels to be in the respective states
  n1H2H=1.0;
  n1L2H=0.0;
  n1H2L=0.0;
  n1L2L=0.0;

  //gXY is the conductance of the channel state, it is VJ dependent
  g1H2H = 0.0;
  g1L2H = 0.0;
  g1H2L = 0.0;
  g1L2L = 0.0;

  V1H2H=0.0;
  V1L2H=0.0;
  V1H2L=0.0;
  V1L2L=0.0;

  a1=0;
  a2=0;
  a3=0;
  a4=0;
  b1=0;
  b2=0;
  b3=0;
  b4=0;

  struct voltage_params params = {G1H, G2H, G1L, G2L, V1H, V2H, V1L, V2L, VJ};

  FDF_hh.f = &voltagehh;
  FDF_hh.df = &voltagehh_deriv;
  FDF_hh.fdf = &voltagehh_fdf;
  FDF_hh.params = &params;
  V1H2H = rootSolver(FDF_hh,V1H2H);//1H2H
  double V1H2H1 = V1H2H; 
  double V1H2H2 = VJ + V1H2H1;

  FDF_ll.f = &voltagell;
  FDF_ll.df = &voltagell_deriv;
  FDF_ll.fdf = &voltagell_fdf;
  FDF_ll.params = &params;
  V1L2L = rootSolver(FDF_ll,V1L2L);//1L2L
  double V1L2L1 = V1L2L; 
  double V1L2L2 = VJ + V1L2L1;

  FDF_hl.f = &voltagehl;
  FDF_hl.df = &voltagehl_deriv;
  FDF_hl.fdf = &voltagehl_fdf;
  FDF_hl.params = &params;
  V1H2L = rootSolver(FDF_hl,V1H2L);//1H2L
  double V1H2L1 = V1H2L; 
  double V1H2L2 = VJ + V1H2L1;

  FDF_lh.f = &voltagelh;
  FDF_lh.df = &voltagelh_deriv;
  FDF_lh.fdf = &voltagelh_fdf;
  FDF_lh.params = &params;
  V1L2H = rootSolver(FDF_lh,V1L2H); //1L2H
  double V1L2H1 = V1L2H;
  double V1L2H2 = VJ + V1L2H1;

  //cout << "Calculate roots" << endl;

  // calculate homotypic conductances (pS)
  g1H2H = (G1H * exp(-V1H2H1 / V1H)) * (G2H * exp(-V1H2H2 / V2H)) / ((G1H * exp(-V1H2H1 / V1H)) + (G2H * exp(-V1H2H2 / V2H)));
  g1L2L = (G1L * exp(-V1L2L1 / V1L)) * (G2L * exp(-V1L2L2 / V2L)) / ((G1L * exp(-V1L2L1 / V1L)) + (G2L * exp(-V1L2L2 / V2L)));

  // calculate heterotypic conductances (pS)
  g1L2H = (G1L * exp(-V1L2H1 / V1L)) * (G2H * exp(-V1L2H2 / V2H)) / ((G1L * exp(-V1L2H1 / V1L)) + (G2H * exp(-V1L2H2 / V2H)));
  g1H2L = (G1H * exp(-V1H2L1 / V1H)) * (G2L * exp(-V1H2L2 / V2L)) / ((G1H * exp(-V1H2L1 / V1H)) + (G2L * exp(-V1H2L2 / V2L)));


  a1 = 2 * alpha_1 / (1 + exp(-V1L2H1 / V1a));
  a2 = 2 * alpha_2 / (1 + exp(-V1H2L2 / V2a));
  a3 = 2 * alpha_1 / (1 + exp(-V1L2L1 / V1a));
  a4 = 2 * alpha_2 / (1 + exp(-V1L2L2 / V2a));
  b1 = beta_1 * exp(-V1H2H1 / V1b);
  b2 = beta_2 * exp(-V1H2H2 / V2b);
  b3 = beta_1 * exp(-V1H2L1 / V1b);
  b4 = beta_2 * exp(-V1L2H2 / V2b);

  n[0] = n1H2H;
  n[1] = n1L2H;
  n[2] = n1H2L;
  n[3] = n1L2L;

  n[0] = n[0] + DT* ( -(b1 + b2) * n[0] + a1 * n[1] + a2 * n[2]);
  n[1] = n[1] + DT* (b1 * n[0] - (a1 + b4) * n[1] + a4* n[3]);
  n[2] = n[2] + DT* (b2 * n[0] - (a2 + b3) * n[2] + a3* n[3]);
  n[3] = n[3] + DT* (b4 * n[1] + b3 * n[2] - (a3 + a4)* n[3]);

  n1H2H = n[0];
  n1L2H = n[1];
  n1H2L = n[2];
  n1L2L = n[3];
  //Originals units are pS, multiply by 1e-3 to convert to nS
  double pico_to_nano_siemens = 1e-3;
  //cout <<"Scale mag"<<endl;
  return (g1H2H * n1H2H + g1L2H * n1L2H + g1H2L * n1H2L + g1L2L * n1L2L)* pico_to_nano_siemens;

}
/*
   void calcAPD(double time, double voltage, int n){
   if ( apd90_start[n] == -1 && voltage > upstroke90[n] ){	
   apd90_t1[n] = time;
   apd90_start[n] = 0;
   }
   if ( apd90_start[n] == 0 && voltage < downstroke90[n] ){
   apd90_t2[n] = time;
   apd90_start[n] = 1;
   apd90[n] = apd90_t2[n] - apd90_t1[n];
   }
   }
   */

double rootSolver(gsl_function_fdf FDF, double x){
  int status;
  int iter = 0, max_iter = 100;
  double x0;

  gsl_root_fdfsolver_set (SolverPtr, &FDF, x);//x is the initial guess
  do
  {
    iter++;
    status = gsl_root_fdfsolver_iterate (SolverPtr);//These functions perform a single iteration of the solver s. 
    x0 = x;
    x = gsl_root_fdfsolver_root (SolverPtr);//These functions return the current estimate of the root for the solver s.
    status = gsl_root_test_delta (x, x0, 0, 1e-8);//search stopping algorithm

    /* if (status == GSL_SUCCESS)
       printf ("Converged:\n");
       */
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  return x;
}
double voltagehh (double vhh, void *params)
{
  struct voltage_params *p  = (struct voltage_params *) params;

  double G1H = p->G1H;
  double G2H = p->G2H;
  double V1H = p->V1H;
  double V2H = p->V2H;
  double VJ = p->VJ;

  return vhh + VJ / (G1H / G2H * exp((vhh * (V1H - V2H) + VJ * V1H) / (V1H * V2H)	+ 1));
  //   return vhh	+ (VJ / (G1H / G2H * exp((vhh * (V1H - V2H) + V1H * VJ) / (V2H * V1H)) + 1));

}
double voltagehh_deriv (double vhh, void *params)
{
  struct voltage_params *p  = (struct voltage_params *) params;

  double G1H = p->G1H;
  double G2H = p->G2H;
  double V1H = p->V1H;
  double V2H = p->V2H;
  double VJ = p->VJ;

  return 1 - 1 * VJ * G2H * (V1H - V2H)* exp(-(VJ * V2H + vhh * (V1H - V2H))/(V1H * V2H) - 1)/(G1H * V1H * V2H);
  //return (VJ * G1H * G2H * (V2H - V1H) * exp((VJ * V1H + vhh * (V1H + V2H))/(V1H * V2H)))/(V1H * V2H * pow(G1H * exp((VJ + vhh)/V2H) + G2H * exp(vhh/V1H), 2)) + 1;

}
void voltagehh_fdf (double vhh, void *params, double *_vhh, double *_dvhh)
{
  //since it is better to calculate both the function and its derivative at the same time
  //this function combines the calculation derivative is in parameter _dvhl
  struct voltage_params *p  = (struct voltage_params *) params;

  double G1H = p->G1H;
  double G2H = p->G2H;
  double V1H = p->V1H;
  double V2H = p->V2H;
  double VJ = p->VJ;

  *_vhh = vhh + VJ / (G1H / G2H * exp((vhh * (V1H - V2H) + VJ * V1H) / (V1H * V2H)	+ 1));
  *_dvhh = 1 - 1 * VJ * G2H * (V1H - V2H)* exp(-(VJ * V2H + vhh * (V1H - V2H))/(V1H * V2H) - 1)/(G1H * V1H * V2H);
  // *_vhh  =   vhh	+ (VJ / (G1H / G2H * exp((vhh * (V1H - V2H) + V1H * VJ) / (V2H * V1H)) + 1));
  // *_dvhh =  (VJ * G1H * G2H * (V2H - V1H) * exp((VJ * V1H + vhh * (V1H + V2H))/(V1H * V2H)))/(V1H * V2H * pow(G1H * exp((VJ + vhh)/V2H) + G2H * exp(vhh/V1H), 2)) + 1;

}
double voltagell (double vll, void *params)
{
  struct voltage_params *p  = (struct voltage_params *) params;

  double G1L = p->G1L;
  double G2L = p->G2L;
  double V1L = p->V1L;
  double V2L = p->V2L;
  double VJ = p->VJ;

  return vll + VJ / (G1L / G2L * exp((vll * (V1L - V2L) + VJ * V1L) / (V1L * V2L)	+ 1));
}
double voltagell_deriv (double vll, void *params)
{
  struct voltage_params *p  = (struct voltage_params *) params;

  double G1L = p->G1L;
  double G2L = p->G2L;
  double V1L = p->V1L;
  double V2L = p->V2L;
  double VJ = p->VJ;

  return 1 - 1 * VJ * G2L * (V1L - V2L)* exp(-(VJ * V2L + vll * (V1L - V2L))/(V1L * V2L) - 1)/(G1L * V1L * V2L);
}

void voltagell_fdf (double vll, void *params, double *_vll, double *_dvll)
{

  struct voltage_params *p  = (struct voltage_params *) params;

  double G1L = p->G1L;
  double G2L = p->G2L;
  double V1L = p->V1L;
  double V2L = p->V2L;
  double VJ = p->VJ;

  *_vll =  vll + VJ / (G1L / G2L * exp((vll * (V1L - V2L) + VJ * V1L) / (V1L * V2L)	+ 1));
  *_dvll = 1 - 1 * VJ * G2L * (V1L - V2L)* exp(-(VJ * V2L + vll * (V1L - V2L))/(V1L * V2L) - 1)/(G1L * V1L * V2L);
}

double voltagehl (double vhl, void *params)
{

  struct voltage_params *p  = (struct voltage_params *) params;//This is typecasting

  double G1H = p->G1H;
  double G2L = p->G2L;
  double V1H = p->V1H;
  double V2L = p->V2L;
  double VJ = p->VJ;


  return vhl + VJ / (G1H / G2L * exp((vhl * (V1H - V2L) + VJ * V1H) / (V1H * V2L)	+ 1));
}
double voltagehl_deriv (double vhl, void *params)
{
  struct voltage_params *p  = (struct voltage_params *) params;

  double G1H = p->G1H;
  double G2L = p->G2L;
  double V1H = p->V1H;
  double V2L = p->V2L;
  double VJ = p->VJ;


  return 1 - 1 * VJ * G2L * (V1H - V2L)* exp(-(VJ * V2L + vhl * (V1H - V2L))/(V1H * V2L) - 1)/(G1H * V1H * V2L);
}

void voltagehl_fdf (double vhl, void *params, double *_vhl, double *_dvhl)
{
  struct voltage_params *p  = (struct voltage_params *) params;

  double G1H = p->G1H;
  double G2L = p->G2L;
  double V1H = p->V1H;
  double V2L = p->V2L;
  double VJ = p->VJ;

  *_vhl = vhl + VJ / (G1H / G2L * exp((vhl * (V1H - V2L) + VJ * V1H) / (V1H * V2L)	+ 1));
  *_dvhl =  1 - 1 * VJ * G2L * (V1H - V2L)* exp(-(VJ * V2L + vhl * (V1H - V2L))/(V1H * V2L) - 1)/(G1H * V1H * V2L);
}

double voltagelh (double vlh, void *params)
{
  struct voltage_params *p  = (struct voltage_params *) params;

  double G2H = p->G2H;
  double G1L = p->G1L;
  double V2H = p->V2H;
  double V1L = p->V1L;
  double VJ = p->VJ;

  return vlh	+ (VJ / (G1L / G2H * exp((vlh * (V1L - V2H) + V1L * VJ) / (V2H * V1L)) + 1));
}
double voltagelh_deriv (double vlh, void *params)
{
  struct voltage_params *p  = (struct voltage_params *) params;

  double G2H = p->G2H;
  double G1L = p->G1L;
  double V2H = p->V2H;
  double V1L = p->V1L;
  double VJ = p->VJ;


  return (VJ * G1L * G2H * (V2H - V1L) * exp((VJ * V1L + vlh * (V1L + V2H))/(V1L * V2H)))/(V1L * V2H * pow(G1L * exp((VJ + vlh)/V2H) + G2H * exp(vlh/V1L), 2)) + 1;
}

void voltagelh_fdf (double vlh, void *params, double *_vlh, double *_dvlh)
{
  struct voltage_params *p  = (struct voltage_params *) params;

  double G2H = p->G2H;
  double G1L = p->G1L;
  double V2H = p->V2H;
  double V1L = p->V1L;
  double VJ = p->VJ;


  *_vlh = vlh	+ (VJ / (G1L / G2H * exp((vlh * (V1L - V2H) + V1L * VJ) / (V2H * V1L)) + 1));
  *_dvlh = (VJ * G1L * G2H * (V2H - V1L) * exp((VJ * V1L + vlh * (V1L + V2H))/(V1L * V2H)))/(V1L * V2H * pow(G1L * exp((VJ + vlh)/V2H) + G2H * exp(vlh/V1L), 2)) + 1;
}
