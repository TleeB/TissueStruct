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
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
using namespace std;

//Model characteristics
const double epi = 0;//Cell type of model is originally endo
const double HF = 0; //To simulate the heart failure model use HF=1

//Global stimulation parameters
int num_stim=5;
double stim_mag=39.5;
double stim_dur=3;
double bcl=1000;
double DT = 0.005;
char output_file_name[30] = "testGrandi.dat";
char output_file_name2[30] = "testGrandi2.dat";
const char *neighbors_myo_file = "neighbors_myo_file.dat";
const char *neighbors_fib_file = "neighbors_fib_file.dat";
char stim_myo_file[30] = "stim_myo_file.dat";
char stim_fib_file[30] = "stim_fib_file.dat";


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
void PDE(double *I, double *V,vector <vector<int> >  &);
vector<vector<int> > fileToVector(const char *name);
int main(){

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


  FILE *output, *output2;
  output = fopen(output_file_name, "w");
  output2 = fopen(output_file_name2, "w");

  vector <vector<int> > neighbors_myo = fileToVector(neighbors_myo_file);
  vector <vector<int> > neighbors_fib = fileToVector(neighbors_fib_file);
  vector <vector<int> > stim_myo = fileToVector(stim_myo_file);
  vector <vector<int> > stim_fib = fileToVector(stim_fib_file);

  //Number of myocyte cells and fibroblast cells
  const int CELLS_myo = stim_myo[0].size();
  const int CELLS_fib = stim_fib[0].size();
  cout << "Fibroblasts: " << CELLS_fib << endl;
  cout << "Myocytes: " << CELLS_myo << endl;


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

  double y[CELLS_myo][41], ydot[CELLS_myo][41];
  double yf[CELLS_fib][3], yfdot[CELLS_fib][3];

  //Fibroblast parameters ----------------------

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


  //Initialize myocyte variables
  for (int node = 0; node<CELLS_myo; node++){
    y[node][0]=-8.1386313e+1;//Duplicate of membrane potential
    y[node][1]=3.8467219e-3;
    y[node][2]=6.2384201e-1;
    y[node][3]=6.2163026e-1;
    y[node][4]=2.9571680e-6;
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
    y[node][31]=5.8241866e-1; //5.545201e-1;
    y[node][32]=8.2712566e+0; //8.80329;
    y[node][33]=8.2702314e+0; //8.80733;
    y[node][34]=8.2704125e+0;  //8.80853;
    y[node][35]=120;
    y[node][36]=1.8131347e-4;
    y[node][37]=1.0826539e-4;
    y[node][38]=9.1363224e-5;
    y[node][39]=-8.1386313e+1;
    y[node][40]=2.8110845e-3;
    y[node][41]=1.6274168e-1;
  }

  //Initialize fibroblast variables
  for (int node = 0; node<CELLS_fib; node++){
    yf[node][0] = -49.21; //CellML from Dr. Fink
    yf[node][1] = 0.0; //r gate
    yf[node][2] = 1.0; //s gate
  }

  // t=tmin;

  //Local pacing & stimulus parameters
  int count, bcl_int, stim_dur_int;
  double time;
  int num = 0;
  stim_dur_int   = stim_dur/DT;
  bcl_int = bcl/DT;



  while ( num < num_stim )	{
    num = num + 1;//Number is the num of beats.
    count = 0;//Integer representation of time

    //for(t=tmin + DT; t<=tmax; t+=DT){
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
        if ( (( count <= stim_dur_int) && (stim_myo[0][x] == 1))&& num > 2 ){
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
       // ydot[x][0] = ydot[x][0]; //duplicate membrane potential

        //Update membrane potential and gates
        for (int cc=0; cc<=41; cc++) y[x][cc] = y[x][cc]+DT*ydot[x][cc];

        if(count%100==0) fprintf(output,"%.12f\t",y[x][0]);


      }//End of myo CELLS
      if(count%100==0) fprintf(output,"\n");

      //Start of the fibroblast dimensions
      ///////////////////////////////////////////////////////////////////////////////////////////
      for (int x = 0; x<CELLS_fib; x++){

        alpha_K1_fib = 0.1  / (1 + exp(0.06  * (yf[x][0] - EK_fib - 200.)));
        beta_K1_fib = (3.  * exp(0.0002 * (yf[x][0] - EK_fib + 100. )) + 1  * exp(0.1 * (yf[x][0] - EK_fib - 10. ))) / (1  + exp(-(0.5 ) * (yf[x][0] - EK_fib)));
        IK1_fib= GK1_fib * alpha_K1_fib / (alpha_K1_fib + beta_K1_fib) * (yf[x][0] - EK_fib);

        IKv_fib = GKv * yf[x][1] * yf[x][2] * (yf[x][0] - EK_fib);

        INaK_fib = Max_INaK * Ko / (KmK_fib + Ko) * pow(Nai_fib, 1.5) / (pow(KmNa_fib, 1.5 ) + pow(Nai_fib, 1.5 )) * (yf[x][0] + 150. ) / (yf[x][0] + 200. );

        IbNa_fib = GbNa_fib * (yf[x][0] - ENa_fib);

        Itotal_fib = IKv_fib + IK1_fib + INaK_fib + IbNa_fib + I_stim_fib;//Removed stimulation of fibroblasts


        //update gates
        r_inf_fib = 1. / (1. + exp(-((yf[x][0]  + 20. )) / 11. ));
        tau_r_fib = 20.3  + 138.  * exp(-(pow((yf[x][0]  + 20.) / 25.9, 2 )));
        yfdot[x][1] = (r_inf_fib - yf[x][1])/tau_r_fib;
        //r_gate_fib = r_inf_fib - (r_inf_fib - r) * exp(-DT / tau_r_fib);

        s_inf_fib = 1./(1. + exp ((yf[x][0]  + 23.0) / 7.0));
        tau_s_fib = 1574. + 5268. * exp(-(pow((yf[x][0]  + 23. ) / 22.7, 2 )));
        //s_gate_fib = s_inf_fib - (s_inf_fib - s) * exp(-DT / tau_s_fib);
        yfdot[x][2] = (s_inf_fib - yf[x][2])/tau_s_fib;

        // Current injection
        if ( (( count <= stim_dur_int)&& (stim_fib[0][x] == 1))&& num>2) {
          I_stim_fib = -stim_mag;
        }
        else {
          I_stim_fib = 0.0;
        }

        yfdot[x][0] = (-Itotal_fib);

        //Update membrane potential and gates
        for (int cc=0; cc<3; cc++) yf[x][cc] = yf[x][cc] + DT * yfdot[x][cc];

        if(count%100==0) fprintf(output2,"%.12f\t",yf[x][0]);

      } //end of cell loop fibroblast
      ////////////////////////////////////////////////////////////////////////////////////////////
      ////END OF CELL LOOP
      // for (int cc=35; cc<=41; cc++) fprintf(output,"%.12f\t",y[x][cc]);
      if(count%100==0) fprintf(output2,"\n");


        //temp storage of 1 d array
  double temp[CELLS_myo], temp2[CELLS_fib], temp3[CELLS_myo], temp4[CELLS_fib];
 // cout<<"1:"<<temp3[0]<<endl;
  // set the contents of buf equal to the contents of the first row of multiarray.
  memcpy(temp, ydot[0], CELLS_myo * sizeof(double)); 
  memcpy(temp2, yfdot[0],CELLS_fib * sizeof(double)); 

  memcpy(temp3, y[0], CELLS_myo * sizeof(double)); 
  memcpy(temp4, yf[0],CELLS_fib * sizeof(double)); 
 // cout<<"2:"<<temp3[0]<<endl;

  //Solve PDE for all neighbors
  PDE(temp,temp3, neighbors_myo);

  PDE(temp2,temp4, neighbors_fib);
 // cout<<"3:"<<temp3[0]<<endl;

  //cout <<y[0][0]<<endl;

    } //end of bcl loop


  }// end of stimulation loop




  //for (int bb=0; bb<=41; bb++) printf("y[x][%d] = %.12f;\n",bb,y[x][bb]);

  fclose(output);
  fclose(output2);

  return 0;
  }
  void PDE( double *I, double *V, vector<vector<int> > & neighbor ){
    double G;
    for (int p = 0; p < neighbor.size(); p++){
      for (int N = 1; N < neighbor[p].size(); N++){//starts at 1 because first entry (0) is self
        if(neighbor[p][0] > 0){//If self is a myocyte, determine neighbor
          if(neighbor[p][N] > 0){ //If neighbor is a myocyte
            if((neighbor[p][N] == p+2)||(neighbor[p][N]== p-0)){//lazy way of indicating end to end connections
              G = 600.;
              //cout << "end to end"<<endl;
            }
            else{
              G = 200.;  //If the neighbor is end to end then conductance should be 600
           }
          }
          else{ //If neighbor is a fibroblast
            G = 3.;
          //   cout << "3" <<endl;

          }
        }
        else {//If self is a fibroblast is is alway small conductance
          G = 3.;

        }
        I[p] += G *DT* (V[ abs(neighbor[p][N]) ] - V[p]);

      }
    }
  }
  vector<vector<int> > fileToVector(const char *name){
    vector<vector<int> > result;
    ifstream input (name);
    string lineData;
    //  cout<<"File empty"<<endl;

    while(getline(input, lineData))
    {
      if(lineData.empty()){
        cout<<"File empty"<<endl;
      }
      else{
        // cout<<"File empty"<<endl;
        double d;
        vector<int> row;
        stringstream lineStream(lineData);

        while (lineStream >> d){
          row.push_back(d);
          //  if( !lineStream.eof() ) cout << "Dammit!"<<endl;
        }

        result.push_back(row);
      }
    }

    return result;
  }

