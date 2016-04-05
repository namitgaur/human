#include <fstream>
#include "MersenneTwister.h"
#include <omp.h>

using namespace std;

#define Nryr 100	//Number of RyR2s in a dyad 
#define Ndhpr 15	//Number of LCCs in a dyad
#define Ntotdyads 10000	//Number of dyads in cell
#define bcl 500		//Basic Cycle Length (ms)
#define beats 107		//Number of beats

#define Nx 100		//Cell dimension x
#define Ny 1		//Cell dimension y
#define Nz 1		//Cell dimension z
#define Nvox 1		//Voxel size
#define Ndyads 100	//Number of dyads employed in the simulation

/* Terms for solution of conductance and Reversal Potential */
double R; //Universal Gas constant (J/mol.K)
double temp; //Temperature (K)
double frdy; //Faraday's Constant (C/mol)

/* Cell Geometry */
double ar; //Radius of the cell (cm)
double al; //Length of the cell (cm)
double pi; //Pi
double ageo; //Geometric membrane area (cm^2)
double acap; //Capacitive membrane area (cm^2)
double vcell; //Cell volume (uL)
double vmyo; //Myoplasmic volume (uL)
double vss;  //Submembrane space volume (uL)
double vnsr; //NSR volume (uL)
double vjsr; //JSR volume (uL)

/* Dyadic Geometry */
double vcyt_dyad; //Volume of local myoplasm (uL)
double vss_dyad; //Volume of local submembrane space (uL)
double vnsr_dyad; //Volume of local NSR (uL)
double vjsr_dyad; //Volume of local JSR (uL)
double vds[Ndyads], vol_ds; //Volume of dyadic space (uL)

/* Time Step Counter */
int idum; // Seed variable for random number generator
int i; // Loop counter over dyads
int k; // Loop variable over channels (RyR2s and LCCs) in a dyad
int icounter;
int steps; // Number of steps
int ibeat; // Stimulation counter
double udt; // Universal Time Step (ms)
double dt; // Time Step (ms)
double dx; // Space step (um)
double dy; // Space step (um)
double dz; // Space step (um)
double dtdiff;
double t; // Time (ms)
double to;// Time during which RyR is open (ms)
double tstim; // Time stimulation is applied (ms)
double stimtime; // Time period during which stimulation is applied

/* Creation of Data File */
double printval;
double printdata;
ofstream f0out ("calcium"); 	 // Output Data file for whole-cell variables
ofstream f1out ("currents");
ofstream f2out ("dyad");    // Output Data file for dyadic variables
ofstream f3out ("cacyt");   // Output Data file for 1D Cai plot (along z axis)
ofstream f4out ("cajsr");   // Output Data file for 1D CaSR plot (along z axis)
ofstream f5out ("initial_conditions");

/* Ion valences */
double zna=1; // Na valence
double zk = 1; // K valence
double zca =2; // Ca valence

/* Free Ion concentrations */
double nai; // Myoplasmic Na concentration (mM)
double dnai; // Change in myoplasmic Na concentration (mM)
double nao; // Extracellular Na concentration (mM)
double ki; // Myoplasmic K concentration (mM)
double ko; // Extracellular K concentration (mM)
double dki; // Change in myoplasmic K concentration (mM)
double cads[Ndyads]; // Dyadic space Ca concentration (uM)
double cass[Ndyads]; // local Submembrane space Ca concentration (mM)
double cansr[Ndyads]; // local NSR concentration (mM)
double cajsr[Ndyads]; // local JSR concentration (mM)
double cacyt[Ndyads]; // local Myoplasmic Ca concentration (mM)
double cao; // Extracellular Ca concentration (mM)
double cai; // Whole-cell myoplasmic Ca concentration (mM)
double casub; // Whole-cell submembrane space Ca concentration (mM)
double cadyad; // Whole-cell (Avg.) dyadic space Ca concentration (uM)
double nsr; // Whole-cell NSR Ca concentration (mM)
double jsr; // Whole-cell JSR Ca concentration (mM)

int Nx_diff = (Nx-1)*Nvox + 1;
int Ny_diff = (Ny-1)*Nvox + 1;
int Nz_diff = (Nz-1)*Nvox + 1;

double casstemp[((Nx-1)*Nvox+1)*((Ny-1)*Nvox+1)*((Nz-1)*Nvox+1)];
double cansrtemp[((Nx-1)*Nvox+1)*((Ny-1)*Nvox+1)*((Nz-1)*Nvox+1)];
double cacyttemp[((Nx-1)*Nvox+1)*((Ny-1)*Nvox+1)*((Nz-1)*Nvox+1)];
double cassdiff[((Nx-1)*Nvox+1)*((Ny-1)*Nvox+1)*((Nz-1)*Nvox+1)];
double cansrdiff[((Nx-1)*Nvox+1)*((Ny-1)*Nvox+1)*((Nz-1)*Nvox+1)];
double cacytdiff[((Nx-1)*Nvox+1)*((Ny-1)*Nvox+1)*((Nz-1)*Nvox+1)];

int iso;

/* Voltage */
double v;

/* Fast Na Current */
double ina; // Fast Na current (uA/uF)
double gna; // Max. Conductance of the Na Channel (mS/uF)
double ena; // Reversal Potential of Na (mV)
double am; // Na alpha-m rate constant (ms^-1)
double bm; // Na beta-m rate constant (ms^-1)
double ah; // Na alpha-h rate constant (ms^-1)
double bh; // Na beta-h rate constant (ms^-1)
double aj; // Na alpha-j rate constant (ms^-1)
double bj; // Na beta-j rate constant (ms^-1)
double mtau; // Na activation time constant
double htau; // Na inactivation time constant
double hftau;
double hstau;
double hsptau;
double jtau; // Na slow inactivation time constant
double mss; // Na steady state activation gate value
double hss; // Na steady state inactivation gate value
double jss; // Na steady state slow inactivation gate value
double m; // Na activation gate
double h; // Na inactivation gate
double j; // Na slow inactivation gate
double Ahf;
double Ahs;
double hf;
double hs;
double hsp;
double hssp;
double hp;
double jp;
double hptau;
double jptau;
double finap;

/* Late Na Current */
double inal;
double gnal;
double mltau;
double mlss;
double ml;
double hltau;
double hlss;
double hl;
double hlptau;
double hlssp;
double hlssptau;
double hlp;
double finalp;

/* Transient outward K current, Ito */
double a,ass,atau;
double i_ito,iftau,istau,delta_epi,AiF,AiS,iF,iS,iss;
double assp,ap,dti_develop,dti_recover,tiFp,tiSp,tiF,tiS,iFp,iSp,ip,fItop;
double Gto,Ito,EK;

/* Background Na current */
double INab, IKb;

/* Slowly Activating Potassium Current */
double xs1,xs2,IKs;

// Rapidly Activating Potassium Current */
double xrf,xrs,IKr;
/* Potassium Current (time-independent) */
double xk1, IK1;
/* Plateau Potassium Current */
double ikp; // Plateau K current (uA/uF)
double gkp; // Channel Conductance of Plateau K Current (mS/uF)
double ekp; // Reversal Potential of Plateau K Current (mV)
double kp; // K plateau gate

/* Sodium-Potassium Pump */
double INaK;

/* Current through L-type Ca Channel */
double pna; // Permeability of membrane to Na (cm/s)
double ganai; // Activity coefficient of Na
double ganao; // Activity coefficient of Na
double pk; // Permeability of membrane to K (cm/s)
double gaki; // Activity coefficient of K
double gako; // Activity coefficient of K
double ibarna; // Max. Na current through Ca channel (uA/uF)
double ibark; // Max. K current through K channel (uA/uF)
double ilcana; // Na current through L-type channel (uA/uF)
double ilcak; // K current through L-type channel (uA/uF)
double ilcatot; // Total current through the L-type channel (uA/uF)
double nc0tot; // Number of channels in C0 state in a dyad
double nc1tot; // Number of channels in C1 state in a dyad
double nc2tot; // Number of channels in C2 state in a dyad
double nc3tot; // Number of channels in C3 state in a dyad
double notot; // Number of channels in O state in a dyad
double nmodeca; // Number of channels in Mode ca in a dyad
int nLCC_dyad_p;
int nLCC_dyad_np;
double nodhpr; // Total number of LCC in O state in the cell
double nmodeca_total; // Total number of LCC in Mode Ca in the cell
double randomnumber; //Random Number used to determine channel state in Monte Carlo simulations
double pca; // Permeability of membrane to Ca (cm/s)
double gacai; // Activity coefficient of Ca
double gacao; // Activity coefficient of Ca
double ibarca; // Max. current through Ca channel (uA/uF)
double ilca; // Ca current through Ca channel (uA/uF)
double ilca_dyad; // Ca current through Ca channel in dyad (uA/uF)
double jlca_dyad[Ndyads]; // Ca flux through Ca channel in dyad (uM/ms)
double d,f,fca,fp,fcap,nca[Ndyads],nca_print,icana_dyad,icak_dyad;
double ff,fs,fcaf,fcas,jca,ffp,fcafp;
double ICaL_print, ICaNa_print, ICaK_print;

/* Total current and stimulus */
double ist; // Stimulus Current (uA/cm^2)
double it; // Total current (uA/cm^2)
double naiont; // Total Na current
double kiont; // Total K current
double caiont; // Total Ca current

/* Ca Cycling Process */
double Km; // Km of activation of RyR2
double Kmcsqn; // Km of activation of RyR2 when inhibited by CSQN
double dryr; // Diffusion rate of Ca through open RyR2 (ms^-1)
double kC1O1; // Transition rate from C1 to O1 of RyR2 (ms^-1)
double kO1C1; // Transition rate from O1 to C1 of RyR2 (ms^-1)
double kC1C2; // Transition rate from C1 to C2 of RyR2 (ms^-1)
double kC2C1; // Transition rate from C2 to C1 of RyR2 (ms^-1)
double kO1O2; // Transition rate from O1 to O2 of RyR2 (ms^-1)
double kO2O1; // Transition rate from O2 to O1 of RyR2 (ms^-1)
double kC2O2; // Transition rate from C2 to O2 of RyR2 (ms^-1)
double kO2C2; // Transition rate from O2 to C2 of RyR2 (ms^-1)
double csqn; // Ca-unbound CSQN (mM)
double csqnbar; // Total concentration of CSQN (mM)
double csqnca; // Ca-bound CSQN (mM)
double kmcsqn; // Half-saturation constant of Ca binding of CSQN (mM)
double trel[Ndyads]; // Activation time of the dyad
double jrel; // SR Ca release flux (mM/s)
double jrel_dyad[Ndyads]; // SR Ca release flux in the dyad (uM/ms)
double nryropen; // Number of open RyR2s in a dyad
int ryrstate[Ndyads][Nryr]; // State of individual RyR2
double nryrc1; // Number of channels in C1 state in a dyad
double nryro1; // Number of channels in O1 state in a dyad
double nryrc2; // Number of channels in C2 state in a dyad
double nryro2; // Number of channels in O2 state in a dyad
double noryr; // Number of open RyR2s in the cell
double tauefflux; // Time constant of efflux from dyadic space to submembrane space
double nactive[Ndyads]; // Activity status of a dyad (active/passive)
double factivedyad; // Fraction of dyads that are active
double jrefill; //local flux from NSR_JSR (uM/ms)
double taurefill[Ndyads]; //Time constant of flux from local NSR to local JSR (ms)
double bjsr; // Buffering factor in JSR
double iup; // SR Uptake flux (mM/ms)
double iupbar; // Max. flux through SR pump (mM/ms)
double kmup; // Half-saturation concentration of iup (mM)
double ileak; // Ca leakage from NSR to myoplasm (mM/ms)
double nsrbar; // Max. Ca in NSR (mM)
double trpnbar; // Max. Ca buffered in troponin (TRPN) (mM)
double cmdnbar; // Max. Ca buffered in calmodulin (mM)
double kmtrpn; // Half-saturation constant of Ca binding to TRPN
double kmcmdn; // Half-saturation constant of Ca binding to CMDN
double bmyo; // Buffering factor in myoplasm

double jefflux; // Ca flux from dyadic space to submembrane space (uM/ms)
 
double bslbar; // Max. conc. of SL buffers (mM)
double bsrbar; // Max. conc. of SR buffers (mM)
double kmbsl; // Half-saturation constant of Ca binding to SL buffers (mM)
double kmbsr; // Half-saturation constant of Ca binding to SR buffers (mM)

double taudiff; // Time constant of Ca diffusion between submembrane space and myoplasm 
double taudiff_nsr_nsr; // Time constant of Ca diffusion between local NSRs (ms)
double taudiff_ds_ds; // Time constant of Ca diffusion (ms)
double jcadiff_nsr_nsr; //  Diffusion flux between local NSR compartments (uM/ms)
double jcadiff_ds_ds[Ndyads]; // Diffusion flux between dyads (uM/ms)
double jcadiff_cass_cass;	//Diffusion between subsarcolemmal space compartments
double jcadiff_cacyt_cacyt;	//Diffusion between myoplasmic compartments
double Dsr;
double Dcyt;

double bsubspace; // Buffering factor in submembrane space
double jcadiff_ds_ss; // Diffusion flux from dyadic space to submembrane space (uM/ms)
double jcadiff_ss_cyt; // Diffusion flux from submembrane space to myoplasm (uM/ms)

/* Na-Ca exchanger current */
double inaca; // Na-Ca exchanger current in myoplasm (uA/uF)
double inacass; // Na-Ca exchanger current in submembrane space (uA/uF)
double c1; // Scaling factor for inaca
double c2; // Half-saturation concentration of NaCa exchanger (mM)
double c3; // Scaling factor for inaca
double gammas; // Position of energy barrier controlling voltage dependence of inaca
double inaca_dyad; // local Na-Ca exchanger current in myoplasm
double inacass_dyad; // local Na-Ca exchanger current in submembrane space

/* Sarcolemmal Ca Pump */
double ipca,ipca_dyad,ibarpca,kmpca;

/* T-type Ca current */
double icat; // Ca current through T-type Ca channel (uA/uF)
double icat_dyad; // local Ca current through T-type Ca channel (uA/uF)
double b; // activation gate of T-type channel 
double bss; // Steady-state value of activation gate b
double taub; // Time constant of gate b (ms^-1)
double g; // Inactivation gate
double gss; // Steady-state value of inactivation gate g
double taug; // Time-constant of gate g (ms^-1)
double eca; // Reversal potential of T-type Ca channel (mV)
double gcat; // Max. conductance of the T-type Ca channel (mS/uF)

/* Background Ca current */
double icab; // Background Ca current (uA/uF)
double icab_dyad; // local background Ca current (uA/uF)
double gcab; // Max. conductance of Ca background (mS/uF)
double ecan; // Nernst Potential for Ca (mV);

//CaMK paramaters
double const aCaMK=0.05;
double const bCaMK=0.00068;
double const CaMKo=0.05;
double const KmCaM=0.0015;
double const KmCaMK=0.15;
double CaMKa[Ndyads],CaMKb[Ndyads],CaMKt[Ndyads];
double CaMKa_cell,CaMKb_cell, CaMKt_cell;

MTRand mtrand1;

/* Functions */
void comp_ina(); // Computes fast Na current
void comp_ito(); // Computes transient outward K current
void comp_inab(); // Computes background Na current
void comp_iks(); // Computes slowly activating K current
void comp_ikr(); // Computes rapidly activating K current
void comp_ik1(); // Computes Time-independent K current
void comp_ikp(); // Computes plateau K current
void comp_it(); // Computes total current and voltage
void comp_ical(); // Computes L-type Ca current
void comp_icalhh();
void comp_inak(); // Computes Na-K pump current
void calc_nai(); // Calculates myoplasmic Na concentration
void calc_ki(); // Calculates myoplasmic K concentration
void comp_calc_dyad();  //Dyadic calculation function//Random Number Generator
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-12
#define RNMX (1.0-EPS)

struct RandData
{
        long myrandomseed;
        long iy;
        long iv[NTAB];
        double padding[16];
};

double ndrand48(struct RandData * pRD)
{
        int j;
        long k;
        double temp;
        if (pRD->myrandomseed <= 0 || ! pRD->iy)
	{
                if (-(pRD->myrandomseed) < 1)
                        pRD->myrandomseed=1;
                else
                        pRD->myrandomseed = -(pRD->myrandomseed);
                for (j=NTAB+7;j>=0;j--) 
		{
                        k=(pRD->myrandomseed)/IQ;
                        pRD->myrandomseed=IA*(pRD->myrandomseed-k*IQ)-IR*k;
                        if (pRD->myrandomseed < 0)
                                pRD->myrandomseed += IM;
                        if (j < NTAB)
                                pRD->iv[j] = pRD->myrandomseed;
                }
                pRD->iy=pRD->iv[0];
        }

        k=(pRD->myrandomseed)/IQ;
        pRD->myrandomseed=IA*(pRD->myrandomseed-k*IQ)-IR*k;
        if (pRD->myrandomseed < 0)
                pRD->myrandomseed += IM;
        j=pRD->iy/NDIV;
        pRD->iy=pRD->iv[j];
        pRD->iv[j] = pRD->myrandomseed;
        if ((temp=AM*pRD->iy) > RNMX)
                return RNMX;
        else
                return temp;
}

struct RandData *pRD = (struct RandData *)malloc(sizeof(struct RandData));
void calc_conc_dyad(int i); // Calculates concentrations in the dyad
void comp_ical_dyad(int i); // Computes L-type current
void comp_irel_dyad(int i); // Computes SR Ca release
void comp_ica_dyad(); // Computes local currents in a dyad
void ca_conc_diff_update();
void prttofile(); // Data File Output function
void prtfinal();

int main(void)
{
	trpnbar = 0.07;
	kmtrpn = 0.0005;
	cmdnbar = 0.05;
	kmcmdn = 0.00238;	
	csqnbar = 10;
	kmcsqn = 0.5;
	taudiff_ds_ds = 200;

	//Diffusion parameters
	Dsr = 0.02;
	Dcyt = 0.04; Dcyt = 0.02;
	R = 8314;
	temp = 310;
	frdy = 96485;
	ar = 0.0011;
	al = 0.01;
	pi = 3.142;
	ageo = 2*pi*ar*ar + 2*pi*ar*al;
	acap = 2*ageo;
	vcell = 1000*pi*ar*ar*al;	
	vmyo = 0.68*vcell;
	vnsr = 0.0552*vcell;
	vjsr = 0.0048*vcell;
	vss = 0.01*vcell;

	vjsr_dyad = vjsr/Ntotdyads;
	vnsr_dyad = vnsr/Ntotdyads;
	vcyt_dyad = vmyo/Ntotdyads;
	vss_dyad = vss/Ntotdyads;
	vol_ds = 1e-13;
	idum = -2;
	t = 0;
	udt = 0.05;
	dt = udt;
	dx = 0.2;
	dy = 0.2;
	dz = 0.2;
	steps = (bcl*beats)/udt;

	printval = 1;
	iso = 1;

	if (iso == 1) kmtrpn = 1.6*kmtrpn;

	cao = 1.8;
	nao = 140;
	ko = 5.4;
	
	t = 0;
	ibeat = 0;
	tstim = 50; // time of stimulus
	stimtime = 100;
		
	pRD->myrandomseed = 5;
	pRD->iy = -1;

	//Initial Conditions
v=-87.4166;
nai=5.88383;
ki=141.712;
cai=9.35235e-05;
casub=9.42268e-05;
nsr=4.27252;
jsr=3.81969;
m=0.0211631;
hf=0.480144;
hs=0.479739;
j=0.473817;
hsp=0.249023;
jp=0.460561;
ml=0.000210404;
hl=0.296362;
hlp=0.134898;
a=0.00104148;
iF=0.999505;
iS=0.20048;
ap=0.000530675;
iFp=0.999505;
iSp=0.221394;
d=7.36138e-08;
ff=1;
fs=0.692324;
fcaf=1;
fcas=0.933007;
jca=0.955313;
ffp=0.999976;
fcafp=0.999998;
xrf=0.0202163;
xrs=0.821959;
xs1=0.524584;
xs2=0.000847115;
xk1=0.996917;






	csqn = (kmcsqn/(kmcsqn+jsr))*csqnbar;

	for (i=0;i<Ndyads;i++)
	{
		cads[i] = casub*1e3;
		cass[i] = casub;
		cacyt[i] = cai;
		cajsr[i] = jsr;
		cansr[i] = nsr;
		trel[i] = -2000;
		CaMKa[i] = CaMKb[i] = CaMKt[i] = 0.1;

		kC1C2 = csqn/csqnbar;
		kC2C1 = 1-csqn/csqnbar;

		for (k = 0; k<Nryr; k++)
		{
			if ( mtrand1()< kC2C1/(kC2C1+kC1C2) ) ryrstate[i][k] = 0; else ryrstate[i][k] = 2;
		}

		vds[i] = vol_ds;
	}
	for (i = 0; i<Nx_diff*Ny_diff*Nz_diff; i++)
	{
		cassdiff[i] = casub;
		cansrdiff[i] = nsr;
		cacytdiff[i] = cai;

	}
	// Time Loop		
//	#pragma omp parallel			
	for (icounter = 0; icounter<steps; icounter++)
	{
		comp_ina();
                comp_ito();
		comp_ikr();
		comp_iks();
		comp_ik1();
		comp_inak();
		comp_inab();
		comp_calc_dyad();
	 	comp_ical();	
		comp_it();
		calc_nai();
		calc_ki();

		stimtime = stimtime+dt;
		v = v-it*dt;
		if (t>(beats-10)*bcl)
		{
			if (icounter % 1 == 0)	
			{	 
				prttofile();
				printdata = 0;
			}
		}
		printdata = printdata+1;

		t = t+udt;
	}
	
	//print final conditions
	prtfinal();

}

void comp_ina ()
{	
	gna = 75;
	if (iso == 1) gna = 2.7*gna;
	
	ena = ((R*temp)/frdy)*log(nao/nai);

	mtau = 1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
	hftau=1.0/(1.432e-5*exp(-(v+1.196)/6.285)+6.149*exp((v+0.5096)/20.27));
	hstau=1.0/(0.009794*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
	hsptau=3.0*hstau;

	jtau=2.038+1.0/(0.02136*exp(-(v+100.6)/8.281)+0.3052*exp((v+0.9941)/38.45));
	jptau=1.46*jtau;
	
	Ahf=0.99;
	Ahs=1.0-Ahf;

	mss =1.0/(1.0+exp((-(v+49.57))/9.871));

	hss =1.0/(1+exp((v+82.90)/6.086));
	hssp=1.0/(1+exp((v+89.1)/6.086));

	if (iso == 1)
	{		
		hss =1.0/(1+exp((v+82.90+5)/6.086));             		
		hssp=1.0/(1+exp((v+89.1+5)/6.086));
	}

	jss = hss;	
	m = mss-(mss-m)*exp(-dt/mtau);
	hf=hss-(hss-hf)*exp(-dt/hftau);
	hs=hss-(hss-hs)*exp(-dt/hstau);
	hsp=hssp-(hssp-hsp)*exp(-dt/hsptau);

	h = Ahf*hf+Ahs*hs;
	hp=Ahf*hf+Ahs*hsp;
	j = jss-(jss-j)*exp(-dt/jtau);
	jp=jss-(jss-jp)*exp(-dt/jptau);
	
	finap=(1.0/(1.0+KmCaMK/CaMKa_cell));

	mlss=1.0/(1.0+exp((-(v+42.85))/5.264));
	mltau = mtau;    
	ml=mlss-(mlss-ml)*exp(-dt/mltau);
    
	hlss=1.0/(1.0+exp((v+87.61)/7.488));
	hlssp=1.0/(1.0+exp((v+93.81)/7.488));
	hltau=200.0;
	hlptau=3.0*hltau;
	hl=hlss-(hlss-hl)*exp(-dt/hltau);
	hlp=hlssp-(hlssp-hlp)*exp(-dt/hlptau);
    
	gnal = 0.0075;
	finalp = finap;
    	
	ina=gna*(v-ena)*m*m*m*((1.0-finap)*h*j+finap*hp*jp);
	inal=gnal*(v-ena)*ml*((1.0-finalp)*hl+finalp*hlp);	
}

void comp_ito()
{
     EK=(R*temp/frdy)*log(ko/ki);
     
     ass=1.0/(1.0+exp((-(v-14.34))/14.82));
	atau=1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814))); 
	a=ass-(ass-a)*exp(-dt/atau);
     
     iss=1.0/(1.0+exp((v+43.94)/5.711));
     delta_epi=1.0;
     
     iftau=4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
     istau=23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
     
     iftau*=delta_epi;
     istau*=delta_epi;
     
     AiF=1.0/(1.0+exp((v-213.6)/151.2));
     AiS=1.0-AiF;
     
     iF=iss-(iss-iF)*exp(-dt/iftau);
     iS=iss-(iss-iS)*exp(-dt/istau);
     
     i_ito=AiF*iF+AiS*iS;
     
     assp=1.0/(1.0+exp((-(v-24.34))/14.82));
     ap=assp-(assp-ap)*exp(-dt/atau);
     
     dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
     dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0));
     
     tiFp=dti_develop*dti_recover*iftau;
     tiSp=dti_develop*dti_recover*istau;
     
     iFp=iss-(iss-iFp)*exp(-dt/tiFp);
     iSp=iss-(iss-iSp)*exp(-dt/tiSp);
     
     ip=AiF*iFp+AiS*iSp;
     
     Gto=0.02;
     fItop=(1.0/(1.0+KmCaMK/CaMKa_cell));

     Ito=Gto*(v-EK)*((1.0-fItop)*a*i_ito+fItop*ap*ip);
}

void comp_inab ()
{
double xkb=1.0/(1.0+exp(-(v-14.48)/18.34));
double GKb=0.003;

IKb=GKb*xkb*(v-EK);

double PNab=3.75e-10;
double vffrt = v*frdy*frdy/(R*temp);
double vfrt = v*frdy/(R*temp);
INab=PNab*vffrt*(nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);
}

void comp_iks ()
{
double EKs=(R*temp/frdy)*log((ko+0.01833*nao)/(ki+0.01833*nai));
double xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
double txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
xs1=xs1ss-(xs1ss-xs1)*exp(-dt/txs1);
double xs2ss=xs1ss;
double txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
xs2=xs2ss-(xs2ss-xs2)*exp(-dt/txs2);
double KsCa=1.0+0.6/(1.0+pow(3.8e-5/3.8e-5,1.4));
double GKs=0.0034;
if (iso == 1) GKs = 3.2*GKs;

IKs=GKs*KsCa*xs1*xs2*(v-EKs);
}

void comp_ikr ()
{
double xrss=1.0/(1.0+exp((-(v+8.337))/6.789));
double txrf=12.98+1.0/(0.3652*exp((v-31.66)/3.869)+4.123e-5*exp((-(v-47.78))/20.38));
double txrs=1.865+1.0/(0.06629*exp((v-34.70)/7.355)+1.128e-5*exp((-(v-29.74))/25.94));
double Axrf=1.0/(1.0+exp((v+54.81)/38.21));
double Axrs=1.0-Axrf;
xrf=xrss-(xrss-xrf)*exp(-dt/txrf);
xrs=xrss-(xrss-xrs)*exp(-dt/txrs);
double xr=Axrf*xrf+Axrs*xrs;
double rkr=1.0/(1.0+exp((v+55.0)/75.0))*1.0/(1.0+exp((v-10.0)/30.0));
double GKr=0.046; GKr = 1*GKr;

// if (t>(beats-17)*bcl && t<(beats-14)*bcl) GKr = 0.3*GKr;

IKr=GKr*sqrt(ko/5.4)*xr*rkr*(v-EK);
}

void comp_ik1 ()
{
double xk1ss=1.0/(1.0+exp(-(v+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
double txk1=122.2/(exp((-(v+127.2))/20.36)+exp((v+236.8)/69.33));
xk1=xk1ss-(xk1ss-xk1)*exp(-dt/txk1);
double rk1=1.0/(1.0+exp((v+105.8-2.6*ko)/9.493));
double GK1=0.1908;
GK1 = 0.4*0.1908;

IK1=GK1*sqrt(ko)*rk1*xk1*(v-EK);

}

void comp_ical()
{
double dss=1.0/(1.0+exp((-(v+3.940))/4.230));
if (iso == 1) dss=1.0/(1.0+exp((-(v+3.940+14))/4.230));

double td=0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
d=dss-(dss-d)*exp(-dt/td);

double fss=1.0/(1.0+exp((v+19.58)/3.696));
if (iso == 1) fss=1.0/(1.0+exp((v+19.58+8)/3.696));

double tff=7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
double tfs=1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
double Aff=0.6;
double Afs=1.0-Aff;
ff=fss-(fss-ff)*exp(-dt/tff);
fs=fss-(fss-fs)*exp(-dt/tfs);
f=Aff*ff+Afs*fs;
double fcass=fss;
double tfcaf=7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
double tfcas=100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));
double Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0));
double Afcas=1.0-Afcaf;
fcaf=fcass-(fcass-fcaf)*exp(-dt/tfcaf);
fcas=fcass-(fcass-fcas)*exp(-dt/tfcas);
fca=Afcaf*fcaf+Afcas*fcas;
double tjca=75.0;
jca=fcass-(fcass-jca)*exp(-dt/tjca);
double tffp=2.5*tff;
ffp=fss-(fss-ffp)*exp(-dt/tffp);
fp=Aff*ffp+Afs*fs;
double tfcafp=2.5*tfcaf;
fcafp=fcass-(fcass-fcafp)*exp(-dt/tfcafp);
fcap=Afcaf*fcafp+Afcas*fcas;

//for plotting ICaL only//
double vffrt = v*frdy*frdy/(R*temp);
double vfrt = v*frdy/(R*temp);
double Kmn=0.002;
double k2n=1000.0;
double km2n=jca*1.0;
double anca=1.0/(k2n/km2n+pow(1.0+Kmn/cadyad,4.0));
nca_print=anca*k2n/km2n-(anca*k2n/km2n-nca_print)*exp(-km2n*dt);
double fICaLp=(1.0/(1.0+KmCaMK/CaMKa_cell));
double PCa=0.0001;
if (iso == 1) PCa = 2.5*PCa;

PCa = PCa;  //ICaL_Block

double PCap=1.1*PCa;
double PCaNa=0.00125*PCa;
double PCaK=3.574e-4*PCa;
double PCaNap=0.00125*PCap;
double PCaKp=3.574e-4*PCap;
double PhiCaL=4.0*vffrt*(1e-3*cadyad*0.01*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
double PhiCaNa=1.0*vffrt*(0.75*nai*exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
double PhiCaK=1.0*vffrt*(0.75*ki*exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);

ICaL_print=(1.0-fICaLp)*PCa*PhiCaL*d*(f*(1.0-nca_print)+jca*fca*nca_print)+fICaLp*PCap*PhiCaL*d*(fp*(1.0-nca_print)+jca*fcap*nca_print);
ICaNa_print=(1.0-fICaLp)*PCaNa*PhiCaNa*d*(f*(1.0-nca_print)+jca*fca*nca_print)+fICaLp*PCaNap*PhiCaNa*d*(fp*(1.0-nca_print)+jca*fcap*nca_print);
ICaK_print=(1.0-fICaLp)*PCaK*PhiCaK*d*(f*(1.0-nca_print)+jca*fca*nca_print)+fICaLp*PCaKp*PhiCaK*d*(fp*(1.0-nca_print)+jca*fcap*nca_print);


}

void comp_inak()
{
double k1p=949.5;
double k1m=182.4;
double k2p=687.2;
double k2m=39.4;
double k3p=1899.0;
double k3m=79300.0;
double k4p=639.0;
double k4m=40.0;
double Knai0=9.073;
double Knao0=27.78;
double delta=-0.1550;
double Knai=Knai0*exp((delta*v*frdy)/(3.0*R*temp));
if (iso == 1) Knai = 0.7*Knai;
double Knao=Knao0*exp(((1.0-delta)*v*frdy)/(3.0*R*temp));
double Kki=0.5;
double Kko=0.3582;
double MgADP=0.05;
double MgATP=9.8;
double Kmgatp=1.698e-7;
double H=1.0e-7;
double eP=4.2;
double Khp=1.698e-7;
double Knap=224.0;
double Kxkur=292.0;
double P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
double a1=(k1p*pow(nai/Knai,3.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
double b1=k1m*MgADP;
double a2=k2p;
double b2=(k2m*pow(nao/Knao,3.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
double a3=(k3p*pow(ko/Kko,2.0))/(pow(1.0+nao/Knao,3.0)+pow(1.0+ko/Kko,2.0)-1.0);
double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
double b4=(k4m*pow(ki/Kki,2.0))/(pow(1.0+nai/Knai,3.0)+pow(1.0+ki/Kki,2.0)-1.0);
double x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
double x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
double x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
double x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
double E1=x1/(x1+x2+x3+x4);
double E2=x2/(x1+x2+x3+x4);
double E3=x3/(x1+x2+x3+x4);
double E4=x4/(x1+x2+x3+x4);
double zk=1.0;
double JnakNa=3.0*(E1*a3-E2*b3);
double JnakK=2.0*(E4*b1-E3*a1);
double Pnak=30;

INaK=Pnak*(zna*JnakNa+zk*JnakK);
}

void comp_it()
{
	if (t>=tstim && t<(tstim+dt))
	{
		stimtime = 0;
		ibeat = ibeat+1;
		tstim = tstim+bcl;
	
		cout<<ibeat<<endl;
	}
	
	if (  (stimtime>=0 && stimtime<=0.5) )
		ist = -80;
	else
		ist = 0;

	if (ibeat>(beats-7)) ist = 0;

	naiont = ina+inal+INab+ilcana+3*INaK+3*inaca+3*inacass;
	kiont = Ito + IKr + IKs + IK1 + IKb + ilcak - 2*INaK + ist;
	caiont = ilca+icab+ipca-2*inaca-2*inacass;		

	it = naiont+kiont+caiont;
}

void calc_nai ()
{
	dnai = -dt*(naiont*acap)/(vmyo*zna*frdy);
	nai = dnai + nai;
}

void calc_ki ()
{
	dki = -dt*((kiont)*acap)/(vmyo*zk*frdy);
	ki = dki + ki;
}

void comp_calc_dyad()
{     
	ilca = nodhpr = jrel = noryr = 0;
	ilcana=ilcak=0;
	inaca = icab = ipca = icat = inacass = 0;
	cai = casub = cadyad = jsr = nsr = 0;
	CaMKa_cell = CaMKb_cell = CaMKt_cell = 0;
	factivedyad  = 0;
	
//Compute nca gate

 #pragma omp parallel for
	for ( i=0; i<Ndyads;i++)
	{
		comp_ical_dyad(i); // Computes local L-type current
		comp_irel_dyad(i); // Computes local SR Ca release
		calc_conc_dyad(i); // Computes local Ca concentrations
	}
 #pragma omp single
	ca_conc_diff_update();	// Updates local Ca concentrations due to inter-dyad coupling

	for (i=0; i<Ndyads; i++)
	{
		//Whole-cell concentrations
		cai += cacyt[i]/Ndyads; 
		casub += cass[i]/Ndyads;
		cadyad += cads[i]/Ndyads;
		jsr += cajsr[i]/Ndyads;
		nsr += cansr[i]/Ndyads;
		factivedyad += nactive[i]/Ndyads;
		CaMKa_cell += CaMKa[i]/Ndyads;
	        CaMKb_cell += CaMKb[i]/Ndyads;
	        CaMKt_cell += CaMKt[i]/Ndyads;
	}
}


void comp_ical_dyad(int i)
{    
     //compute nca gate
     double Kmn=0.002;
     double k2n=1000.0;
     double km2n=jca*1.0;
     double anca=1.0/(k2n/km2n+pow(1.0+Kmn/(0.001*cads[i]),4.0));
     nca[i]=anca*k2n/km2n-(anca*k2n/km2n-nca[i])*exp(-km2n*dt);     
     nLCC_dyad_np = 0;
     nLCC_dyad_p = 0;         
        
double a1 = (1-(1.0/(1.0+KmCaMK/CaMKa[i])))*(1-nca[i])*d*f;
double b1 = (1-(1.0/(1.0+KmCaMK/CaMKa[i])))*(nca[i])*d*jca*fca;
double c1 = (1.0/(1.0+KmCaMK/CaMKa[i]))*(1-nca[i])*d*fp;
double d1 =  (1.0/(1.0+KmCaMK/CaMKa[i]))*nca[i]*d*jca*fcap;

      
    for (int k = 0;k<Ndhpr;k++)
    {
        if ( ndrand48(pRD) > (1.0/(1.0+KmCaMK/CaMKa[i])) ) //Channel is not CaMK phosphorylated
        {
           if ( ndrand48(pRD) > nca[i] ) //Channel is not in Ca-CaMK dependent state
           {
              if (ndrand48(pRD) <d*f ) //Channel is Open
              {
                 nLCC_dyad_np ++;
              }                                    
           }           
           else                                    //Channel is in Ca-CaMK dependent state
           {                                                   
              if (ndrand48(pRD) < d*jca*fca)
              {
                 nLCC_dyad_np ++;
              }
           }                                          
        }
        else                                          //Channel is CaMK phosphorylated
        {
           if ( ndrand48(pRD) > nca[i] ) //Channel is not in Ca-CaMK dependent state
           {
              if (ndrand48(pRD) <d*fp ) //Channel is Open
              {
                 nLCC_dyad_p ++;
              }                                    
           }           
           else                                    //Channel is in Ca-CaMK dependent state
           {                                                   
              if (ndrand48(pRD) < d*jca*fcap)
              {
                 nLCC_dyad_p ++;
              }
           }                                          
        }       
    }   



    double vffrt = v*frdy*frdy/(R*temp);
    double vfrt = v*frdy/(R*temp);
    double PhiCaL=4.0*vffrt*(1e-3*cads[i]*0.01*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
    double PhiCaNa=1.0*vffrt*(0.75*nai*exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
    double PhiCaK=1.0*vffrt*(0.75*ki*exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
    double PCa=0.0001;     
    if (iso == 1) PCa = 2.5*PCa;

    PCa = PCa; //ICaL Block
		
    double PCap=1.1*PCa;
    double PCaNa=0.00125*PCa;
    double PCaK=3.574e-4*PCa;
    double PCaNap=0.00125*PCap;
    double PCaKp=3.574e-4*PCap;

    ilca_dyad = PCa*PhiCaL*(nLCC_dyad_np*1.0/Ndhpr) + PCap*PhiCaL*(nLCC_dyad_p*1.0/Ndhpr);	//uA/uF
    icana_dyad = PCaNa*PhiCaNa*(nLCC_dyad_np*1.0/Ndhpr) + PCaNap*PhiCaNa*(nLCC_dyad_p*1.0/Ndhpr);	//uA/uF
    icak_dyad = PCaK*PhiCaK*(nLCC_dyad_np*1.0/Ndhpr) + PCaKp*PhiCaK*(nLCC_dyad_p*1.0/Ndhpr);	//uA/uF

    ilca += ilca_dyad/Ndyads;	//uA/uF
    ilcana += icana_dyad/Ndyads;
    ilcak += icak_dyad/Ndyads;
    nodhpr +=  nLCC_dyad_np + nLCC_dyad_p;

    jlca_dyad[i] = -acap*ilca_dyad*1e3/(2*vds[i]*frdy*Ntotdyads);	//uM/ms
}

void comp_irel_dyad(int i)
{
	Km = 5;		// Cytosolic Ca
	Kmcsqn = 150;
	dryr = 4;

	nryropen = 0;
	nactive[i] = 0;

	for (int k=0; k<Nryr; k++)
	{
		if(ryrstate[i][k] == 1) nryropen++;
	}

	if (  (1.0*nryropen/Nryr)>0.2) {trel[i] = t; nactive[i] = 1;}	
		
	kC1O1 = 3.*(pow(cads[i],4)/(pow(cads[i],4) + pow(Km,4)));
	kO1C1 = 0.5;	
	kC2O2 = 3.*(pow(cads[i],4)/(pow(cads[i],4) + pow(Km,4)));
	kO2C2 = kO1C1;

	csqn = csqnbar*pow(kmcsqn,8)/(pow(cajsr[i],8) + pow(kmcsqn,8) );
	csqnca = csqnbar - csqn;

	kC1C2 = 1*(csqn/csqnbar);

	double tref;
	tref = t-trel[i]-0.95*bcl;
//	tref = 1;

	kC2C1 = 1*(csqnca/csqnbar)*1/(1 + exp(-tref/0.001));
	kO1O2 = kC1C2;
	kO2O1 = kC2C1;

	nryrc1 = nryro1 = nryrc2 = nryro2 = nryropen = 0;

	for ( int k=0;k<Nryr;k++)
	{		
		// Following calculations determine the state of the RyR2 channel
		randomnumber = ndrand48(pRD);
		switch (ryrstate[i][k])
		{
			case 0:
		  		if (randomnumber <= kC1O1*dt)
			    	{ryrstate[i][k] = 1; nryro1++; break;}
				if ((randomnumber>kC1O1*dt) && (randomnumber<= (kC1O1+kC1C2)*dt))
			    	{ryrstate[i][k]=2; nryrc2++;break;}
				if ((randomnumber> (kC1O1+kC1C2)*dt) && (randomnumber <= 1))
			    	{ryrstate[i][k]=0; nryrc1++;break;}
			case 1:
				if (randomnumber <= kO1C1*dt)
				{ryrstate[i][k] = 0; nryrc1++; break;}
				if ((randomnumber>kO1C1*dt) && (randomnumber<=(kO1C1+kO1O2)*dt))
				{ryrstate[i][k] = 3; nryro2++; break;}
				if ((randomnumber>(kO1C1+kO1O2)*dt) && (randomnumber<=1))
				{ryrstate[i][k] = 1; nryro1++; break;}
			case 2:
				if (randomnumber <= kC2C1*dt)
				{ryrstate[i][k] = 0; nryrc1++; break;}
				if ((randomnumber>kC2C1*dt) && (randomnumber<=(kC2C1+kC2O2)*dt))
				{ryrstate[i][k] = 3; nryro2++; break;}
				if ((randomnumber>(kC2C1+kC2O2)*dt) && (randomnumber<=1))
				{ryrstate[i][k] = 2; nryrc2++; break;}
			case 3:
				if (randomnumber <= kO2C2*dt)
				{ryrstate[i][k] = 2; nryrc2++; break;}
				if ((randomnumber>kO2C2*dt) && (randomnumber<=(kO2C2+kO2O1)*dt))
				{ryrstate[i][k] = 1; nryro1++; break;}
				if ((randomnumber>(kO2C2+kO2O1)*dt) && (randomnumber<=1))
				{ryrstate[i][k] = 3; nryro2++; break;}
		}
	}

//	nryropen = nryro1 + nryro2;
	nryropen = nryro1;	

	jrel_dyad[i] = dryr*nryropen*(cajsr[i]*1e3-cads[i]); //uM/ms
	
	jrel += jrel_dyad[i]*(vds[i]/vmyo)*(Ntotdyads/Ndyads); //mM/s

	noryr += nryropen;

}

void calc_conc_dyad(int i)
{    
	tauefflux = 7e-4;
	jcadiff_ds_ss = (cads[i]-cass[i]*1e3)/tauefflux;
	jefflux = jcadiff_ds_ss;

	//Dyadic space concentration
	cads[i] =  (jrel_dyad[i] + jlca_dyad[i] + cass[i]*1e3/tauefflux)*tauefflux;

	bsrbar = 0.047;
	bslbar = 1.124;
	kmbsr = 0.00087;
	kmbsl = 0.0087;
	taudiff = 0.1;
	bsubspace = 1/(1 + ((bsrbar*kmbsr)/(pow((kmbsr+cass[i]),2))) + ((bslbar*kmbsl)/(pow((kmbsl+cass[i]),2))));
	jcadiff_ss_cyt = (cass[i]-cacyt[i])/taudiff;

	// local submembrane space concentration 
	cass[i] = cass[i] + (-bsubspace*(((-2*inacass_dyad)*(acap/(Ntotdyads*zca*frdy*vss_dyad))) + jcadiff_ss_cyt - jcadiff_ds_ss*1e-3*(vds[i]/vss_dyad) ))*dt;


	taurefill[i] = 100;

	jrefill = (cansr[i] - cajsr[i])/taurefill[i];

	bjsr = 1/( 1 + (csqnbar*kmcsqn)/pow((kmcsqn + cajsr[i]),2));

	// local JSR concentration		
	cajsr[i] = cajsr[i] + bjsr*(-(jrel_dyad[i]*(vds[i]/vjsr_dyad))/1000. + jrefill)*dt;
    double Jupnp=0.004375*cacyt[i]/(cacyt[i]+0.00092);
    double Jupp=2.75*0.004375*cacyt[i]/(cacyt[i]+0.00092-0.00017);

  if (iso == 1)
    {
	Jupnp=0.004375*cacyt[i]/(cacyt[i]+0.54*0.00092);
    	Jupp=2.75*0.004375*cacyt[i]/(cacyt[i]+0.54*(0.00092-0.00017));
    }

    double fJupp=(1.0/(1.0+KmCaMK/CaMKa[i]));
    
    

	ileak=0.0039375*cansr[i]/15.0;
	iup=((1.0-fJupp)*Jupnp+fJupp*Jupp);

	// local NSR concentration
	cansr[i] = cansr[i] + (iup - ileak - jrefill*(vjsr_dyad/vnsr_dyad))*dt;	

    // local CaMK bound, active and trapped concentration
    CaMKb[i] = CaMKo*(1.0-CaMKt[i] )/(1.0+KmCaM/(0.001*cads[i]));
    CaMKa[i] = CaMKb[i] + CaMKt[i];
    CaMKt[i] += dt*(aCaMK*CaMKb[i]*(CaMKb[i]+CaMKt[i])-bCaMK*CaMKt[i]);

	// calculate currents that depend on local Ca concentrations
	comp_ica_dyad();	
	bmyo = 1/(1 + ((trpnbar*kmtrpn)/(pow((kmtrpn+cacyt[i]),2))) + ((cmdnbar*kmcmdn)/(pow((kmcmdn+cacyt[i]),2))));

	// local myoplasmic Ca concentration
	cacyt[i] = cacyt[i] + (-bmyo*(((icab_dyad+ipca_dyad-2*inaca_dyad)*(acap/(Ntotdyads*zca*frdy*vcyt_dyad ))) + (iup-ileak)*(vnsr_dyad/vcyt_dyad) - jcadiff_ss_cyt*(vss_dyad/vcyt_dyad) ))*dt;
}

void ca_conc_diff_update()
{
	for (int i = 0; i<Ndyads; i++) //xfer ca from reaction to diffusion space
	{
		int ix,iy,iz;
	    	iz= i/(Nx*Ny);    				
	     	iy = (i - iz*(Nx*Ny))/Nx; 			
	     	ix = i -iz*(Nx*Ny) - iy*Nx;		

		cacytdiff[(Nvox*iz)*Nx_diff*Ny_diff + (Nvox*iy)*Nx_diff + (Nvox*ix)] = cacyt[i];
  		cansrdiff[(Nvox*iz)*Nx_diff*Ny_diff + (Nvox*iy)*Nx_diff + (Nvox*ix)] = cansr[i];
		cassdiff[ (Nvox*iz)*Nx_diff*Ny_diff + (Nvox*iy)*Nx_diff + (Nvox*ix)] = cass[i];
	}


	for (int i = 0; i<Nx_diff*Ny_diff*Nz_diff; i++)
	{
	     // temp variables to compute ca in diffusion space
	     cacyttemp[i] = cacytdiff[i]; 
	     cansrtemp[i] = cansrdiff[i]; 
	     casstemp[i] = cassdiff[i];
	
		int iz,iy,ix;

	     iz= i/(Nx_diff*Ny_diff);    				//3D mapping
	     iy = (i - iz*(Nx_diff*Ny_diff))/Nx_diff; 			//3D mapping
	     ix = i -iz*(Nx_diff*Ny_diff) - iy*Nx_diff;		//3D mapping

	     bsubspace = 1/(1 + ((bsrbar*kmbsr)/(pow((kmbsr+cassdiff[i]),2))) + ((bslbar*kmbsl)/(pow((kmbsl+cassdiff[i]),2))));
	     bmyo = 1/(1 + ((trpnbar*kmtrpn)/(pow((kmtrpn+cacytdiff[i]),2))) + ((cmdnbar*kmcmdn)/(pow((kmcmdn+cacytdiff[i]),2))));

		if (ix>0 && ix<Nx_diff-1)
		{ 
			cacyttemp[i] = cacyttemp[i] + Dcyt*bmyo*(cacytdiff[i-1] + cacytdiff[i+1] -2*cacytdiff[i])*dt/(dx*dx);
			casstemp[i] = casstemp[i] + Dcyt*bsubspace*(cassdiff[i-1] + cassdiff[i+1] - 2*cassdiff[i])*dt/(dx*dx);
			cansrtemp[i] = cansrtemp[i] + Dsr*(cansrdiff[i-1] + cansrdiff[i+1] -2*cansrdiff[i])*dt/(dx*dx);
		}

		if (iy>0 && iy<Ny_diff-1)
		{ 
			cacyttemp[i] = cacyttemp[i] + Dcyt*bmyo*(cacytdiff[i-Nx_diff] + cacytdiff[i+Nx_diff] -2*cacytdiff[i])*dt/(dy*dy);
			casstemp[i] = casstemp[i] + Dcyt*bsubspace*(cassdiff[i-Nx_diff] + cassdiff[i+Nx_diff] - 2*cassdiff[i])*dt/(dy*dy);
			cansrtemp[i] = cansrtemp[i] + Dsr*(cansrdiff[i-Nx_diff] + cansrdiff[i+Nx_diff] -2*cansrdiff[i])*dt/(dy*dy);
		}
	
		if (iz>0 && iz<Nz_diff-1)
		{ 
			cacyttemp[i] = cacyttemp[i] + 2*Dcyt*bmyo*(cacytdiff[i-Nx_diff*Ny_diff] + cacytdiff[i+Nx_diff*Ny_diff] -2*cacytdiff[i])*dt/(dz*dz);
			casstemp[i] = casstemp[i] + 2*Dcyt*bsubspace*(cassdiff[i-Nx_diff*Ny_diff] + cassdiff[i+Nx_diff*Ny_diff] - 2*cassdiff[i])*dt/(dz*dz);
			cansrtemp[i] = cansrtemp[i] + Dsr*(cansrdiff[i-Nx_diff*Ny_diff] + cansrdiff[i+Nx_diff*Ny_diff] -2*cansrdiff[i])*dt/(dz*dz);
		}

		if (ix==0 && ix!=(Nx_diff-1))
		{
			cacyttemp[i] = cacyttemp[i] + Dcyt*bmyo*( 2*cacytdiff[i+1] -2*cacytdiff[i])*dt/(dx*dx);
			casstemp[i] = casstemp[i] + Dcyt*bsubspace*( 2*cassdiff[i+1] - 2*cassdiff[i])*dt/(dx*dx);
			cansrtemp[i] = cansrtemp[i] + Dsr*( 2*cansrdiff[i+1] -2*cansrdiff[i])*dt/(dx*dx);
		}
		if (ix == Nx_diff-1 && ix != 0)
		{
			cacyttemp[i] = cacyttemp[i] + Dcyt*bmyo*( 2*cacytdiff[i-1] -2*cacytdiff[i])*dt/(dx*dx);
			casstemp[i] = casstemp[i] + Dcyt*bsubspace*( 2*cassdiff[i-1] - 2*cassdiff[i])*dt/(dx*dx);
			cansrtemp[i] = cansrtemp[i] + Dsr*( 2*cansrdiff[i-1] -2*cansrdiff[i])*dt/(dx*dx);
		}
		if (iy == 0 && iy != (Ny_diff-1))
		{
			cacyttemp[i] = cacyttemp[i] + Dcyt*bmyo*( 2*cacytdiff[i+Nx_diff] -2*cacytdiff[i])*dt/(dy*dy);
			casstemp[i] = casstemp[i] + Dcyt*bsubspace*( 2*cassdiff[i+Nx_diff] - 2*cassdiff[i])*dt/(dy*dy);
			cansrtemp[i] = cansrtemp[i] + Dsr*( 2*cansrdiff[i+Nx_diff] -2*cansrdiff[i])*dt/(dy*dy);
		}
		if (iy == Ny_diff-1 && iy != 0)
		{
			cacyttemp[i] = cacyttemp[i] + Dcyt*bmyo*( 2*cacytdiff[i-Nx_diff] -2*cacytdiff[i])*dt/(dy*dy);
			casstemp[i] = casstemp[i] + Dcyt*bsubspace*( 2*cassdiff[i-Nx_diff] - 2*cassdiff[i])*dt/(dy*dy);
			cansrtemp[i] = cansrtemp[i] + Dsr*( 2*cansrdiff[i-Nx_diff] -2*cansrdiff[i])*dt/(dy*dy);
		}

		if (iz == 0 && iz != (Nz_diff-1))
		{
			cacyttemp[i] = cacyttemp[i] + 2*Dcyt*bmyo*( 2*cacytdiff[i+Nx_diff*Ny_diff] -2*cacytdiff[i])*dt/(dz*dz);
			casstemp[i] = casstemp[i] + 2*Dcyt*bsubspace*( 2*cassdiff[i+Nx_diff*Ny_diff] - 2*cassdiff[i])*dt/(dz*dz);
			cansrtemp[i] = cansrtemp[i] + Dsr*( 2*cansrdiff[i+Nx_diff*Ny_diff] -2*cansrdiff[i])*dt/(dz*dz);
		}
		if (iz == Nz_diff-1 && iz != 0)
		{
			cacyttemp[i] = cacyttemp[i] + 2*Dcyt*bmyo*( 2*cacytdiff[i-Nx_diff*Ny_diff] -2*cacytdiff[i])*dt/(dz*dz);
			casstemp[i] = casstemp[i] + 2*Dcyt*bsubspace*( 2*cassdiff[i-Nx_diff*Ny_diff] - 2*cassdiff[i])*dt/(dz*dz);
			cansrtemp[i] = cansrtemp[i] + Dsr*( 2*cansrdiff[i-Nx_diff*Ny_diff] -2*cansrdiff[i])*dt/(dz*dz);
		}
	
	}

	for (int i = 0; i<Nx_diff*Ny_diff*Nz_diff; i++)
	{
		cacytdiff[i] = cacyttemp[i]; cassdiff[i] = casstemp[i]; cansrdiff[i] = cansrtemp[i];
	}


	for (int i = 0; i<Ndyads; i++) //xfer ca from diffusion to reaction space
	{
		int ix,iy,iz;
	    	iz= i/(Nx*Ny);    				
	     	iy = (i - iz*(Nx*Ny))/Nx; 			
	     	ix = i -iz*(Nx*Ny) - iy*Nx;		

		cacyt[i] = cacytdiff[(Nvox*iz)*Nx_diff*Ny_diff + (Nvox*iy)*Nx_diff + (Nvox*ix)];
  		cansr[i] = cansrdiff[(Nvox*iz)*Nx_diff*Ny_diff + (Nvox*iy)*Nx_diff + (Nvox*ix)];
		cass[i] = cassdiff[ (Nvox*iz)*Nx_diff*Ny_diff + (Nvox*iy)*Nx_diff + (Nvox*ix)];
	}
}


void prttofile ()
{
	//Whole-cell variables
	f0out<<t<<"	"<<v<<"	"<<cadyad<<"	"<<casub<<"	"<<cai<<"	"<<nsr<<"	"<<jsr<<endl;
	f1out<<t<<"	"<<ina<<"	"<<ICaL_print<<"	"<<Ito<<"	"<<IKs<<"	"<<IKr<<"	"<<IK1<<"	"<<INaK<<"	"<<inaca+inacass<<"	"<<jrel<<endl;
	f2out<<t<<"	"<<cads[Ndyads/2]<<"	"<<cass[Ndyads/2]<<"	"<<cacyt[Ndyads/2]<<"	"<<cansr[Ndyads/2]<<"	"<<cajsr[Ndyads/2]<<endl;

		if (icounter%50==0)
		{
			for (int i = 0; i<Nx; i++)
			{
				f3out<<cacyt[i*Nz*Ny + Nz*Ny/2]<<"	";
				f4out<<cajsr[i*Nz*Ny + Nz*Ny/2]<<"	";
			}
			f3out<<endl;
			f4out<<endl;
		}

}

void prtfinal()
{
f5out<<"v="<<v<<";"<<endl;
f5out<<"nai="<<nai<<";"<<endl;
f5out<<"ki="<<ki<<";"<<endl;
f5out<<"cai="<<cai<<";"<<endl;
f5out<<"casub="<<casub<<";"<<endl;
f5out<<"nsr="<<nsr<<";"<<endl;
f5out<<"jsr="<<jsr<<";"<<endl;
f5out<<"m="<<m<<";"<<endl;
f5out<<"hf="<<hf<<";"<<endl;
f5out<<"hs="<<hs<<";"<<endl;
f5out<<"j="<<j<<";"<<endl;
f5out<<"hsp="<<hsp<<";"<<endl;
f5out<<"jp="<<jp<<";"<<endl;
f5out<<"ml="<<ml<<";"<<endl;
f5out<<"hl="<<hl<<";"<<endl;
f5out<<"hlp="<<hlp<<";"<<endl;
f5out<<"a="<<a<<";"<<endl;
f5out<<"iF="<<iF<<";"<<endl;
f5out<<"iS="<<iS<<";"<<endl;
f5out<<"ap="<<ap<<";"<<endl;
f5out<<"iFp="<<iFp<<";"<<endl;
f5out<<"iSp="<<iSp<<";"<<endl;
f5out<<"d="<<d<<";"<<endl;
f5out<<"ff="<<ff<<";"<<endl;
f5out<<"fs="<<fs<<";"<<endl;
f5out<<"fcaf="<<fcaf<<";"<<endl;
f5out<<"fcas="<<fcas<<";"<<endl;
f5out<<"jca="<<jca<<";"<<endl;
f5out<<"ffp="<<ffp<<";"<<endl;
f5out<<"fcafp="<<fcafp<<";"<<endl;
f5out<<"xrf="<<xrf<<";"<<endl;
f5out<<"xrs="<<xrs<<";"<<endl;
f5out<<"xs1="<<xs1<<";"<<endl;
f5out<<"xs2="<<xs2<<";"<<endl;
f5out<<"xk1="<<xk1<<";"<<endl;

}

void comp_ica_dyad()
{
double kna1=15.0;
double kna2=5.0;
double kna3=88.12;
double kasymm=12.5;
double wna=6.0e4;
double wca=6.0e4;
double wnaca=5.0e3;
double kcaon=1.5e6;
double kcaoff=5.0e3;
double qna=0.5224;
double qca=0.1670;
double hca=exp((qca*v*frdy)/(R*temp));
double hna=exp((qna*v*frdy)/(R*temp));
double h1=1+nai/kna3*(1+hna);
double h2=(nai*hna)/(kna3*h1);
double h3=1.0/h1;
double h4=1.0+nai/kna1*(1+nai/kna2);
double h5=nai*nai/(h4*kna1*kna2);
double h6=1.0/h4;
double h7=1.0+nao/kna3*(1.0+1.0/hna);
double h8=nao/(kna3*hna*h7);
double h9=1.0/h7;
double h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
double h11=nao*nao/(h10*kna1*kna2);
double h12=1.0/h10;
double k1=h12*cao*kcaon;
double k2=kcaoff;
double k3p=h9*wca;
double k3pp=h8*wnaca;
double k3=k3p+k3pp;
double k4p=h3*wca/hca;
double k4pp=h2*wnaca;
double k4=k4p+k4pp;
double k5=kcaoff;
double k6=h6*cacyt[i]*kcaon;
double k7=h5*h2*wna;
double k8=h8*h11*wna;
double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
double E1=x1/(x1+x2+x3+x4);
double E2=x2/(x1+x2+x3+x4);
double E3=x3/(x1+x2+x3+x4);
double E4=x4/(x1+x2+x3+x4);
double KmCaAct=150.0e-6;
double allo=1.0/(1.0+pow(KmCaAct/cacyt[i],2.0));
double zna=1.0;
double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
double JncxCa=E2*k2-E1*k1;
double Gncx=0.0008;


inaca_dyad=0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa);

h1=1+nai/kna3*(1+hna);
h2=(nai*hna)/(kna3*h1);
h3=1.0/h1;
h4=1.0+nai/kna1*(1+nai/kna2);
h5=nai*nai/(h4*kna1*kna2);
h6=1.0/h4;
h7=1.0+nao/kna3*(1.0+1.0/hna);
h8=nao/(kna3*hna*h7);
h9=1.0/h7;
h10=kasymm+1.0+nao/kna1*(1+nao/kna2);
h11=nao*nao/(h10*kna1*kna2);
h12=1.0/h10;
k1=h12*cao*kcaon;
k2=kcaoff;
k3p=h9*wca;
k3pp=h8*wnaca;
k3=k3p+k3pp;
k4p=h3*wca/hca;
k4pp=h2*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6*cass[i]*kcaon;
k7=h5*h2*wna;
k8=h8*h11*wna;
x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
KmCaAct=150.0e-6;
allo=1.0/(1.0+pow(KmCaAct/cass[i],2.0));
JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
JncxCa=E2*k2-E1*k1;

inacass_dyad=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);

double PCab=2.5e-8;
double vffrt = (v*frdy*frdy)/(R*temp);
double vfrt = v*frdy/(R*temp);
icab_dyad=PCab*4.0*vffrt*(cacyt[i]*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);

double GpCa=0.0005;
ipca_dyad=GpCa*cacyt[i]/(0.0005+cacyt[i]);

inaca += (inaca_dyad)/Ndyads;
inacass += (inacass_dyad)/Ndyads;
ipca += (ipca_dyad)/Ndyads;
icab += (icab_dyad)/Ndyads;
}
 
