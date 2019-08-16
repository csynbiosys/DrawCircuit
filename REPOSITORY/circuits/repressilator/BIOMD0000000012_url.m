#ifdef SIZE_DEFINITIONS
#define N_METABS 6
#define N_ODE_METABS 0
#define N_INDEP_METABS 6
#define N_COMPARTMENTS 1
#define N_GLOBAL_PARAMS 16
#define N_KIN_PARAMS 0
#define N_REACTIONS 12

#define N_ARRAY_SIZE_P  8	// number of parameters
#define N_ARRAY_SIZE_X  6	// number of initials
#define N_ARRAY_SIZE_Y  9	// number of assigned elements
#define N_ARRAY_SIZE_XC 6	// number of x concentration
#define N_ARRAY_SIZE_PC 0	// number of p concentration
#define N_ARRAY_SIZE_YC 0	// number of y concentration
#define N_ARRAY_SIZE_DX 6	// number of ODEs 
#define N_ARRAY_SIZE_CT 0	// number of conserved totals

#endif // SIZE_DEFINITIONS

#ifdef TIME
#define T  <set here a user name for the time variable> 
#endif // TIME

#ifdef NAME_ARRAYS
const char* p_names[] = {"cell", "translation efficiency", "n", "KM", "mRNA half life", "protein half life", "tps_active", "tps_repr",  "" };
const char* x_names[] = {"LacI protein", "TetR protein", "cI protein", "LacI mRNA", "TetR mRNA", "cI mRNA",  "" };
const char* y_names[] = {"beta", "alpha0", "alpha", "average mRNA life time", "kd_mRNA", "kd_prot", "k_tl", "a_tr", "a0_tr",  "" };
const char* xc_names[] = {"LacI protein", "TetR protein", "cI protein", "LacI mRNA", "TetR mRNA", "cI mRNA",  "" };
const char* pc_names[] = { "" };
const char* yc_names[] = { "" };
const char* dx_names[] = {"ODE LacI protein", "ODE TetR protein", "ODE cI protein", "ODE LacI mRNA", "ODE TetR mRNA", "ODE cI mRNA",  "" };
const char* ct_names[] = { "" };
#endif // NAME_ARRAYS

#ifdef INITIAL
x[0] = 0;	//metabolite 'LacI protein': reactions
x[1] = 0;	//metabolite 'TetR protein': reactions
x[2] = 0;	//metabolite 'cI protein': reactions
x[3] = 0;	//metabolite 'LacI mRNA': reactions
x[4] = 20;	//metabolite 'TetR mRNA': reactions
x[5] = 0;	//metabolite 'cI mRNA': reactions
#endif /* INITIAL */

#ifdef FIXED
p[0] = 1;	//compartment 'cell':fixed
p[1] = 20;	//global quantity 'translation efficiency':fixed
p[2] = 2;	//global quantity 'n':fixed
p[3] = 40;	//global quantity 'KM':fixed
p[4] = 2;	//global quantity 'mRNA half life':fixed
p[5] = 10;	//global quantity 'protein half life':fixed
p[6] = 0.5;	//global quantity 'tps_active':fixed
p[7] = 0.0005;	//global quantity 'tps_repr':fixed
#endif /* FIXED */

#ifdef ASSIGNMENT
y[0] = p[4]/p[5];	//model entity 'beta':assignment
y[1] = y[8]*p[1]*p[5]/(log(2.00000000000000000)*p[3]);	//model entity 'alpha0':assignment
y[2] = y[7]*p[1]*p[5]/(log(2.00000000000000000)*p[3]);	//model entity 'alpha':assignment
y[3] = p[4]/log(2.00000000000000000);	//model entity 'average mRNA life time':assignment
y[4] = log(2.00000000000000000)/p[4];	//model entity 'kd_mRNA':assignment
y[5] = log(2.00000000000000000)/p[5];	//model entity 'kd_prot':assignment
y[6] = p[1]/y[3];	//model entity 'k_tl':assignment
y[7] = (p[6]-p[7])*60.00000000000000000;	//model entity 'a_tr':assignment
y[8] = p[7]*60.00000000000000000;	//model entity 'a0_tr':assignment
x_c[0] = x[0]/p[0];	//concentration of metabolite 'LacI protein': reactions
x_c[1] = x[1]/p[0];	//concentration of metabolite 'TetR protein': reactions
x_c[2] = x[2]/p[0];	//concentration of metabolite 'cI protein': reactions
x_c[3] = x[3]/p[0];	//concentration of metabolite 'LacI mRNA': reactions
x_c[4] = x[4]/p[0];	//concentration of metabolite 'TetR mRNA': reactions
x_c[5] = x[5]/p[0];	//concentration of metabolite 'cI mRNA': reactions
#endif /* ASSIGNMENT */

#ifdef FUNCTIONS_HEADERS
double FunctionForTranslationOfLacI(double modif_0, double param_0); 
double FunctionForTranslationOfTetR(double modif_0, double param_0); 
double FunctionForTranslationOfCI(double modif_0, double param_0); 
double FunctionForTranscriptionOfLacI(double param_0, double modif_0, double param_1, double param_2, double volume_0, double param_3); 
double FunctionForTranscriptionOfTetR(double param_0, double modif_0, double param_1, double param_2, double volume_0, double param_3); 
double FunctionForTranscriptionOfCI(double param_0, double modif_0, double param_1, double param_2, double volume_0, double param_3); 
#endif /* FUNCTIONS_HEADERS */

#ifdef FUNCTIONS
double FunctionForTranslationOfLacI(double modif_0, double param_0) 	//Function for translation of LacI
{return  param_0*modif_0;} 
double FunctionForTranslationOfTetR(double modif_0, double param_0) 	//Function for translation of TetR
{return  param_0*modif_0;} 
double FunctionForTranslationOfCI(double modif_0, double param_0) 	//Function for translation of CI
{return  param_0*modif_0;} 
double FunctionForTranscriptionOfLacI(double param_0, double modif_0, double param_1, double param_2, double volume_0, double param_3) 	//Function for transcription of LacI
{return  (param_1+param_2*pow(param_0,param_3)/(pow(param_0,param_3)+pow((modif_0*volume_0),param_3)))/volume_0;} 
double FunctionForTranscriptionOfTetR(double param_0, double modif_0, double param_1, double param_2, double volume_0, double param_3) 	//Function for transcription of TetR
{return  (param_1+param_2*pow(param_0,param_3)/(pow(param_0,param_3)+pow((modif_0*volume_0),param_3)))/volume_0;} 
double FunctionForTranscriptionOfCI(double param_0, double modif_0, double param_1, double param_2, double volume_0, double param_3) 	//Function for transcription of CI
{return  (param_1+param_2*pow(param_0,param_3)/(pow(param_0,param_3)+pow((modif_0*volume_0),param_3)))/volume_0;} 
#endif /* FUNCTIONS */

#ifdef ODEs
dx[0] = FunctionForTranslationOfLacI(x_c[3], y[6])*p[0]-(y[5] * x_c[0]) *p[0];
dx[1] = FunctionForTranslationOfTetR(x_c[4], y[6])*p[0]-(y[5] * x_c[1]) *p[0];
dx[2] = FunctionForTranslationOfCI(x_c[5], y[6])*p[0]-(y[5] * x_c[2]) *p[0];
dx[3] = -(y[4] * x_c[3]) *p[0]+FunctionForTranscriptionOfLacI(p[3], x_c[2], y[8], y[7], p[0], p[2])*p[0];
dx[4] = -(y[4] * x_c[4]) *p[0]+FunctionForTranscriptionOfTetR(p[3], x_c[0], y[8], y[7], p[0], p[2])*p[0];
dx[5] = -(y[4] * x_c[5]) *p[0]+FunctionForTranscriptionOfCI(p[3], x_c[1], y[8], y[7], p[0], p[2])*p[0];
#endif /* ODEs */
