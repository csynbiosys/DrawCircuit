#include <amigoRHS.h>
#include <math.h>
#include <amigoJAC.h>
#include <amigoSensRHS.h>

double FunctionForTranscriptionOfCI(double param_0, double modif_0, double param_1, double param_2, double volume_0, double param_3);

#define LacIprot  Ith(y,0) 	//concentration of metabolite 'LacI protein': reactions
#define TetRprot  Ith(y,1) 	//concentration of metabolite 'TetR protein': reactions
#define cIprot  Ith(y,2) 	//concentration of metabolite 'cI protein': reactions
#define LacImRNA  Ith(y,3) 	//concentration of metabolite 'LacI mRNA': reactions
#define TetRmRNA  Ith(y,4) 	//concentration of metabolite 'TetR mRNA': reactions
#define cImRNA  Ith(y,5) 	//concentration of metabolite 'cI mRNA': reactions

int amigoRHS(realtype t, N_Vector y, N_Vector ydot, void *data){
    AMIGO_model* amigo_model=(AMIGO_model*)data;
    
    
    double translationEfficiency = 20;	//global quantity 'translation efficiency':fixed
    double n = 2;	//global quantity 'n':fixed
    double KM = 40;	//global quantity 'KM':fixed
    double mRNAhalfLife = 2;	//global quantity 'mRNA half life':fixed
    double protHalfLife = 10;	//global quantity 'protein half life':fixed
    double tps_active = 0.5;	//global quantity 'tps_active':fixed
    double tps_repr = 0.0005;	//global quantity 'tps_repr':fixed
    double a_tr = (tps_active-tps_repr)*60.00000000000000000;	//model entity 'a_tr':assignment
    double a0_tr = tps_repr*60.00000000000000000;	//model entity 'a0_tr':assignment
    double beta = mRNAhalfLife/protHalfLife;	//model entityL 'beta':assignment
    double alpha0 = a0_tr*translationEfficiency*protHalfLife/(log(2.00000000000000000)*KM);	//model entity 'alpha0':assignment
    double alpha = a_tr*translationEfficiency*protHalfLife/(log(2.00000000000000000)*KM);	//model entity 'alpha':assignment
    double averagemRNAlifetime = mRNAhalfLife/log(2.00000000000000000);	//model entity 'average mRNA life time':assignment
    double kd_mRNA = log(2.00000000000000000)/mRNAhalfLife;	//model entity 'kd_mRNA':assignment
    double kd_prot = log(2.00000000000000000)/protHalfLife;	//model entity 'kd_prot':assignment
    double k_tl = translationEfficiency/averagemRNAlifetime;	//model entity 'k_tl':assignment
    
    Ith(ydot,0)  = LacImRNA *  k_tl -  kd_prot * LacIprot;
    Ith(ydot,1)  = TetRmRNA * k_tl - kd_prot * TetRprot ;
    Ith(ydot,2)  = cImRNA * k_tl -kd_prot * cIprot;
    Ith(ydot,3)  = -kd_mRNA * LacImRNA+a0_tr+a_tr*pow(KM,n)/(pow(KM,n)+pow(cIprot,n));
    Ith(ydot,4)  = -kd_mRNA * TetRmRNA+a0_tr+a_tr*pow(KM,n)/(pow(KM,n)+pow(LacIprot,n));
    Ith(ydot,5)  = -kd_mRNA * cImRNA+a0_tr+a_tr*pow(KM,n)/(pow(KM,n)+pow(TetRprot,n));
    
    return(0);
}

double FunctionForTranscriptionOfCI(double param_0, double modif_0, double param_1, double param_2, double volume_0, double param_3) 	//Function for transcription of CI
{return  param_1+param_2*pow(param_0,param_3)/(pow(param_0,param_3)+pow(modif_0,param_3));}


/* Jacobian of the system (dfdx)*/
int amigoJAC(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    AMIGO_model* amigo_model=(AMIGO_model*)user_data;
    
    return(0);
}

/* R.H.S of the sensitivity dsi/dt = (df/dx)*si + df/dp_i */
int amigoSensRHS(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *data, N_Vector tmp1, N_Vector tmp2){
    AMIGO_model* amigo_model=(AMIGO_model*)data;
    
    return(0);
}

void amigoRHS_get_OBS(void* data){
}

void amigoRHS_get_sens_OBS(void* data){
}


void amigo_Y_at_tcon(void* data,realtype t, N_Vector y){
    AMIGO_model* amigo_model=(AMIGO_model*)data;
}


