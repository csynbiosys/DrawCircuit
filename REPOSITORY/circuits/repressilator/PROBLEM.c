#include <amigoRHS.h>

#include <math.h>

#include <amigoJAC.h>

#include <amigoSensRHS.h>


	/* *** Definition of the states *** */

#define	LacI_CDS     Ith(y,0)
#define	TetR_CDS     Ith(y,1)
#define	cI_CDS       Ith(y,2)
#define	LacI_protein Ith(y,3)
#define	TetR_protein Ith(y,4)
#define	cI_protein   Ith(y,5)
#define iexp amigo_model->exp_num

	/* *** Definition of the sates derivative *** */

#define	dLacI_CDS     Ith(ydot,0)
#define	dTetR_CDS     Ith(ydot,1)
#define	dcI_CDS       Ith(ydot,2)
#define	dLacI_protein Ith(ydot,3)
#define	dTetR_protein Ith(ydot,4)
#define	dcI_protein   Ith(ydot,5)

	/* *** Definition of the parameters *** */

#define	n_LacI                           (*amigo_model).pars[0]
#define	KM_LacI                          (*amigo_model).pars[1]
#define	a_tr_LacI                        (*amigo_model).pars[2]
#define	a0_tr_LacI                       (*amigo_model).pars[3]
#define	n_TetR                           (*amigo_model).pars[4]
#define	KM_TetR                          (*amigo_model).pars[5]
#define	a_tr_TetR                        (*amigo_model).pars[6]
#define	a0_tr_TetR                       (*amigo_model).pars[7]
#define	n_cI                             (*amigo_model).pars[8]
#define	KM_cI                            (*amigo_model).pars[9]
#define	a_tr_cI                          (*amigo_model).pars[10]
#define	a0_tr_cI                         (*amigo_model).pars[11]
#define	translationEffiency_LacIRep_PROM (*amigo_model).pars[12]
#define	translationEffiency_TetRRep_PROM (*amigo_model).pars[13]
#define	translationEffiency_cIRep_PROM   (*amigo_model).pars[14]
#define	k_tl_LacI                        (*amigo_model).pars[15]
#define	kdeg_mRNA_LacI                   (*amigo_model).pars[16]
#define	k_deg_LacI_protein               (*amigo_model).pars[17]
#define	k_tl_TetR                        (*amigo_model).pars[18]
#define	kdeg_mRNA_TetR                   (*amigo_model).pars[19]
#define	k_deg_TetR_protein               (*amigo_model).pars[20]
#define	k_tl_cI                          (*amigo_model).pars[21]
#define	kdeg_mRNA_cI                     (*amigo_model).pars[22]
#define	k_deg_cI_protein                 (*amigo_model).pars[23]
#define	Y1                               (*amigo_model).pars[24]
#define	Y2                               (*amigo_model).pars[25]
#define	Y3                               (*amigo_model).pars[26]
#define	Y4                               (*amigo_model).pars[27]
#define	Y5                               (*amigo_model).pars[28]
#define	Y6                               (*amigo_model).pars[29]
#define	Y7                               (*amigo_model).pars[30]
#define	Y8                               (*amigo_model).pars[31]
#define	Y9                               (*amigo_model).pars[32]
/* Right hand side of the system (f(t,x,p))*/
int amigoRHS(realtype t, N_Vector y, N_Vector ydot, void *data){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	/* *** Definition of the algebraic variables *** */

	double	LacIRep_PROM;
	double	TetRRep_PROM;
	double	cIRep_PROM;

	/* *** Equations *** */

	LacIRep_PROM=a0_tr_LacI+a_tr_LacI*pow(KM_LacI,n_LacI)/(pow(KM_LacI,n_LacI)+pow(LacI_protein,n_LacI));
	TetRRep_PROM=a0_tr_TetR+a_tr_TetR*pow(KM_TetR,n_TetR)/(pow(KM_TetR,n_TetR)+pow(TetR_protein,n_TetR));
	cIRep_PROM=a0_tr_cI+a_tr_cI*pow(KM_cI,n_cI)/(pow(KM_cI,n_cI)+pow(cI_protein,n_cI));
	dLacI_CDS=+translationEffiency_LacIRep_PROM*Y1*LacIRep_PROM+translationEffiency_LacIRep_PROM*Y2*LacIRep_PROM+translationEffiency_LacIRep_PROM*Y3*LacIRep_PROM;
	dTetR_CDS=+translationEffiency_TetRRep_PROM*Y4*TetRRep_PROM+translationEffiency_TetRRep_PROM*Y5*TetRRep_PROM+translationEffiency_TetRRep_PROM*Y6*TetRRep_PROM;
	dcI_CDS=+translationEffiency_cIRep_PROM*Y7*cIRep_PROM+translationEffiency_cIRep_PROM*Y8*cIRep_PROM+translationEffiency_cIRep_PROM*Y9*cIRep_PROM;
	dLacI_protein=k_tl_LacI*LacI_CDS-k_deg_LacI_protein*LacI_protein;
	dTetR_protein=k_tl_TetR*TetR_CDS-k_deg_TetR_protein*TetR_protein;
	dcI_protein=k_tl_cI*cI_CDS-k_deg_cI_protein*cI_protein;

	return(0);

}


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

#define	 LacI_CDS     (amigo_model->sim_results[0][j]) 
#define	 TetR_CDS     (amigo_model->sim_results[1][j]) 
#define	 cI_CDS       (amigo_model->sim_results[2][j]) 
#define	 LacI_protein (amigo_model->sim_results[3][j]) 
#define	 TetR_protein (amigo_model->sim_results[4][j]) 
#define	 cI_protein   (amigo_model->sim_results[5][j]) 



void amigoRHS_get_OBS(void* data){

	int j;

	double t;
	AMIGO_model* amigo_model=(AMIGO_model*)data;


	 switch (amigo_model->exp_num){

		#define	 LacI_CDS_obs     amigo_model->obs_results[0][j] 
		#define	 TetR_CDS_obs     amigo_model->obs_results[1][j] 
		#define	 cI_CDS_obs       amigo_model->obs_results[2][j] 
		#define	 LacI_protein_obs amigo_model->obs_results[3][j] 
		#define	 TetR_protein_obs amigo_model->obs_results[4][j] 
		#define	 cI_protein_obs   amigo_model->obs_results[5][j] 

		 case 0:


			 for (j = 0; j < amigo_model->n_times; ++j){

				t=amigo_model->t[j];
				LacI_CDS_obs=LacI_CDS        ;
				TetR_CDS_obs=TetR_CDS        ;
				cI_CDS_obs=cI_CDS            ;
				LacI_protein_obs=LacI_protein;
				TetR_protein_obs=TetR_protein;
				cI_protein_obs=cI_protein    ;

			}

		 break;

	}

	

}

#define	 LacI_CDS     (amigo_model->sens_results[0][j][k]) 
#define	 TetR_CDS     (amigo_model->sens_results[1][j][k]) 
#define	 cI_CDS       (amigo_model->sens_results[2][j][k]) 
#define	 LacI_protein (amigo_model->sens_results[3][j][k]) 
#define	 TetR_protein (amigo_model->sens_results[4][j][k]) 
#define	 cI_protein   (amigo_model->sens_results[5][j][k]) 



void amigoRHS_get_sens_OBS(void* data){
	int j,k;

	double t;

	AMIGO_model* amigo_model=(AMIGO_model*)data;


	 switch (amigo_model->exp_num){


		 case 0:

		#define	 LacI_CDS_obs     amigo_model->sens_obs[0][j][k] 
		#define	 TetR_CDS_obs     amigo_model->sens_obs[1][j][k] 
		#define	 cI_CDS_obs       amigo_model->sens_obs[2][j][k] 
		#define	 LacI_protein_obs amigo_model->sens_obs[3][j][k] 
		#define	 TetR_protein_obs amigo_model->sens_obs[4][j][k] 
		#define	 cI_protein_obs   amigo_model->sens_obs[5][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){

				t=amigo_model->t[k];
					LacI_CDS_obs=LacI_CDS        ;
					TetR_CDS_obs=TetR_CDS        ;
					cI_CDS_obs=cI_CDS            ;
					LacI_protein_obs=LacI_protein;
					TetR_protein_obs=TetR_protein;
					cI_protein_obs=cI_protein    ;
				}
			}
		 break;
	}
}


void amigo_Y_at_tcon(void* data,realtype t, N_Vector y){
	AMIGO_model* amigo_model=(AMIGO_model*)data;


}
