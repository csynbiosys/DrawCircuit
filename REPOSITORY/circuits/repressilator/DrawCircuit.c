#include <amigoRHS.h>

#include <math.h>

#include <amigoJAC.h>

#include <amigoSensRHS.h>


	/* *** Definition of the states *** */

#define	LACI Ith(y,0)
#define	TETR Ith(y,1)
#define	IPTG Ith(y,2)
#define	ATC  Ith(y,3)
#define iexp amigo_model->exp_num

	/* *** Definition of the sates derivative *** */

#define	dLACI Ith(ydot,0)
#define	dTETR Ith(ydot,1)
#define	dIPTG Ith(ydot,2)
#define	dATC  Ith(ydot,3)

	/* *** Definition of the parameters *** */

#define	k_promoter1_TETR (*amigo_model).pars[0]
#define	n_promoter1_TETR (*amigo_model).pars[1]
#define	k_promoter1_ATC  (*amigo_model).pars[2]
#define	n_promoter1_ATC  (*amigo_model).pars[3]
#define	k_promoter2_LACI (*amigo_model).pars[4]
#define	n_promoter2_LACI (*amigo_model).pars[5]
#define	k_promoter2_IPTG (*amigo_model).pars[6]
#define	n_promoter2_IPTG (*amigo_model).pars[7]
#define	IPTG_kdiff       (*amigo_model).pars[8]
#define	ATC_kdiff        (*amigo_model).pars[9]
#define	w_promoter1_LACI (*amigo_model).pars[10]
#define	w_promoter2_LACI (*amigo_model).pars[11]
#define	tau_LACI         (*amigo_model).pars[12]
#define	w_promoter1_TETR (*amigo_model).pars[13]
#define	w_promoter2_TETR (*amigo_model).pars[14]
#define	tau_TETR         (*amigo_model).pars[15]
#define IPTGi	((*amigo_model).controls_v[0][(*amigo_model).index_t_stim]+(t-(*amigo_model).tlast)*(*amigo_model).slope[0][(*amigo_model).index_t_stim])
#define ATCi 	((*amigo_model).controls_v[1][(*amigo_model).index_t_stim]+(t-(*amigo_model).tlast)*(*amigo_model).slope[1][(*amigo_model).index_t_stim])
/* Right hand side of the system (f(t,x,p))*/
int amigoRHS(realtype t, N_Vector y, N_Vector ydot, void *data){
	AMIGO_model* amigo_model=(AMIGO_model*)data;

	/* *** Definition of the algebraic variables *** */

	double	promoter1;
	double	promoter2;

	/* *** Equations *** */

	promoter1=(1-(pow(ATC,n_promoter1_ATC)*(pow(k_promoter1_ATC,n_promoter1_ATC)+1))/(pow(k_promoter1_ATC,n_promoter1_ATC)+pow(ATC,n_promoter1_ATC)))*(1-(pow(TETR,n_promoter1_TETR)*(pow(k_promoter1_TETR,n_promoter1_TETR)+1))/(pow(k_promoter1_TETR,n_promoter1_TETR)+pow(TETR,n_promoter1_TETR)));
	promoter2=(1-(pow(LACI,n_promoter2_LACI)*(pow(k_promoter2_LACI,n_promoter2_LACI)+1))/(pow(k_promoter2_LACI,n_promoter2_LACI)+pow(LACI,n_promoter2_LACI)))*(1-(pow(IPTG,n_promoter2_IPTG)*(pow(k_promoter2_IPTG,n_promoter2_IPTG)+1))/(pow(k_promoter2_IPTG,n_promoter2_IPTG)+pow(IPTG,n_promoter2_IPTG)));
	dLACI=(w_promoter1_LACI*promoter1+w_promoter2_LACI*promoter2-LACI)*tau_LACI;
	dTETR=(w_promoter1_TETR*promoter1+w_promoter2_TETR*promoter2-TETR)*tau_TETR;
	dIPTG=(IPTGi-IPTG)*IPTG_kdiff;
	dATC=(ATCi-ATC)*ATC_kdiff;

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

#define	 LACI (amigo_model->sim_results[0][j]) 
#define	 TETR (amigo_model->sim_results[1][j]) 
#define	 IPTG (amigo_model->sim_results[2][j]) 
#define	 ATC  (amigo_model->sim_results[3][j]) 



void amigoRHS_get_OBS(void* data){

	int j;

	double t;
	AMIGO_model* amigo_model=(AMIGO_model*)data;


	 switch (amigo_model->exp_num){

		#define	 LACI_o amigo_model->obs_results[0][j] 
		#define	 TETR_o amigo_model->obs_results[1][j] 

		 case 0:


			 for (j = 0; j < amigo_model->n_times; ++j){

				t=amigo_model->t[j];
				LACI_o=LACI;
				TETR_o=TETR;

			}

		 break;
		#define	 LACI_o amigo_model->obs_results[0][j] 
		#define	 TETR_o amigo_model->obs_results[1][j] 

		 case 1:


			 for (j = 0; j < amigo_model->n_times; ++j){

				t=amigo_model->t[j];
				LACI_o=LACI;
				TETR_o=TETR;

			}

		 break;
		#define	 LACI_o amigo_model->obs_results[0][j] 
		#define	 TETR_o amigo_model->obs_results[1][j] 

		 case 2:


			 for (j = 0; j < amigo_model->n_times; ++j){

				t=amigo_model->t[j];
				LACI_o=LACI;
				TETR_o=TETR;

			}

		 break;
		#define	 LACI_o amigo_model->obs_results[0][j] 
		#define	 TETR_o amigo_model->obs_results[1][j] 

		 case 3:


			 for (j = 0; j < amigo_model->n_times; ++j){

				t=amigo_model->t[j];
				LACI_o=LACI;
				TETR_o=TETR;

			}

		 break;
		#define	 LACI_o amigo_model->obs_results[0][j] 
		#define	 TETR_o amigo_model->obs_results[1][j] 

		 case 4:


			 for (j = 0; j < amigo_model->n_times; ++j){

				t=amigo_model->t[j];
				LACI_o=LACI;
				TETR_o=TETR;

			}

		 break;
		#define	 LACI_o amigo_model->obs_results[0][j] 
		#define	 TETR_o amigo_model->obs_results[1][j] 

		 case 5:


			 for (j = 0; j < amigo_model->n_times; ++j){

				t=amigo_model->t[j];
				LACI_o=LACI;
				TETR_o=TETR;

			}

		 break;
		#define	 LACI_o amigo_model->obs_results[0][j] 
		#define	 TETR_o amigo_model->obs_results[1][j] 

		 case 6:


			 for (j = 0; j < amigo_model->n_times; ++j){

				t=amigo_model->t[j];
				LACI_o=LACI;
				TETR_o=TETR;

			}

		 break;

	}

	

}

#define	 LACI (amigo_model->sens_results[0][j][k]) 
#define	 TETR (amigo_model->sens_results[1][j][k]) 
#define	 IPTG (amigo_model->sens_results[2][j][k]) 
#define	 ATC  (amigo_model->sens_results[3][j][k]) 



void amigoRHS_get_sens_OBS(void* data){
	int j,k;

	double t;

	AMIGO_model* amigo_model=(AMIGO_model*)data;


	 switch (amigo_model->exp_num){


		 case 0:

		#define	 LACI_o amigo_model->sens_obs[0][j][k] 
		#define	 TETR_o amigo_model->sens_obs[1][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){

				t=amigo_model->t[k];
					LACI_o=LACI;
					TETR_o=TETR;
				}
			}
		 break;

		 case 1:

		#define	 LACI_o amigo_model->sens_obs[0][j][k] 
		#define	 TETR_o amigo_model->sens_obs[1][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){

				t=amigo_model->t[k];
					LACI_o=LACI;
					TETR_o=TETR;
				}
			}
		 break;

		 case 2:

		#define	 LACI_o amigo_model->sens_obs[0][j][k] 
		#define	 TETR_o amigo_model->sens_obs[1][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){

				t=amigo_model->t[k];
					LACI_o=LACI;
					TETR_o=TETR;
				}
			}
		 break;

		 case 3:

		#define	 LACI_o amigo_model->sens_obs[0][j][k] 
		#define	 TETR_o amigo_model->sens_obs[1][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){

				t=amigo_model->t[k];
					LACI_o=LACI;
					TETR_o=TETR;
				}
			}
		 break;

		 case 4:

		#define	 LACI_o amigo_model->sens_obs[0][j][k] 
		#define	 TETR_o amigo_model->sens_obs[1][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){

				t=amigo_model->t[k];
					LACI_o=LACI;
					TETR_o=TETR;
				}
			}
		 break;

		 case 5:

		#define	 LACI_o amigo_model->sens_obs[0][j][k] 
		#define	 TETR_o amigo_model->sens_obs[1][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){

				t=amigo_model->t[k];
					LACI_o=LACI;
					TETR_o=TETR;
				}
			}
		 break;

		 case 6:

		#define	 LACI_o amigo_model->sens_obs[0][j][k] 
		#define	 TETR_o amigo_model->sens_obs[1][j][k] 

			 for (j = 0; j < amigo_model->n_total_x; ++j){
				 for (k = 0; k < amigo_model->n_times; ++k){

				t=amigo_model->t[k];
					LACI_o=LACI;
					TETR_o=TETR;
				}
			}
		 break;
	}
}


void amigo_Y_at_tcon(void* data,realtype t, N_Vector y){
	AMIGO_model* amigo_model=(AMIGO_model*)data;


}
