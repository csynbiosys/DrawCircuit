parsNamesTranscription={'n_cI','KM_cI','a_tr_cI','a0_tr_cI'};
parsValsTranslation=[2 40 29.9700 0.0300];
interactions={'cI_protein'};
transcriptionFunction='a0_tr_cI+a_tr_cI*pow(KM_cI,n_cI)/(pow(KM_cI,n_cI)+pow(cI_protein,n_cI))';