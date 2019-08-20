parsNamesTranscription={'n_LacI','KM_LacI','a_tr_LacI','a0_tr_LacI'};
parsValsTranslation=[2 40 29.9700 0.0300];
interactions={'LacI_protein'};
transcriptionFunction='a0_tr_LacI+a_tr_LacI*pow(KM_LacI,n_LacI)/(pow(KM_LacI,n_LacI)+pow(LacI_protein,n_LacI))';