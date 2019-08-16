parsNamesTranscription={'n_TetR','KM_TetR','a_tr_TetR','a0_tr_TetR'};
parsValsTranslation=[2 40 29.9700 0.0300];
interactions={'TetR_prot'};
transcriptionFunction='a0_tr_TetR+a_tr_TetR*pow(KM_TetR,n_TetR)/(pow(KM_TetR,n_TetR)+pow(TetR_prot,n_TetR))';