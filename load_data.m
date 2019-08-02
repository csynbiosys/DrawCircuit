function [exps] = load_data(file,normalize)

DATA = importdata(file);

EXPERIMENT_COL=find(startsWith(DATA.textdata,'Experiment'));
TIME_COL=find(startsWith(DATA.textdata,'Time'));
STIMULI_COLS=find(startsWith(DATA.textdata,'TR_'));
READOUT_COLS=find(startsWith(DATA.textdata,'READOUT_'));
STD_COLS=find(startsWith(DATA.textdata,'STD_'));

EXPERIMENTS=unique(sort(DATA.data(:,EXPERIMENT_COL)));
N_EXPERIMENTS=length(EXPERIMENTS);
exps.n_exp=N_EXPERIMENTS;

if normalize
    vec=[READOUT_COLS STIMULI_COLS];
    for i=1:length([READOUT_COLS ])
       DATA.data(:,vec(i))=(DATA.data(:,vec(i))-min(DATA.data(:,vec(i))))/( max(DATA.data(:,vec(i)))-min(DATA.data(:,vec(i))) );
    end
end


for iexp=1:N_EXPERIMENTS
    
    index=DATA.data(:,EXPERIMENT_COL)==EXPERIMENTS(iexp);
    TIME=DATA.data(index,TIME_COL);
    exps.n_obs{iexp} = length(READOUT_COLS);
    exps.exp_type{iexp} = 'fixed';
    
    st_names=regexprep(DATA.textdata(READOUT_COLS),'READOUT_','');
    obs_names = regexprep(DATA.textdata(READOUT_COLS),'READOUT_','');
    
    for jobs=1:length(obs_names)
        obs_names{jobs}=[obs_names{jobs} '_o'];
    end
    
    exps.obs_names{iexp}=char(obs_names);
    exps.obs{iexp}=[];
    
    for jobs=1:exps.n_obs{iexp}
        
        exps.obs{iexp} =[exps.obs{iexp}; [obs_names{jobs} '=' st_names{jobs}]];
        
    end
   
    exps.t_f{iexp} = TIME(end);
    exps.n_s{iexp} = length(TIME);
    exps.t_s{iexp} = TIME';
    exps.u_interp{iexp} = 'step';
    exps.t_con{iexp} = TIME';
    exps.n_steps{iexp} = length(TIME)-1;
    exps.u{iexp} = DATA.data(index,STIMULI_COLS)';
    exps.data_type = 'real';
    exps.noise_type = 'homo'; 
    exps.exp_data{iexp} = DATA.data(index,READOUT_COLS);
    exps.error_data{iexp} = DATA.data(index,STD_COLS);
    exps.exp_y0{iexp} = exps.exp_data{iexp}(1,:);

end

end

