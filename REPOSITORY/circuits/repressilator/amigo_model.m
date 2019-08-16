%% Basic model definition
AMIGOModel.n_st = 6;
AMIGOModel.n_par = 0;
AMIGOModel.n_stimulus = 0;
AMIGOModel.stimulus_names = {};
AMIGOModel.eqns = '';
AMIGOModel.st_names = ['x1';'x2';'x3';'x1';'x2';'x3'];
AMIGOModel.par = [];
AMIGOModel.par_names = '';

model_name='repressilator';
%% Model execution
AMIGOModel.names_type='custom';
AMIGOModel.AMIGOsensrhs=0;
AMIGOModel.odes_file=[model_name '.c'];
AMIGOModel.mexfile=[model_name 'CostMex'];
AMIGOModel.exe_type='costMex';
AMIGOModel.overwrite_model=0;
AMIGOModel.compile_model=1;
AMIGOModel.input_model_type = 'charmodelC';
AMIGOModel.names_type = 'custom';
AMIGOModel.cvodes_include=[];
AMIGOModel.debugmode=0;
AMIGOModel.shownetwork=0;

iexp=1;
exps.n_exp=1;
exps.n_obs{iexp} = 0;
exps.exp_type{iexp} = 'fixed';
exps.obs_names{iexp}='';
exps.obs{iexp}='';
TIME=0:10:200;
exps.t_f{iexp} = TIME(end);
exps.n_s{iexp} = length(TIME);
exps.t_s{iexp} = TIME;

% exps.u_interp{iexp} = 'step';
% exps.t_con{iexp} = TIME';
% exps.n_steps{iexp} = length(TIME)-1;
% exps.u{iexp} = DATA.data(index,STIMULI_COLS)';
exps.data_type = 'real';
exps.noise_type = 'homo';
exps.exp_data{iexp} = nan(length(TIME),3);
exps.error_data{iexp} = nan(length(TIME),3);
exps.exp_y0{iexp} = [0 0 0 0 20 0];
inputs.PEsol.id_global_theta='none';

inputs.model=AMIGOModel;
inputs.exps=exps;
[inputs privstruct]=AMIGO_Prep(inputs);

feval(inputs.model.mexfunction,'sim_CVODES');