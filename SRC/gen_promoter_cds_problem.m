function model = gen_promoter_cds_problem(promoters,CDSs,RBSs,expData,problem_name)

N_PROMOTERS=length(promoters);
N_CDS=length(CDSs);

eqns={};
promEqns={};
parNames={};
parValues=[];
states={};

for i=1:length(promoters)
    
    eval(promoters{i});
    promEqns={promEqns{:} [promoters{i} '=' transcriptionFunction ';']};
    parNames={parNames{:},parsNamesTranscription{:}};
    parValues=[parValues parsValsTranslation];
    
end

promEqns=promEqns';

counter=0;
for i=1:length(CDSs)
    
    eval(CDSs{i});
    eval(RBSs{i})
    eqn=['d' CDSs{i} '='  ] ;
    
    states={states{:} CDSs{i}};
    
    parNames={parNames{:},['translationEffiency_' promoters{i}]};
    parValues=[parValues translationEffiency];
    
    for j=1:length(promoters)
        eqn=sprintf('%s\n',eqn);
        counter=counter+1;
        eqn =[eqn sprintf('\t%s',[ '+ translationEffiency_' promoters{i} '* Y' num2str(counter) ' *' promoters{i}])];
    end
    
    eqn=[eqn ';'];
    eqns={eqns{:} eqn};
    
end

counter=0;
for i=1:length(CDSs)
    
    eval(CDSs{i});
    
    for j=1:length(codes)
        
        if strcmp(codes{j},'protein')
            
            eqn=['d' names{j} '='  translationParameter '*' CDSs{i}] ;
            states={states{:} names{j}};
            eval(names{j});
            parNames={parNames{:} parNamesCDS{:}};
            parValues=[parValues parValuesCDS];
            eqn=[eqn '-' parName '*' transcriptDegradationParameter];
               
        else
            
            stop('Invalid sequece type!');
            
        end
        
    end
    
    eqn=[eqn ';'];
    eqns={eqns{:} eqn};
end

NBINS=counter;
NCONT=length(parValues);
parValues=[parValues zeros(1,NBINS)];

for i=1:NBINS
    parNames={parNames{:} ['Y' numstr(i)]};
end

eqns=eqns';
eqns={promEqns{:} eqns{:}};
states=states';

model.parNames=parNames;
model.parValues=parValues;
model.eqns=eqns;
model.states=states;

%% Basic model definition
AMIGOModel.n_st =length(states);
AMIGOModel.n_par = length(parValues);
AMIGOModel.n_stimulus = 0;
AMIGOModel.stimulus_names = '';
AMIGOModel.eqns = char(eqns);
AMIGOModel.st_names =char(states);
AMIGOModel.par = parValues;
AMIGOModel.par_names = char(parNames');
model_name='PROBLEM';

%% Model execution
AMIGOModel.names_type='custom';
AMIGOModel.AMIGOsensrhs=0;
AMIGOModel.odes_file=[model_name '.c'];
AMIGOModel.mexfile=[model_name 'CostMex'];
AMIGOModel.exe_type='costMex';
AMIGOModel.overwrite_model=1;
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


end
