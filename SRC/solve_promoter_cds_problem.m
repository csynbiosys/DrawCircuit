function [promoterCDSmatrix] = solve_promoter_cds_problem(promoters,CDSs,exp_data,error_data,TIME)

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

    eqn=['d' CDSs{i} '='  ] ;
    
    states={states{:} CDSs{i}};

    for j=1:length(promoters)
        eqn=sprintf('%s\n',eqn);
        counter=counter+1;
        eqn =[eqn sprintf('\t%s',[ '+  Y' num2str(counter) ' *' promoters{j}])];
    end
    
    eqn=[eqn '-' transcriptDegradationParameter '*' CDSs{i}];
    
   
    eqns={eqns{:} eqn};
    
end

for i=1:length(CDSs)
    
    eval(CDSs{i});
    
    for j=1:length(codes)
        
        if strcmp(codes{j},'protein')
            
            eqn=['d' names{j} '='  translationParameter '*' CDSs{i}] ;
            states={states{:} names{j}};
            eval(names{j});
            parNames={parNames{:} parNamesCDS{:}};
            parValues=[parValues parValuesCDS];
            parNames={parNames{:} parName};
            parValues=[parValues parValue];
            eqn=[eqn '-' parName '*' names{j}];
               
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
    parNames={parNames{:} ['Y' num2str(i)]};
end

eqns=eqns';
eqns={promEqns{:} eqns{:}};
states=states';

workdir;

%% Basic model definition
AMIGOModel.n_st =length(states);
AMIGOModel.n_par = length(parValues);
AMIGOModel.n_stimulus = 0;
AMIGOModel.stimulus_names = '';
AMIGOModel.eqns = char(eqns);
AMIGOModel.st_names =char(states);
AMIGOModel.par = parValues;
AMIGOModel.par_names = char(parNames');

model_name='CDSPROBLEM';

obs_names={};

for i=1:length(states)
   obs_names{i}=[states{i} '_obs'];
   obs{i}=[obs_names{i} '=' states{i}];
end

%% Model execution
AMIGOModel.names_type='custom';
AMIGOModel.AMIGOsensrhs=0;
AMIGOModel.odes_file=fullfile(draw_circuit_work_dir,'TEMP',[model_name '.c']);
AMIGOModel.mexfile=fullfile(draw_circuit_work_dir,'TEMP',[model_name 'costMex']);
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
exps.n_obs{iexp} = length(states);
exps.exp_type{iexp} = 'fixed';
exps.exp_data{iexp}=exp_data;
exps.error_data{iexp}=exp_data;
exps.obs_names{iexp}=char(obs_names);
exps.obs{iexp}=char(obs);
exps.t_f{iexp} = TIME(end);
exps.n_s{iexp} = length(TIME);
exps.t_s{iexp} = TIME;

exps.data_type = 'real';
exps.noise_type = 'homo';

exps.exp_y0{iexp} = exp_data(1,:);

inputs.PEsol.id_global_theta='none';
inputs.model=AMIGOModel;
inputs.exps=exps;

[inputs privstruct]=AMIGO_Prep(inputs);

inputs.ivpsol.nthread=1;
problem.x_0= zeros(1,NBINS);
problem.x_L= zeros(1,NBINS);
problem.x_U= ones(1,NBINS);
problem.int_var=NBINS;
problem.neq=N_PROMOTERS;
opts.maxeval=10000;
algorithm='eSS';
problem.f='promoter_cds_obj';
opts.local.solver='misqp';

[Results]=MEIGO(problem,opts,algorithm,inputs,privstruct,NCONT,NBINS,N_PROMOTERS,N_CDS);
promoterCDSmatrix=reshape(Results.xbest,N_PROMOTERS,N_CDS)';


end
