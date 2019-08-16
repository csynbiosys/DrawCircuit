function [AMIGOModel,index_k,index_n,index_tau,index_w, index_diff] = load_model()

DATA = importdata('MODEL/PROMOTER_LIST.csv');

N_PROMOTERS=size(DATA.data,1);
EQNS={};
INPUTS={};
PARAMETERS={};
PROMOTER_NAME={};
par_names={};
PARAMETER_NAMES={};
par_values=[];

index_k=[];
index_n=[];
index_tau=[];
index_w=[];
index_diff=[];
counter=1;

for PRO=1:N_PROMOTERS
    
    PROMOTER_NAME{PRO}=DATA.textdata{PRO+1,1};
    EQNS{PRO}=[PROMOTER_NAME{PRO} '=' DATA.textdata{PRO+1,3}];
    INPUTS{PRO}=split(DATA.textdata{PRO+1,2},' ');
    PARAMETERS{PRO}=DATA.data(PRO,:);
    
    for input=1:length(INPUTS{PRO})
        
        k_name=['k_' PROMOTER_NAME{PRO} '_' INPUTS{PRO}{input}];
        n_name=['n_' PROMOTER_NAME{PRO} '_' INPUTS{PRO}{input}];
        input_name=[INPUTS{PRO}{input}];
        str='(xxxx^nnnn*(kkkk^nnnn + 1))/(kkkk^nnnn + xxxx^nnnn)';
        str=strrep(str,'xxxx',[input_name]);
        str=strrep(str,'kkkk',k_name);
        str=strrep(str,'nnnn',n_name);
        EQNS{PRO}=strrep(EQNS{PRO},['NOT(' INPUTS{PRO}{input} ')'],['(1-' str ')']);
        PARAMETER_NAMES{PRO}={k_name,n_name};
        par_names{end+1}= k_name;
        par_names{end+1}= n_name;
        index_k=[index_k counter];
        counter=counter+1;
        index_n=[index_n counter];
        counter=counter+1;
        
    end
    par_values=[par_values  PARAMETERS{PRO}];
end


DATA = importdata('MODEL/INPUT_LIST.csv');
INPUT_NAMES=DATA.textdata;

TRANSPORT={};
TRANSPORT_STATES={};
for i =1:length(INPUT_NAMES)
    
    parName=[INPUT_NAMES{i} '_kdiff'];
    TRANSPORT{i}=['d' INPUT_NAMES{i} '=(' [INPUT_NAMES{i} 'i'] ' - '  INPUT_NAMES{i}  ') * ' parName];
    TRANSPORT_STATES{i}=[INPUT_NAMES{i}] ;
    par_values=[par_values  0.01];
    par_names=[par_names parName];
    index_diff=[index_diff counter];
    counter=counter+1;
    INPUT_NAMES{i}=[INPUT_NAMES{i} 'i'];
    
end


DATA = importdata('MODEL/EXPRESSED_LIST.csv');
EXPRESSED_NAMES=DATA.textdata;
N_EXPRESSED=length(EXPRESSED_NAMES);

STATE={};


for i=1:N_EXPRESSED
    
    tau_name=['tau_' EXPRESSED_NAMES{i}];
    ODES{i}='';
    
    for j=1:N_PROMOTERS
        
        w_name=['w_' PROMOTER_NAME{j} '_' EXPRESSED_NAMES{i}];
        
        ODES{i}=[ODES{i} w_name '*' PROMOTER_NAME{j}];
        par_values=[par_values  1];
        par_names=[par_names w_name];
        index_w=[index_w counter];
        counter=counter+1;
        
        if j==N_PROMOTERS
            
        else
            ODES{i}=[ODES{i} '+'];
        end
        
    end
    
    par_names{end+1}=tau_name;
    
    index_tau=[index_tau counter];
    counter=counter+1;
    par_values=[par_values  DATA.data(i)];
    
    ODES{i}=['(' ODES{i} '-' EXPRESSED_NAMES{i} ')*'   tau_name];
    ODES{i}=['d' EXPRESSED_NAMES{i} '=' ODES{i}];
    
end

%% Basic model definition
AMIGOModel.n_st = N_EXPRESSED+length(TRANSPORT_STATES);
AMIGOModel.n_par = length(par_names);
AMIGOModel.n_stimulus = length(INPUT_NAMES);
AMIGOModel.stimulus_names = char(INPUT_NAMES);
AMIGOModel.eqns = char({EQNS{:} ODES{:} TRANSPORT{:}}');
AMIGOModel.st_names = char({EXPRESSED_NAMES{:} TRANSPORT_STATES{:}});
AMIGOModel.par = par_values;
AMIGOModel.par_names = char(par_names');

model_name='DrawCircuit';
%% Model execution
AMIGOModel.names_type='custom';
AMIGOModel.AMIGOsensrhs=0;
AMIGOModel.odes_file=['MODEL/' model_name '.c'];
AMIGOModel.mexfile=['MODEL/' model_name 'CostMex'];
AMIGOModel.exe_type='costMex';
AMIGOModel.overwrite_model=1;
AMIGOModel.compile_model=1;
AMIGOModel.input_model_type = 'charmodelC';
AMIGOModel.names_type = 'custom';
AMIGOModel.cvodes_include=[];
AMIGOModel.debugmode=0;
AMIGOModel.shownetwork=0;


end