%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TITLE: Find circuit with optimal switch-like behaviour starting from a
% library of 4 transcripts and 8 promoters. Hill Kinetics. 

% REF:  Otero-Muras I and Banga JR 2014, BMC Systems Biology

% Single Objective
% Constrained 3D problem. 
% Real and Binary variables. 
% Solver: ESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%================================
% MIXED INTEGER MODEL FRAMEWORK
%================================
inputs.model.lib_type ='HL_Library';  %Choose 'MA_Library'|'HL_Library'|'CN_Library'
inputs.model.ode_name ='HL_odefile_c';
inputs.model.n_integer_var = 0;
inputs.model.n_real_var = 2;
inputs.model.n_binary_var = 32;
inputs.model.def_param_function= 'HL_default_parameters_1';
inputs.model.def_state_function= 'HL_default_states_1';
inputs.model.transc_promot_function= 'HL_transcripts_and_promoters';
inputs.model.u_values=[40;0];

%============================
% DESIGN PROBLEM OPTIONS
%============================S
inputs.design.idx ={1:20};
%inputs.design.par_x={{'kdeg_tetR','kdeg_lacI','kdeg_cI','kdeg_araC'},{'kf_tetRaTc'}};


inputs.design.var_L = [zeros(1,20) zeros(1,32)];
inputs.design.var_U = [10.*ones(1,20) ones(1,32)];
inputs.design.var_0 = zeros(1,52);
inputs.design.D_max = 5;  % only applies in MITS, ESS, ACO
inputs.design.D_min = 5;  % only applies in MITS, ESS, ACO



% SINGLE OBJECTIVE 
%============================
inputs.design.objective='min_line_obj';
    

%============================
% SIMULATE OPTIONS
%============================ 
%inputs.simulate.var_circuit =   [0.05 0.1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0]; 
inputs.simulate.var_circuit =   [0.0277598	0.0401563	0.000962291	0.136446	0.0335545	0.000962558	0.00731231	0.327307	0.0165788	4.99E-18	0.00627083	0.00154912	0.0187579	0.279542	0.00243506	0.521152	0.190127	0.372098	0.00108568	0.0838595	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	1	1	0	0	0	0];
inputs.simulate.tspan = 0:1:1000;
inputs.simulate.objective={'min_line_obj'};

%==================================
% MINLP SOLVER OPTIONS
%==================================
inputs.optsol.optsolver='ESS';       % Choose MINLP solver 'ESS'|'MITS'|'ACO'|'VNS'
inputs.optsol.maxtime = 50;
inputs.optsol.maxeval = [];
inputs.optsol.ess.local.solver = 'misqp';

%==================================
% IVP SOLVER OPTIONS
%==================================
inputs.ivpsol.rtol=1.0D-7;           % [] IVP solver integration tolerances
inputs.ivpsol.atol=1.0D-7;


