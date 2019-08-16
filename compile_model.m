function [inputs, privstruct]=compile_model()
normalize=1;
exps=load_data('DATA/experimental_data_loop_1.csv',normalize);

[AMIGOModel,index_k index_n, index_tau, index_w, index_diff]=load_model();

try
    pars=importdata('pars.csv');
    AMIGOModel.par=pars';
catch
end

inputs.PEsol.global_theta_guess=AMIGOModel.par;
inputs.PEsol.global_theta_min(index_k)=0.005;
inputs.PEsol.global_theta_max(index_k)=1;
inputs.PEsol.global_theta_min(index_n)=0.1;
inputs.PEsol.global_theta_max(index_n)=60;
inputs.PEsol.global_theta_min(index_tau)=1e-4;
inputs.PEsol.global_theta_max(index_tau)=1e-1;
inputs.PEsol.global_theta_min(index_w)=0;
inputs.PEsol.global_theta_max(index_w)=1;

inputs.PEsol.global_theta_max(index_diff)=10;
inputs.PEsol.global_theta_min(index_diff)=1e-5;

inputs.PEsol.id_global_theta=AMIGOModel.par_names;
inputs.ivpsol.atol=1e-6;
inputs.ivpsol.rtol=1e-6;

inputs.model=AMIGOModel;

inputs.exps=exps;
[inputs,privstruct]=AMIGO_Prep(inputs);

end
