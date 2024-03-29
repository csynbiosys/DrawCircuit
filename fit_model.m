function[inputs,privstruct,res_ssm]=fit_model(inputs,privstruct)

vguess= inputs.PEsol.global_theta_guess;

inputs.ivpsol.nthread=8;

problem.x_0=vguess;
problem.x_U=inputs.PEsol.global_theta_max;
problem.x_L=inputs.PEsol.global_theta_min;

problem.f=@davidObj;

inputs.nlpsol.eSS.dim_refset=2;
inputs.nlpsol.eSS.ndiverse=100;
inputs.nlpsol.eSS.local.n1=1;
inputs.nlpsol.eSS.local.n2=3;
inputs.nlpsol.eSS.maxeval=10000;
inputs.nlpsol.eSS.maxtime=inf;
inputs.nlpsol.eSS.local.solver='fmincon';
inputs.nlpsol.eSS.local.finish=0;

[res_ssm]=ess_kernel(problem,inputs.nlpsol.eSS,inputs,privstruct);

x=res_ssm.xbest;

inputs.model.par(inputs.PEsol.index_global_theta)=x;

end

function [f,g,r]=davidObj(x,inputs,privstruct,INITIALU)

inputs.model.par(inputs.PEsol.index_global_theta)=x(1:length(inputs.PEsol.index_global_theta));

feval(inputs.model.mexfunction,'cost_LSQ');

f=outputs.f;
r=outputs.w_res;

g=[];

if(f==0),f=Inf;end

for iexp=1:length(inputs.exps.exp_data)
   if outputs.sim_stats{iexp}.flag<0
      f=1e10;
   end
end

end
