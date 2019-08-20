
function [f,g,r]=promoter_cds_obj(x,inputs,privstruct,NCONT, NBINS,NPROMOTERS,NCDS)

inputs.model.par(NCONT+1:NCONT+NBINS)=x;
x=reshape(x,NPROMOTERS,NCDS)';
g=sum(x');

g(g==1)=0;
g(g==0)=0;
g(g>1)=g(g>1)-1;

feval(inputs.model.mexfunction,'cost_LSQ');

f=outputs.f;
r=outputs.w_res;

if(f==0)
    f=Inf;
end

for iexp=1:length(inputs.exps.exp_data)
   if outputs.sim_stats{iexp}.flag<0
      f=1e10;
   end
end

end

