function [f,g] = min_line_obj(vy,opstr,xs,ys)


global mim_line_x;
global min_line_y;

idx = opstr.idx;
par = opstr.par;

n_real = opstr.n_real;
n_int = opstr.n_int;

x = vy(1:n_real);
y_int = vy(n_real+1:n_real+n_int);
y_bin = vy(n_real+n_int+1:end);
y = [y_int;y_bin];

k = par.value;

for ii=1:1:size(idx,2)
    k(idx{ii})=x(ii);
end

D_max = opstr.D_max;
D_min = opstr.D_min;

eval(sprintf('states= %s;',opstr.def_states));

z0=states.z0;
parstr.k = k;
parstr.y = y;

tspan = 0:0.1:1000;
parstr.u(1) = 40;
parstr.u(2) = 0;

[t,Y]=SYNBAD_CVODES('HL_odefile_c',tspan,z0,parstr,opstr.ivpsol_rtol,opstr.ivpsol_atol,Inf);

f=sum((min_line_y-interp1(t,Y(:,1),mim_line_x)).^2);

if isnan(f)
    f = 1e10;
end

g(1) = +D_max - sum(y);
g(2) = -D_min + sum(y);





