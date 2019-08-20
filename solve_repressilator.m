
%% Check our local REPOSITORY for available promoters and coding sequences
promoters=get_promoter_list();
CDSs=get_CDS_list();

%% Load data 
load('repressilator_data.mat')

%% Plot the known promoter coding sequence structure
subplot(1,2,1);
imagesc([0 0 1;1 0 0; 0 1 0]);
title('True solution');
yticks([1 2 3]);
yticklabels(promoters);
xticks([1 2 3]);
xticklabels(CDSs);
xlabel('CDS/PROTEIN/REPRESSOR');
ylabel('PROMOTER');

%% solve the problem
[promoterCDSmatrix]=solve_promoter_cds_problem(promoters,CDSs,exp_data,error_data,TIME);

%% Plot solution found by solver
subplot(1,2,2);
imagesc(promoterCDSmatrix);
title('Solution found by solver');
yticks([1 2 3]);
yticklabels(promoters);
xticks([1 2 3]);
xticklabels(CDSs);
xlabel('CDS/PROTEIN/REPRESSOR');
ylabel('PROMOTER');

