clear;
close all;
plot(0);
ylim([0 1]);
xlim([0 100]);

disp('Use right click to stop draw');

global mim_line_x
global min_line_y;

[myobj,mim_line_x,min_line_y] = freehanddraw(gca,'color','r','linewidth',3); 

SYNBAD_Design_SO('min_line')

load RESULTS_DESIGN

eval(sprintf('min_line'))
inputs.simulate.var_circuit = results.xbest;
SYNBAD_Simulate(inputs)
