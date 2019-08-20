addpath(genpath(pwd));
fid=fopen('workdir.m','w');
fprintf(fid,'draw_circuit_work_dir=''%s'';',pwd);
fclose(fid);
