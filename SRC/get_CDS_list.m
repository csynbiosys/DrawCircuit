function available_CDS = get_CDS_list()
    
    workdir;
    available_CDS=what(fullfile(draw_circuit_work_dir,'REPOSITORY/parts/CDS/'));
    available_CDS=available_CDS.m;
    available_CDS=strrep(available_CDS,'.m','');

end