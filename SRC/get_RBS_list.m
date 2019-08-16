function available_CDS = get_RBS_list()
    
    workdir;
    available_CDS=what(fullfile(work_dir,'REPOSITORY/parts/RBS/'));
    available_CDS=available_CDS.m;
    available_CDS=strrep(available_CDS,'.m','');

end