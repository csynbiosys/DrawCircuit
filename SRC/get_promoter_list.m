function available_promoters = get_promoter_list()
    
    workdir;
    available_promoters=what(fullfile(work_dir,'REPOSITORY/parts/PROM/'));
    available_promoters=available_promoters.m;
    available_promoters=strrep(available_promoters,'.m','');

end