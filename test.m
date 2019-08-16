promoters=get_promoter_list();
CDSs=get_CDS_list();
RBSs=get_RBS_list();

gen_promoter_cds_problem(promoters,CDSs,{RBSs{3},RBSs{3},RBSs{3}})