cd GSEA_4.3.2
gsea-cli.bat GSEAPreranked -gmx ftp.broadinstitute.org://pub/gsea/msigdb/_species_/gene_sets/_gene_set_ -collapse _collapse_ -mode Abs_max_of_probes -norm meandiv -nperm 1000 -rnd_seed timestamp -rnk _rnk_file_ -scoring_scheme weighted -rpt_label _the_label_ -chip ftp.broadinstitute.org://pub/gsea/msigdb/mouse/annotations/_chip_ -create_svgs false -include_only_symbols true -make_sets true -plot_top_x _numberplot_ -set_max 500 -set_min 15 -zip_report false -out _save_dir_
