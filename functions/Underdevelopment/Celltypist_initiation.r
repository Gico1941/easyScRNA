install_miniconda(path = miniconda_path(), update = TRUE, force = FALSE)
conda_create('celltypist')
conda_install('celltypist',list("scanpy","celltypist","pandas","numpy"))