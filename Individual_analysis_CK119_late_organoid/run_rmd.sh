# Define $PATH for generating the Rmarkdown document output
export RSTUDIO_PANDOC=/usr/lib/rstudio/bin/pandoc
# Generate reports for individual analysis for each sample
R -e "rmarkdown::render('1_initial_clustering_CK119_organoid.Rmd',output_file='1_initial_clustering_CK119_organoid.md')"
R -e "rmarkdown::render('2_cell_assignment_CK119_organoid.Rmd',output_file='2_cell_assignment_CK119_organoid.md')"
