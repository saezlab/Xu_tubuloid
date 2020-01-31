# Define $PATH for generating the Rmarkdown document output
export RSTUDIO_PANDOC=/usr/lib/rstudio/bin/pandoc
# Generate reports for individual analysis for each sample
R -e "rmarkdown::render('1_initial_clustering_CK120_CD13.Rmd',output_file='1_initial_clustering_CK120_CD13.md')"
R -e "rmarkdown::render('2_cell_assignment_CK120_CD13.Rmd',output_file='2_cell_assignment_CK120_CD13.md')"
R -e "rmarkdown::render('3_refine_clustering_CK120_CD13.Rmd',output_file='3_refine_clustering_CK120_CD13.md')"
R -e "rmarkdown::render('4_final_assignment_CK120_CD13.Rmd',output_file='4_final_assignment_CK120_CD13.md')"
