# Define $PATH for generating the Rmarkdown document output
export RSTUDIO_PANDOC=/usr/lib/rstudio/bin/pandoc
# Generate reports for individual analysis for each sample
R -e "rmarkdown::render('CK121_CD24_individual_analysis.Rmd',output_file='CK121_CD24_individual_analysis.html')"
R -e "rmarkdown::render('CK120_CD13_individual_analysis.Rmd',output_file='CK120_CD13_individual_analysis.html')"
R -e "rmarkdown::render('CK119_organoid_individual_analysis.Rmd',output_file='CK119_organoid_individual_analysis.html')"
R -e "rmarkdown::render('CK5_organoid_individual_analysis.Rmd',output_file='CK5_organoid_individual_analysis.html')"
R -e "rmarkdown::render(input='CK120_CK121_GO.Rmd',output_file='CK120_CK121_GO.pdf')"
# Generate main figures with a merged analysis (not integrated, but the four samples in same color scheme and scale)
R -e "rmarkdown::render('merged_analysis.Rmd',output_file='merged_analysis.html')"
