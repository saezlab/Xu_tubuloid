export RSTUDIO_PANDOC=/usr/lib/rstudio/bin/pandoc
R -e "rmarkdown::render('CK121_CD24_individual_analysis.Rmd',output_file='CK121_CD24_individual_analysis.html')"
R -e "rmarkdown::render('CK120_CD13_individual_analysis.Rmd',output_file='CK120_CD13_individual_analysis.html')"
R -e "rmarkdown::render('CK119_organoid_individual_analysis.Rmd',output_file='CK119_organoid_individual_analysis.html')"
R -e "rmarkdown::render('CK5_organoid_individual_analysis.Rmd',output_file='CK5_organoid_individual_analysis.html')"
