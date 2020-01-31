### TITLE : Make file for reproducible results
### AUTHOR : Javier Perales-Paton - javier.perales@bioquant.uni-heidelberg.de

## Please note: this script is supposed to be run in the root folder of the repository

# Store the root folder of the project
echo -e "[INFO] : Store current working directory as root of the project.";
PWD_project="$PWD"

# Define $PATH for generating the Rmarkdown document output
export RSTUDIO_PANDOC=/usr/lib/rstudio/bin/pandoc

# Get all bash scripts within the project to be run one by one
echo -e "[INFO] : Finding bash scripts to build markdowns for reproducible analysis.";
runs=`find $PWD_project/Individual_analysis_* -name run_rmd.sh`;

# Run one by one all the 
echo -e "[INFO] : Running individual analysis...";
for run_sh in "${runs[@]}";do
	echo -e "[INFO] \t Starting with '$(dirname $(basename $run_sh))' at $(date)";
	cd $(dirname $run_sh);
	bash $(basename $run_sh);
	cd $PWD_project;
done

# Finally run the merged analysis
echo -e "[INFO] : Creating figures for manuscript in a merged analysis with the four samples.";
cd Merged_analysis;
bash run_rmd.sh

# Completed
echo -e "[INFO] : Finished! ";
echo -e "[INFO] : \t You can find reports of every step in the analysis in every *.md genereted file.";
echo -e "[INFO] : \t Final figures are included at './Merged_analysis/output/figs'";

