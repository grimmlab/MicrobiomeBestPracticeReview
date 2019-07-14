#!/usr/bin/env bash

### Workflow for whole genome shotgun metagenomics analysis

# Root folder name"
NAME=Test_Metagenomic

# Raw data folder path
SRC_RAWDATA='/data1/Active_Projects/Metagenomic_QC/rawdata/'
LINKPATH_DB='/data1/Active_Projects/paper_scripts/reference/'



metagenomics_analysis_main(){
   create_folders
   set_variables # -> Never comment this function 
   copy_rawdata
   #fetch_example_data
   run_qc
   run_assembly
   run_coassembly
   run_reference_analysis
   run_comparative_analysis
   run_coverage_and_bining
   run_binrefinement
   run_bin_taxonomic_classification
   run_bin_functional_classification
   echo $LINKPATH_DB
}

create_folders(){

   echo "Creating sub-folders..."

   # Sub-folders in the root folder
   for FOLDER in analysis tools rawdata reference bin
   do
      mkdir -p $NAME/$FOLDER
   done
   echo "DONE creating sub-folders!"
}

# setting variable path
set_variables(){
   echo "Setting variables for paths..."

   export ROOT_FOLDER_NAME=$NAME
   export TOOLS_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/tools
   export RAWDATA_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/rawdata
   export ANALYSIS_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/analysis
   export REFERENCE_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/reference
   export BIN_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/bin
   export LINKPATH_DB=$LINKPATH_DB
   echo "DONE setting variables for paths!"
}


# copy raw data from source folder to analysis folder structure 
copy_rawdata(){
   lst=$(ls -d $SRC_RAWDATA/*.fastq.gz)
   for file in $lst 
   do
      echo "Copying ${file}"
      cp ${file} ${RAWDATA_FOLDER}/
   done
   echo "DONE copying rawdata!"
}


# fetch raw data from web to the analysis folder structure 
fetch_example_data(){
   lst=$(ls -d $SRC_RAWDATA/*.fastq.gz)

   cd ${RAWDATA_FOLDER}

   wget ftp://public-ftp.hmpdacc.org/Illumina/HMDEMO/SRP002430/plasma/affected/SRS072364.tar.bz2
   wget ftp://public-ftp.hmpdacc.org/Illumina/HMDEMO/SRP002430/plasma/affected/SRS072331.tar.bz2

   cd -
}

# Run QC analysis 
run_qc(){
   echo "Running Quality Control" 

   . ./run_qc.sh

   echo "DONE Quality Control!" 
}

# Run assembly 
run_assembly(){
   echo "Running Assembly" 

   . ./run_assembly.sh

   echo "DONE Running Assembly!" 
}

# Run assembly 
run_coassembly(){
   echo "Running Coassembly" 

   . ./run_coassembly.sh

   echo "DONE Running Coassembly!" 
}

run_reference_analysis(){
   echo "Running Reference Analysis" 

   . ./run_reference_analysis.sh

   echo "DONE Running Reference Analysis!" 
}

run_comparative_analysis(){
   echo "Running Comparative Analysis" 

   . ./run_comparative_analysis.sh

   echo "DONE Running Comparative Analysis!" 
}

run_coverage_and_bining(){
   echo "Running Coverage and Bining" 

   . ./run_coverage_bining.sh

   echo "DONE running Coverage and Bining!" 
}

run_binrefinement(){
   echo "Running Binrefinement" 

   . ./run_binrefinement.sh

   echo "DONE running Binrefinement!" 
}

run_bin_taxonomic_classification(){
   echo "Running Bin Taxonomic Classification" 

   . ./run_bin_taxonomic_classification.sh

   echo "DONE running Bin Taxonomic Classification!" 
}

run_bin_functional_classification(){
   echo "Running Bin Functional Classification" 

   . ./run_bin_functional_classification.sh

   echo "DONE running Bin Functional Classification!" 
}

metagenomics_analysis_main
