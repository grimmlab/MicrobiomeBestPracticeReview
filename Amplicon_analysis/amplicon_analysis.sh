#!/usr/bin/env bash

### Workflow for amplicon analysis for bacterial data

# Root folder name"
NAME=example_run

# Raw data folder path
SRC_RAWDATA='/home/richa/Desktop/MicrobiomeBestPracticeReview/Amplicon_analysis'
#LINKPATH_DB='/data1/Active_Projects/paper_scripts/reference/'

amplicon_analysis_main(){
   create_folders
   set_variables # -> Never comment this function
   fetch_example_data # -> Uncomment this function if you want to run it on an example data
   copy_rawdata
   download_reference_database
   run_mothur_workflow
   run_dada2_workflow
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
   export BIN_FOLDER=$(pwd)/bin
   #export LINKPATH_DB=$LINKPATH_DB

   export SRC_RAWDATA=$NAME/example_data
   export CURRDIR=$(pwd)
   echo "DONE setting variables for paths!"
}

fetch_example_data(){
   mkdir -p $NAME/example_data
   mkdir -p $NAME/metadata   
   cd $NAME/example_data
   wget https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip
   #wget www.mothur.org/w/images/d/d6/MiSeqSOPData.zip
   echo "Richa is testing"
   unzip miseqsopdata.zip
   cd MiSeq_SOP
   rm Mock_S280_*
   rm HMP_MOCK.v35.fasta
   mv *.fastq ../
   cp *.metadata ../../metadata
   cd ..
   rm MiSeq_SOP -rf
   rm -rf miseqsopdata.zip
   rm -rf __MACOSX/ 
   export SRC_RAWDATA=$PWD/

   # Download test data
   #To test this we are using public available dataset from Multi-omics differentially Classify Disease State and Treatment Outcome in Pediatric Crohn’s Disease. In this study total 115 sample are submitted put of which 40 are 16S rRNA samples and 75 are metagenomics samples (https://www.ebi.ac.uk/ena/data/view/PRJEB21933) https://www.ncbi.nlm.nih.gov//bioproject/PRJEB21933.</br>
   #We will use 10 these samples (5 treatment-naive CD and 5 control pediatric patients) as an example dataset for demonstrate MOTHUR 16S rRNA workflow </p>

   #Specs of the data:
   #- Human gut data
   #- Two conditions: Disease State and Treatment Outcome
   #- 16S rRNA form V6-V8 region
   #- Illumina MiSeq paired end data
   #~ 574 bp in length

   #Using only 5 samples for example purpose .....
   #ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR204/ERR2042042/S1_R1_001.fastq.gz
   #S=0
   #for FILE in ERR20420{42..46}
   #do
   #wget -r -l2 -A.gz ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR204/
   #S=$(($S+1))
   #wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR204/${FILE}/S${S}_R1_001.fastq.gz
   #wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR204/${FILE}/S${S}_R2_001.fastq.gz
   #done

}

# copy raw data from source folder to analysis folder structure
copy_rawdata(){
   lst=$(ls -d $SRC_RAWDATA/*.fastq)
   for file in $lst
   do
      echo "Copying ${file}"
      cp ${file} ${RAWDATA_FOLDER}/
   done
   echo "DONE copying rawdata!"
}

download_reference_database(){
   echo "Running download reference database"

   ## Download the references and taxonomy


   #Greengenes reference files can be downloaded from  [here](https://www.mothur.org/wiki/Greengenes-formatted_databases)
   cd ${REFERENCE_FOLDER}
   mkdir -p ${REFERENCE_FOLDER}/greengenes
   cd  ${REFERENCE_FOLDER}/greengenes
   wget https://mothur.s3.us-east-2.amazonaws.com/wiki/gg_13_8_99.taxonomy.tgz
   #wget http://www.mothur.org/w/images/6/68/Gg_13_8_99.taxonomy.tgz
   tar -xzvf gg_13_8_99.taxonomy.tgz
   rm gg_13_8_99.taxonomy.tgz

   cd $REFERENCE_FOLDER
   REFERENCE_VERSION=v132
   mkdir -p ${REFERENCE_FOLDER}/silva
   cd  ${REFERENCE_FOLDER}/silva
   #wget https://mothur.org/w/images/3/32/Silva.nr_${REFERENCE_VERSION}.tgz
   wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_${REFERENCE_VERSION}.tgz
   tar -xzvf silva.nr_${REFERENCE_VERSION}.tgz
   rm silva.nr_${REFERENCE_VERSION}.tgz

   # Silva database for DADA2
   cd ${REFERENCE_FOLDER}
   mkdir -p ${REFERENCE_FOLDER}/dada2
   cd  ${REFERENCE_FOLDER}/dada2
   wget https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz
   gunzip silva_nr_v132_train_set.fa.gz

   echo "DONE running download reference database"
}

# Run mothur analysis
run_mothur_workflow(){
   echo "Running mothur workflow"

   cd ${CURRDIR}
	
	 . ./mothur_workflow.sh
	
   echo "DONE running mothur workflow"
}

# Run dada2 analysis
run_dada2_workflow(){
   echo "Running dada2 workflow"

   cd ${CURRDIR}
   mkdir -p ${ANALYSIS_FOLDER}/dada2/outdir/plots/
	echo $PWD
   Rscript dada2_workflow.R \
     --input ${RAWDATA_FOLDER}/ \
	  --output ${ANALYSIS_FOLDER}/dada2/outdir/ \
	  --reference ${REFERENCE_FOLDER}/dada2/ \
     --plots ${ANALYSIS_FOLDER}/dada2/outdir/plots/
   echo "DONE running dada2 workflow"
}

amplicon_analysis_main
