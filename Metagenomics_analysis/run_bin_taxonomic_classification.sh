#!/usr/bin/evn bash

### coverage and bining 

bin_taxonomic_classification_main(){
   check_and_install
   run_binclassification
}



# check and install missing packages for the qc pipeline
check_and_install(){
   echo "Checking and installing packages for bin taxonomic classification"
   install_kraken
   download_krakendb
   echo "DONE checking and installing packages for bin taxonomic classification!"
}


install_kraken(){
   echo "Checking and installing kraken"
   if [ -d "${TOOLS_FOLDER}/kraken2-2.0.8-beta" ]; then
      echo "kraken already installed"
   else
      #Installing kraken 
      echo "Installing kraken"
      cd ${TOOLS_FOLDER}
      wget http://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz
      tar xzf v2.0.8-beta.tar.gz
      cd kraken2-2.0.8-beta
      sh install_kraken2.sh ${TOOLS_FOLDER}/kraken2-2.0.8-beta
      cd ..
      rm v2.0.8-beta.tar.gz
   fi
   echo "DONE checking and installing kraken!"
}

download_krakendb(){
   echo "Checking and downloading kraken database"
   if [ -d "${REFERENCE_FOLDER}/reference_database/kraken2DB" ]; then
      echo "kraken database already installed"
   else
      mkdir -p ${REFERENCE_FOLDER}/reference_database
      ln -s ${LINKPATH_DB}/reference_database/kraken2DB ${REFERENCE_FOLDER}/reference_database/kraken2DB
   fi
   echo "DONE checking and downloading kraken database!"
}

run_binclassification(){
   echo "Running Bin Classification"
   
   mkdir -p ${ANALYSIS_FOLDER}/bin_taxonomic_classification/kraken_out/metabat
   mkdir -p ${ANALYSIS_FOLDER}/bin_taxonomic_classification/kraken_out/maxbin
   # Run gene prediction using prodigal on
   #finallist=$(ll -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')

   echo "Running Bin Classification on Metabat"
   metabatbinlist=$(ls -d ${ANALYSIS_FOLDER}/binning/metabat/*.fa | awk '{print $NF}')
   for bin in ${metabatbinlist}
   do
      bin_name=$(basename ${bin%.*})
      ${TOOLS_FOLDER}/kraken2-2.0.8-beta/kraken2 \
      --db ${REFERENCE_FOLDER}/reference_database/kraken2DB \
      --threads 16 \
      --use-names \
      ${bin} \
      --report ${ANALYSIS_FOLDER}/bin_taxonomic_classification/kraken_out/metabat/${bin_name}.kreport > \
      ${ANALYSIS_FOLDER}/bin_taxonomic_classification/kraken_out/metabat/${bin_name}.kraken
   done

   echo "Running Bin Classification on Maxbin"
   maxbinlist=$(ls -d ${ANALYSIS_FOLDER}/binning/maxbin/*.fasta | awk '{print $NF}')
   for bin in ${maxbinlist}
   do
      bin_name=$(basename ${bin%.*})
      ${TOOLS_FOLDER}/kraken2-2.0.8-beta/kraken2 \
      --db ${REFERENCE_FOLDER}/reference_database/kraken2DB \
      --threads 16 \
      --use-names \
      ${bin} \
      --report ${ANALYSIS_FOLDER}/bin_taxonomic_classification/kraken_out/maxbin/${bin_name}.kreport > \
      ${ANALYSIS_FOLDER}/bin_taxonomic_classification/kraken_out/maxbin/${bin_name}.kraken
   done

   echo "DONE running bin classification!"
}

bin_taxonomic_classification_main
