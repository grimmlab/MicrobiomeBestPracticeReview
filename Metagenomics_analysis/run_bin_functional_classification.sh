#!/usr/bin/evn bash

### coverage and bining 

bin_functional_classification_main(){
   check_and_install
   run_binclassification
}



# check and install missing packages for the qc pipeline
check_and_install(){
   echo "Checking and installing packages for bin functinal classification"
   install_prokka
   echo "DONE checking and installing packages for bin functinal classification!"
}

run_binclassification(){
   echo "Running Bin Classification"
   
   mkdir -p ${ANALYSIS_FOLDER}/bin_functional_annotation/prokka_out/metabat
   mkdir -p ${ANALYSIS_FOLDER}/bin_functional_annotation/prokka_out/maxbin
   # Run gene prediction using prodigal on
   #finallist=$(ll -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')

   echo "Running prokka on maxbin"
   maxbinlist=$(ls -d ${ANALYSIS_FOLDER}/binning/maxbin/*.fasta | awk '{print $NF}')
   for i in $(ls ${maxbinlist})
   do
      bin_name=$(basename ${i%.*})

      ${TOOLS_FOLDER}/prokka/bin/prokka ${i} \
      --quiet \
      --cpus 24 \
      --outdir ${ANALYSIS_FOLDER}/bin_functional_annotation/prokka_out/maxbin/$bin_name \
      --prefix $bin_name \
      --metagenome \
      --kingdom 'Bacteria' \
      --locustag 'PROKKA' \
      --force
      
   done  

   echo "Running prokka on metabat"
   metabatbinlist=$(ls -d ${ANALYSIS_FOLDER}/binning/metabat/*.fa | awk '{print $NF}')
   for bin in ${metabatbinlist}
   do
      bin_name=$(basename ${i%.*})

      ${TOOLS_FOLDER}/prokka/bin/prokka ${i} \
      --quiet \
      --cpus 24 \
      --outdir ${ANALYSIS_FOLDER}/bin_functional_annotation/prokka_out/metabat/$bin_name \
      --prefix $bin_name \
      --metagenome \
      --kingdom 'Bacteria' \
      --locustag 'PROKKA' \
      --force
   done  

   echo "DONE running bin classification!"
}

install_prokka(){
   echo "Checking and installing prokka"
   if [ -d "${TOOLS_FOLDER}/prokka" ]; then
      echo "prokka already installed"
   else
      #Installing prokka 
      echo "Installing prokka"
      #prokka
      cd ${TOOLS_FOLDER}
      git clone https://github.com/tseemann/prokka.git 
      cd prokka/bin/
      prokka --setupdb
   fi

   echo "DONE checking and installing prokka!"
}

bin_functional_classification_main
