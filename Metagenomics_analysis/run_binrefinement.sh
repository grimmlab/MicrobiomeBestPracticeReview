#!/usr/bin/evn bash

### binrefinement


binrefinement_main(){
   check_and_install
   run_binrefiner
}



# check and install missing packages for the qc pipeline
check_and_install(){
   echo "Checking and installing packages for binrefinement"
   install_binrefiner
   install_checkm
   echo "DONE checking and installing packages for binrefinement!"
}


run_binrefiner(){
   echo "Running Binrefiner"

   mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/binrefiner_output/
   mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/input/maxbin
   mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/input/metabat

   cp ${ANALYSIS_FOLDER}/binning/maxbin/*.fasta ${ANALYSIS_FOLDER}/bin_refinement/input/maxbin
   cp ${ANALYSIS_FOLDER}/binning/metabat/*.fa ${ANALYSIS_FOLDER}/bin_refinement/input/metabat

   cd  ${ANALYSIS_FOLDER}/bin_refinement/
   #Binning_refiner \
   -i ${ANALYSIS_FOLDER}/bin_refinement/input/ \
   -p BR

   #mv BR* ${ANALYSIS_FOLDER}/bin_refinement/binrefiner_output/

   #running CheckM
   mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/checkM/maxbin/
   #running for maxbin
   python2 /usr/local/bin/checkm lineage_wf \
   -t 8 \
   -x fasta ${ANALYSIS_FOLDER}/binning/maxbin/ \
   ${ANALYSIS_FOLDER}//bin_refinement/checkM/maxbin/

   python2 /usr/local/bin/checkm taxonomy_wf \
   domain \
   Bacteria \
   ${ANALYSIS_FOLDER}/binning/maxbin/ \
   ${ANALYSIS_FOLDER}/bin_refinement/checkM/maxbin/ \
   -x fasta

   mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/checkM/maxbin/plots/
   python2 /usr/local/bin/checkm bin_qa_plot \
   -x fasta \
   ${ANALYSIS_FOLDER}/bin_refinement/checkM/ \
   ${ANALYSIS_FOLDER}/binning/maxbin/ \
   ${ANALYSIS_FOLDER}/bin_refinement/checkM/maxbin/plots


   mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/checkM/metabat
   #runrning for metabat
   python2 /usr/local/bin/checkm lineage_wf \
   -t 8 \
   -x fa ${ANALYSIS_FOLDER}/binning/metabat/ \
   ${ANALYSIS_FOLDER}/bin_refinement/checkM/metabat

   mkdir -p ${ANALYSIS_FOLDER}/bin_refinement/checkM/metabat/plots/
   python2 /usr/local/bin/checkm taxonomy_wf \
   domain \
   Bacteria \
   ${ANALYSIS_FOLDER}/binning/metabat/ \
   ${ANALYSIS_FOLDER}/bin_refinement/checkM/metabat/ \
   -x fasta

   python2 /usr/local/bin/checkm bin_qa_plot \
   -x fasta \
   ${ANALYSIS_FOLDER}/bin_refinement/checkM/ \
   ${ANALYSIS_FOLDER}/binning/metabat/ \
   ${ANALYSIS_FOLDER}/bin_refinement/checkM/metabat/plots

   echo "DONE running Binrefiner!"

}

install_binrefiner(){
   echo "Checking and installing binrefiner"

   M_TMP=$(pip list | grep "Binning-refiner")
   if [ -z "$M_TMP" ]; then 
      echo "Installing binrefiner"
      pip install Binning_refiner
   else
      echo "binrefiner already installed"
   fi

   echo "DONE checking and installing binrefiner!"
}

install_checkm(){
   echo "Checking and installing checkm"

   M_TMP=$(pip list | grep "checkm-genome")
   if [ -z "$M_TMP" ]; then 
      echo "Installing checkm"

      pip install checkm-genome

      mkdir -p $REFRENCE_FOLDER/checkm
      cd $REFRENCE_FOLDER/checkm
      wget https://data.ace.uq.edu.au/public/CheckM_databases/
      tar -zxvf checkm_data_2015_01_16.tar.gz
      sudo python2 /usr/local/bin/checkm data setRoot $REFRENCE_FOLDER/checkm

      cd $TOOLS_FOLDER
      wget https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip
      unzip pplacer-Linux-v1.1.alpha19.zip
   else
      echo "checkm already installed"
   fi

   export PATH=${TOOLS_FOLDER}/hmmer-3.2.1/:$PATH
   export PATH=${TOOLS_FOLDER}/pplacer-Linux-v1.1.alpha19/:$PATH
   export PATH=${TOOLS_FOLDER}/Prodigal-2.6.1/:$PATH
   echo "DONE checking and installing checkm!"
}

binrefinement_main
