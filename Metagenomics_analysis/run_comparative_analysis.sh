#!/usr/bin/evn bash

### Comparative Analysis


comparative_analysis_main(){
   check_and_install
   contigs_classification_with_kraken
   comparative_functional_annotation_prokka
   gene_classification_with_diamond_megan
   predict_metabolic_pathway_minpath
}

# Check and install missing packages for the QC pipeline
check_and_install(){
   echo "Checking and installing packages for Comparative Analysis"
   install_diamond
   download_nrdb
   install_kraken
   download_krakendb
   install_prokka
   install_megan
   download_megandb
   install_minpath
   echo "DONE checking and installing packages for Comparative Analysis!"
}


contigs_classification_with_kraken(){

   echo "Running contigs classification with kraken"
   #run_kranken

   # running kraken
   mkdir -p ${ANALYSIS_FOLDER}/taxonomic_classification/kranken_output
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      ${TOOLS_FOLDER}/kraken2-2.0.8-beta/kraken2 \
      --db ${REFERENCE_FOLDER}/reference_database/kraken2DB \
      --threads 16 \
      --use-names \
      ${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs_1000_filtered.fasta \
      --report ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/kranken_output/$(basename ${i%.1*}.kreport) > \
      ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/kranken_output/$(basename ${i%.1*}.kraken)
   done

   echo "DONE running contigs classification with kraken!"
}

gene_classification_with_diamond_megan(){
   echo "Running contigs classification with diamond"

   mkdir -p ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/diamond_output
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      ${TOOLS_FOLDER}/diamond/diamond blastp \
      --threads 24 \
      --query ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/prokka_output/$(basename ${i%.1*}).faa \
      --db ${REFERENCE_FOLDER}/reference_database/nr \
      --daa ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/diamond_output/$(basename ${i%.1*}.daa)


      $HOME/megan/tools/daa2rma \
      --in ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/diamond_output/$(basename ${i%.1*}.daa) \
      --acc2taxa ${REFERENCE_FOLDER}/reference_database/megan_ref/prot_acc2tax-Nov2018X1.abin \
      --acc2interpro2go ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2interpro-June2018X.abin \
      --acc2seed  ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2seed-May2015XX.abin \
      --acc2eggnog ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2eggnog-Oct2016X.abin \
      -fwa true \
      --out ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/diamond_output/$(basename ${i%.1*}.rma)
   done
   
   echo "DONE running contigs classification with diamond!"
}

comparative_functional_annotation_prokka(){
   echo "Running gene functional annotation using prokka"

   #Run prokka
   mkdir -p ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/prokka_output

   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      #ln -s ~/mg-workshop/results/assembly/$SAMPLE/${SAMPLE}_${kmer}/contigs.fa .
      ${TOOLS_FOLDER}/prokka/bin/prokka ${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs_1000_filtered.fasta \
      --outdir ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/prokka_output \
      --locustag 'PROKKA' \
      --kingdom 'Bacteria' \
      --prefix $(basename ${i%.1*}) \
      --norrna \
      --notrna \
      --metagenome \
      --cpus 16 \
      --addgenes \
      --centre X \
      --force

      grep \
         "eC_number=" ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/prokka_output/$(basename ${i%.1*}).gff | cut -f9 | cut -f1,3 -d ';' | sed "s/ID=//g" | sed "s/;eC_number=/\t/g" > ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/prokka_output/$(basename ${i%.1*}).ec
   done

   echo "DONE running gene functional annotation using prokka!"

}

predict_metabolic_pathway_minpath()
{
   echo "Running predict metabolic pathway"

   #awk -F "\t" '{print $1 "\t" $3}' 121832.assembled.EC > for_MinPath.ec
   #sed -i "s/EC://g" for_MinPath.ec
   #MinPath1.2.py -any for_MinPath.ec -map /home/rprops/MinPath/data/ec2path -report report.metacyc.minpath

   mkdir -p ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/minpath_output

   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      ${TOOLS_FOLDER}/MinPath/MinPath1.4.py \
         -any  ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/prokka_output/$(basename ${i%.1*}).ec\
         -map ${TOOLS_FOLDER}/MinPath/data/ec2path \
         -report ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/functional_classification/minpath_output/$(basename ${i%.1*}).report.metacyc.minpath
   done

   echo "DONE running predict metabolic pathway!"
}


#install krken2 and kraken2 database
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
      #export PATH=`pwd`/kraken2/:$PATH
      #export PATH=`pwd`/kraken2-build/:$PATH
      #export PATH=`pwd`/kraken2-inspect/:$PATH
   fi
   echo "DONE checking and installing kraken!"
}

download_krakendb(){
   echo "Checking and downloading kraken database"
   if [ -d "${REFERENCE_FOLDER}/reference_database/kraken2DB" ]; then
      echo "kraken database already installed"
   else
      mkdir -p ${REFERENCE_FOLDER}/reference_database/
      ln -s ${LINKPATH_DB}/reference_database/kraken2DB ${REFERENCE_FOLDER}/reference_database/kraken2DB
   fi
   echo "DONE checking and downloading kraken database!"
}

download_megandb(){
   echo "Checking and downloading megan database"
   if [ -d "${REFERENCE_FOLDER}/reference_database/megan_ref" ]; then
      echo "megan database already installed"
   else
      mkdir -p ${REFERENCE_FOLDER}/reference_database
      ln -s ${LINKPATH_DB}/reference_database/megan_ref ${REFERENCE_FOLDER}/reference_database/megan_ref
   fi
   echo "DONE checking and downloading megan database!"
}

install_diamond(){
   echo "Checking and installing diamond"
   if [ -d "${TOOLS_FOLDER}/diamond" ]; then
      echo "diamond already installed"
   else
      #Installing diamond 
      echo "Installing diamond"
      cd ${TOOLS_FOLDER}
      mkdir -p  ${TOOLS_FOLDER}/diamond/
      cd diamond
      wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz
      tar xzf diamond-linux64.tar.gz
   fi
   echo "DONE checking and installing diamond!"
}


download_nrdb(){
   echo "Checking and downloading NR database"
   if [ -f "${REFERENCE_FOLDER}/reference_database/nr.dmnd" ]; then
      echo "NR database already installed"
   else
      mkdir -p ${REFERENCE_FOLDER}/reference_database
      ln -s ${LINKPATH_DB}/reference_database/nr.dmnd ${REFERENCE_FOLDER}/reference_database/nr.dmnd
   fi
   echo "DONE checking and downloading NR database!"
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

#install megan and download dependency files
install_megan(){
   echo "Checking and installing megan"
   if [ -d "${TOOLS_FOLDER}/megan" ]; then
      echo "megan already installed"
   else
      #Installing megan 
      echo "Installing megan"
      cd ${TOOLS_FOLDER}
      mkdir -p ${TOOLS_FOLDER}/megan 
      cd megan
      wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/MEGAN_Community_unix_6_16_4.sh
      sh MEGAN_Community_unix_6_16_4.sh
      cd ..
   fi
   echo "DONE checking and installing megan!"
}

install_minpath(){
   echo "Checking and installing minpath"
   if [ -d "${TOOLS_FOLDER}/MinPath" ]; then
      echo "minpath already installed"
   else
      #Installing minpath 
      echo "Installing minpath"
      cd ${TOOLS_FOLDER}
      wget "http://omics.informatics.indiana.edu/mg/get.php?justdoit=yes&software=minpath1.4.tar.gz" -O minpath1.4.tar.gz
      tar zxvf minpath1.4.tar.gz
      export MinPath=${TOOLS_FOLDER}/MinPath
      rm minpath1.4.tar.gz
   fi
   export MinPath=${TOOLS_FOLDER}/MinPath
   echo "DONE checking and installing minpath!"
}

install_samtools(){
   echo "Checking and installing samtools"

   if [ -d "${TOOLS_FOLDER}/samtools-1.9" ]; then
      echo "samtools already installed"
   else
      echo "Installing samtools"
      cd $TOOLS_FOLDER
      #install samtools  and dependencies
      git clone git://github.com/samtools/htslib.git
      cd htslib
      make
      #git clone git://github.com/samtools/bcftools.git
      cd ..
      wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
      tar -xvvf samtools-1.9.tar.bz2
      cd samtools-1.9   # and similarly for bcftools and htslib
      make
      export PATH=${TOOLS_FOLDER}/samtools-1.9:$PATH
      cd ..
      rm samtools-1.9.tar.bz2
   fi

   echo "DONE checking and installing samtools!"
}


comparative_analysis_main
