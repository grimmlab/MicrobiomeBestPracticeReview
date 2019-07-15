#!/usr/bin/evn bash

### Comparative Analysis


comparative_analysis_main(){
   check_and_install
   contigs_classification_with_kraken
   contigs_classification_with_diamond_megan
   #gene_based_classification -> EXPERIMENT! UNCOMMENT AT YOUR OWN RISK
   comparative_functional_annotation_prokka
   predict_metabolic_pathway_minpath
}

# Check and install missing packages for the QC pipeline
check_and_install(){
   echo "Checking and installing packages for Comparative Analysis"
   install_diamond
   download_nrdb
   install_kraken
   download_krakendb
   install_prodigal
   install_prokka
   install_megan
   download_megandb
   install_minpath
   install_support_scripts
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

#contigs classification and to visualize it in megan
contigs_classification_with_diamond_megan(){
   echo "Running contigs classification with diamond"

   mkdir -p ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/diamond_output
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      ${TOOLS_FOLDER}/diamond/diamond blastx \
      --threads 24 \
      --query ${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs_1000_filtered.fasta \
      --db ${REFERENCE_FOLDER}/reference_database/nr \
      --daa ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/diamond_output/$(basename ${i%.1*}.daa)

      $HOME/megan/tools/daa-meganizer \
      --in ${ANALYSIS_FOLDER}/contigs_taxonomic_classification/diamond_output/$(basename ${i%.1*}.daa) \
      --acc2taxa ${REFERENCE_FOLDER}/reference_database/megan_ref/prot_acc2tax-Nov2018X1.abin \
      --acc2kegg ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2kegg-Dec2017X1-ue.abin \
      --acc2interpro2go ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2interpro-June2018X.abin \
      --acc2eggnog ${REFERENCE_FOLDER}/reference_database/megan_ref/acc2eggnog-Oct2016X.abin \
      -lg \
      -fwa true
   done
   
   echo "DONE running contigs classification with diamond!"
}

#run_diamond
#https://github.com/rprops/MetaG_analysis_workflow/wiki/11.-Genome-annotation
gene_based_classification(){

   echo "Running gene based taxonomic classification"

   #filterout shortreads
   # Search for genes on contigs
   mkdir -p ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/gene_prediction
   mkdir -p ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      #echo "Running prodigal"
      #${TOOLS_FOLDER}/Prodigal-2.6.1/prodigal \
      #-p meta \
      #-a ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/gene_prediction/$(basename ${i%.1*}.genes.faa) \
      #-d ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/gene_prediction/$(basename ${bin%.1*}.genes.fna) \
      #-f gff \
      #-o ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/gene_prediction/$(basename ${i%.1*}.genes.gff) \
      #-i ${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs_1000_filtered.fasta


      #echo "Running diamond blastp"
      # Blastp to NR database
      #${TOOLS_FOLDER}/diamond/diamond blastp \
      #-p 24 \
      #-d ${REFERENCE_FOLDER}/reference_database/nr \
      #-q ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/gene_prediction/$(basename ${i%.1*}.genes.faa) \
      #-a ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}.daa) > \
      #${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/d.out

      #echo "Running diamond view"
      #${TOOLS_FOLDER}/diamond/diamond view \
      #-a ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}.daa) \
      #-o ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*})_nr.m8


      # Start classification
      #echo "Running Lengths"
      #python ${BIN_FOLDER}/Lengths.py \
      #-i ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/gene_prediction/$(basename ${i%.1*}.genes.faa) \
      #> ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}.length)

      #echo "Modifying m8 file for compatibility"
      #Rscript ${BIN_FOLDER}/modify_m8.R \
      #   -i ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}_nr.m8) \
      #   -o ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}_mod_nr.m8)

      #echo "Running ClassifyContigNR"
      #python ${BIN_FOLDER}/ClassifyContigNR.py \
      #${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}_mod_nr.m8) \
      #${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}.length) \
      #-o ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}_nr) \
      #-l ${REFERENCE_FOLDER}/all_taxa_lineage_notnone.tsv \
      #-g ${REFERENCE_FOLDER}/gi_taxid_prot.dmp
      

      # Extracting specific classification level
      echo "Running Filter"
      perl ${BIN_FOLDER}/Filter.pl \
      4 < \
      ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}_nr_contigs.csv) > \
      ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}_classification_contigs_4.csv)

      # Make tsv
      cat ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}_classification_contigs_4.csv) | sed 's/,/\t/g' >   ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}_classification_contigs_4.tsv)

      # Get contig list
      grep ">" ${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs_1000_filtered.fasta | sed "s/>//g" > ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}.txt)

      # Format to annotation file for vizbin
      echo "Running diamond2vizbin"
      #Rscript ${BIN_FOLDER}/diamond2vizbin.R \
      #   ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}_classification_contigs_4.tsv) \
      #   ${ANALYSIS_FOLDER}/Gene_based_analysis_onContigs/taxonomic_classification/assignTaxa/$(basename ${i%.1*}.txt)
   done

   echo "DONE running gene based taxonomic classification!"
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

#install krken2 and kraken2 database
install_prodigal(){
   echo "Checking and installing prodigal"
   if [ -d "${TOOLS_FOLDER}/Prodigal-2.6.1" ]; then
      echo "prodigal already installed"
   else
      #Installing prodigal 
      echo "Installing prodigal"
      cd ${TOOLS_FOLDER}
      wget https://github.com/hyattpd/Prodigal/archive/v2.6.1.tar.gz -O prodigal.tar.gz
      tar xzf prodigal.tar.gz
      cd Prodigal-2.6.1 
      make
      cd ..
      rm Prodigal.tar.gz
   fi
   export PATH=${TOOLS_FOLDER}/Prodigal-2.6.1/:$PATH
   echo "DONE checking and installing prodigal!"
}

install_support_scripts(){

   echo "Checking and installing support scripts"
   if [ -d "${BIN_FOLDER}/DESMAN-master"  -a -d "${BIN_FOLDER}/MetaG_analysis_workflow-master" ]; then
      echo "support scripts already installed"
   else
      #Installing support scripts 
      cd ${BIN_FOLDER}
      echo "Installing support scripts"
      wget https://github.com/rprops/DESMAN/archive/master.zip -O DESMAN.zip
      unzip DESMAN.zip
      cp DESMAN-master/scripts/Lengths.py .
      cp DESMAN-master/scripts/ClassifyContigNR.py .
      cp DESMAN-master/scripts/Filter.pl .
      rm DESMAN.zip

      wget https://github.com/rprops/MetaG_analysis_workflow/archive/master.zip -O Metagenomics_WorkFlow.zip 
      unzip Metagenomics_WorkFlow.zip
      cp MetaG_analysis_workflow-master/diamond2vizbin.R .
      rm Metagenomics_WorkFlow.zip

      cd ${REFERENCE_FOLDER}
      wget https://desmandatabases.s3.climb.ac.uk/all_taxa_lineage_notnone.tsv
      wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
      wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz
      gunzip prot.accession2taxid.gz
      gunzip gi_taxid_prot.dmp.gz
   fi

   echo "DONE checking and installing support scripts!"
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

comparative_analysis_main
