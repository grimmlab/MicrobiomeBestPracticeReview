#!/usr/bin/env bash
# Use this script only once initially and save the path

prepare_databases_main(){
   create_folders
   #prepare_hostgenome_db
   #prepare_megandb
   #prepare_nrdb
   #prepare_humanndb
   #prepare_krakendb
   #prepare_metaphlan
   cleanup_tmp
   echo "Use this path in the workflow scripts"
   echo "LINKPATH_DB="$REFERENCE_FOLDER 
}



create_folders(){

   echo "Creating sub-folders..."

   # Sub-folders in the root folder
   for FOLDER in reference 
   do
      mkdir -p $PWD/$FOLDER
   done
   export REFERENCE_FOLDER=$PWD/$FOLDER
   mkdir -p ${REFERENCE_FOLDER}/tmp
   export TMP=${REFERENCE_FOLDER}/tmp
   echo "DONE creating sub-folders!"
}

#### Downloading and indexing the Host genome (human and mouse) using prinseq for cleaning and BBMAP for indexing the genome
#Creating human reference database
prepare_hostgenome_db(){
   echo "Downloading and indexing genome"

   if [ -d "${REFERENCE_FOLDER}/human" ]; then
      echo "human genome already exists"
   else
      echo "Downloading human genome"
      mkdir -p ${REFERENCE_FOLDER}/human
      cd ${REFERENCE_FOLDER}/human
      for i in {1..22} X Y MT
      do
         wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38.p12_chr$i.fa.gz
      done

      for i in {1..22} X Y MT
      do
         gzip -dvc hs_ref_GRCh38.p12_chr$i.fa.gz >>${REFERENCE_FOLDER}/human/hsref_GRCh38_p12.fa
      done

      cd ${TMP}
      wget https://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz
      tar xvf prinseq-lite-0.20.4.tar.gz
      rm prinseq-lite-0.20.4.tar.gz

      cd ${TMP}
      wget http://downloads.sourceforge.net/project/bbmap/BBMap_38.44.tar.gz
      tar zvxf BBMap_38.44.tar.gz
      rm BBMap_38.44.tar.gz

      #Clean
      echo "Cleaning human genome"
      perl ${TMP}/prinseq-lite-0.20.4/prinseq-lite.pl \
      -log \
      -verbose \
      -fasta ${REFERENCE_FOLDER}/human/hsref_GRCh38_p12.fa \
      -min_len 200 \
      -ns_max_p 10 \
      -derep 12345 \
      -out_good ${REFERENCE_FOLDER}/human/hsref_GRCh38_p12_clean \
      -seq_id hsref_GRCh38_p12_ \
      -rm_header \
      -out_bad null

      rm ${REFERENCE_FOLDER}/human/hsref_GRCh38_p12.fa
      rm ${REFERENCE_FOLDER}/human/*fa.gz

      #indexing the genome
      echo "Indexing human genome"
      ${TMP}/bbmap/bbmap.sh \
      ref=${REFERENCE_FOLDER}/human/hsref_GRCh38_p12_clean.fasta \
      path=${REFERENCE_FOLDER}/human
   fi


   if [ -d "${REFERENCE_FOLDER}/mouse" ]; then
      echo "mouse genome already exists"
   else
      #Creating mouse reference database
      echo "Downloading mouse genome"
      mkdir -p $REFERENCE_FOLDER/mouse
      cd $REFERENCE_FOLDER/mouse
      for i in {1..19} X Y MT
      do
         wget ftp://ftp.ncbi.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/mm_ref_GRCm38.p4_chr$i.fa.gz
      done

      for i in {1..19} X Y MT
      do
         gzip -dvc mm_ref_GRCm38.p4_chr$i.fa.gz >>$REFERENCE_FOLDER/mouse/mmref_GRCm38_p4.fa
      done

      #Clean
      echo "Cleaning mouse genome"
      perl ${TMP}/prinseq-lite-0.20.4/prinseq-lite.pl \
      -log \
      -verbose \
      -fasta ${REFERENCE_FOLDER}/mouse/mmref_GRCm38_p4.fa \
      -min_len 200 \
      -ns_max_p 10 \
      -derep 12345 \
      -out_good ${REFERENCE_FOLDER}/mouse/mmref_GRCm38_p4_clean \
      -seq_id mmref_GRCm38_p4_ \
      -rm_header \
      -out_bad null


      rm ${REFERENCE_FOLDER}/mouse/mmref_GRCm38_p4.fa
      rm ${REFERENCE_FOLDER}/mouse/*fa.gz

      echo "Indexing mouse genome"
      ${TMP}/bbmap/bbmap.sh \
      ref=${REFERENCE_FOLDER}/mouse/mmref_GRCm38_p4_clean.fasta \
      path=${REFERENCE_FOLDER}/mouse
   fi

   echo "DONE downloading and indexing genome!"
}

prepare_megandb(){
   echo "Checking and downloading megan database"
   if [ -d "${REFERENCE_FOLDER}/reference_database/megan_ref" ]; then
      echo "megan database already installed"
   else
      echo "Downloading megan database (Takes a long time, be patient)"
      mkdir -p  ${REFERENCE_FOLDER}/reference_database/megan_ref
      cd ${REFERENCE_FOLDER}/reference_database/megan_ref
      wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/prot_acc2tax-Nov2018X1.abin.zip
      unzip prot_acc2tax-Nov2018X1.abin.zip
      wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/acc2kegg-Dec2017X1-ue.abin.zip
      unzip acc2kegg-Dec2017X1-ue.abin.zip
      wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/acc2interpro-June2018X.abin.zip
      unzip acc2interpro-June2018X.abin.zip
      wget http://ab.inf.uni-tuebingen.de/data/software/megan6/download/acc2eggnog-Oct2016X.abin.zip
      unzip acc2eggnog-Oct2016X.abin.zip
   fi
   echo "DONE checking and downloading megan database!"
}

prepare_nrdb(){
   echo "Checking and downloading NR database"
   if [ -f "${REFERENCE_FOLDER}/reference_database/nr.dmnd" ]; then
      echo "NR database already installed"
   else
      echo "Downloading NR database (Takes a long time, be patient)"
      wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
      mkdir -p  ${TMP}/diamond/
      cd ${TMP}/diamond
      wget http://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz
      tar xzf diamond-linux64.tar.gz
      ./diamond makedb --in nr.gz --db ${REFERENCE_FOLDER}/reference_database/nr
   fi
   echo "DONE checking and downloading NR database!"
}

prepare_humanndb(){
   echo "Checking and downloading humann2 database"
   if [ -d "${REFERENCE_FOLDER}/reference_database/humann2_databases" ]; then
      echo "humann2 database already installed"
   else
      cd ${TMP}
      hg clone https://bitbucket.org/biobakery/humann2
      cd humann2
      python3 setup.py install --user
      echo "Downloading humann2 database (Takes a long time, be patient)"
      mkdir -p ${REFERENCE_FOLDER}/reference_database/humann2_databases
      humann2_databases --download chocophlan full ${REFERENCE_FOLDER}/reference_database/humann2_databases
      #To download the full UniRef90 database (11.0GB, recommended):
      humann2_databases --download uniref uniref90_diamond ${REFERENCE_FOLDER}/reference_database/humann2_databases
      #To download the EC-filtered UniRef90 database (0.8GB):
      humann2_databases --download uniref uniref90_ec_filtered_diamond ${REFERENCE_FOLDER}/reference_database/humann2_databases
      #To download the full UniRef50 database (4.6GB):
      humann2_databases --download uniref uniref50_diamond ${REFERENCE_FOLDER}/reference_database/humann2_databases
      #To download the EC-filtered UniRef50 database (0.2GB):
      humann2_databases --download uniref uniref50_ec_filtered_diamond ${REFERENCE_FOLDER}/reference_database/humann2_databases
      humann2_config --update database_folders nucleotide ${REFERENCE_FOLDER}/reference_database/humann2_databases/chocophlan/
      humann2_config --update database_folders protein ${REFERENCE_FOLDER}/reference_database/humann2_databases/uniref/
   fi
   echo "DONE checking and downloading humann2 database!"
}

prepare_krakendb(){
   echo "Checking and downloading kraken database"
   if [ -d "${REFERENCE_FOLDER}/reference_database/kraken2DB" ]; then
      echo "kraken database already installed"
   else
      cd ${TMP}
      wget http://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz
      tar xzf v2.0.8-beta.tar.gz
      cd kraken2-2.0.8-beta
      sh install_kraken2.sh ${TMP}/kraken2-2.0.8-beta
      cd ..
      rm v2.0.8-beta.tar.gz

      # building kraken database
      ${TMP}/kraken2-2.0.8-beta/kraken2-build \
         --standard \
         --db ${REFERENCE_FOLDER}/reference_database/kraken2DB \
         --threads 16 \
         --use-ftp
   fi
   echo "DONE checking and downloading kraken database!"
}

prepare_metaphlan(){
   echo "Checking and downloading metaphlan database"
   if [ -d "${REFERENCE_FOLDER}/reference_database/metaphlan2" ]; then
      echo "metaphlan database already installed"
   else
      cd $TMP
      wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3.1/bowtie2-2.3.3.1-linux-x86_64.zip/download
      unzip download
      rm download
      export PATH=$PATH"${TMP}/bowtie2-2.3.3.1-linux-x86_64:"

      mkdir -p ${REFERENCE_FOLDER}/reference_database/metaphlan2
      cd ${REFERENCE_FOLDER}/reference_database/metaphlan2
      wget https://bitbucket.org/biobakery/metaphlan2/downloads/mpa_v20_m200.tar
      tar -xvvf mpa_v20_m200.tar
      bunzip2 mpa_v20_m200.fna.bz2
      bowtie2-build --threads 20 mpa_v20_m200.fna mpa_v20_m200
      rm mpa_v20_m200.tar
   fi
   echo "DONE checking and downloading metaphlan database!"
}

cleanup_tmp(){

   if [ -d "${TMP}" ]; then
      rm -rf ${TMP}
   fi

}


prepare_databases_main
