#!/usr/bin/evn bash

### coverage and bining 


coverage_and_bining_main(){
   check_and_install
   create_coveragefile
   binning_maxbin
   binning_metabat
}

# check and install missing packages for the qc pipeline
check_and_install(){
   echo "checking and installing packages for coverage and bining"
   #install_metabat # needs boost libraries (may need to be installed separately)
   install_bowtie
   install_fraggenescan
   install_hmmer
   install_idba
   install_maxbin
   install_samtools
   install_bbmap
   echo "DONE checking and installing packages for coverage and bining!"
}


create_coveragefile(){
   echo "Running Create coverage file"

   mkdir -p ${ANALYSIS_FOLDER}/coverage
   mkdir -p ${ANALYSIS_FOLDER}/coverage/bam
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      echo "Running bbwrap in coverage file"
      ${TOOLS_FOLDER}/bbmap/bbwrap.sh \
      ref=${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs_1000_filtered.fasta \
      nodisk \
      in=${i},${i%1*}u.final.clean.fq \
      in2=${i%1*}2.unmerged.final.clean.fq,NULL \
      t=24 \
      kfilter=22 \
      subfilter=15 \
      maxindel=80 \
      out=${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}sam) \
      covstats=${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}coverage)

      #converting sam file to bam
      ${TOOLS_FOLDER}/samtools-1.9/samtools view \
      -S \
      -b \
      -u \
      ${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}sam) \
      > ${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}bam)

      #sorting the bam file
      ${TOOLS_FOLDER}/samtools-1.9/samtools sort \
      -m 16G \
      -@ 3 \
      ${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}bam) \
      -o ${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}sorted.bam)

      #indexing the sorted bam file
      ${TOOLS_FOLDER}/samtools-1.9/samtools index \
      ${ANALYSIS_FOLDER}/coverage/bam/$(basename ${i%1*}sorted.bam)

   done

   jgi_summarize_bam_contig_depths \
   --outputDepth ${ANALYSIS_FOLDER}/coverage/depth.txt \
   --pairedContigs ${ANALYSIS_FOLDER}/coverage/paired.txt \
   --minContigLength 1000 \
   --minContigDepth 2 \
   ${ANALYSIS_FOLDER}/coverage/bam/*sorted.bam

   tail -n+2 ${ANALYSIS_FOLDER}/coverage/depth.txt | cut -f 1,3 > ${ANALYSIS_FOLDER}/coverage/maxbin.cov

   echo "DONE running create coverage file!"

}

binning_maxbin(){
   echo "Running Binning MaxBin"

   #chose the contigs
   #samplename=${i%1*}
   mkdir -p ${ANALYSIS_FOLDER}/binning/
   mkdir -p ${ANALYSIS_FOLDER}/binning/maxbin

   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      ${TOOLS_FOLDER}/MaxBin-2.2.6/run_MaxBin.pl \
      -thread 16 \
      -contig ${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs_1000_filtered.fasta \
      -abund ${ANALYSIS_FOLDER}/coverage/maxbin.cov \
      -out ${ANALYSIS_FOLDER}/binning/maxbin/$(basename ${i%.1*})_maxbin
   done
   echo "DONE running Binning MaxBin!"
}

binning_metabat(){
   echo "Running Binning Metabat"

   mkdir -p ${ANALYSIS_FOLDER}/binning/metabat
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
     # ${TOOLS_FOLDER}/metabat/metabat \
      metabat \
      -i ${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs.fasta \
      -a ${ANALYSIS_FOLDER}/coverage/depth.txt \
      -o ${ANALYSIS_FOLDER}/binning/metabat/$(basename ${i%.1*})_metabat \
      -t 16 \
      -m 1500 \
      -v \
      --unbinned
   done

   #In order to check which contigs were grouped together into separate bins by MetaBat example
   #grep ">" E01452_L001_to_L004_metabat.1.fa | sed 's/>//g' > E01452_L001_to_L004_metabat.1.contigNames

   echo "DONE running Binning Metabat!"

}


install_bbmap(){
   echo "Checking and installing bbmap"
   if [ -d "${TOOLS_FOLDER}/bbmap" ]; then
      echo "bbmap already installed"
   else
      echo "Installing bbmap"
      cd $TOOLS_FOLDER
      wget http://downloads.sourceforge.net/project/bbmap/BBMap_38.44.tar.gz
      tar zvxf BBMap_38.44.tar.gz
      rm BBMap_38.44.tar.gz
   fi
   echo "DONE checking and installing bbmap!"
}

install_metabat(){
   echo "Checking and installing metabat"
   if [ -d "${TOOLS_FOLDER}/metabat" ]; then
      echo "metabat already installed"
   else
      echo "Installing metabat"
      cd $TOOLS_FOLDER

      git clone https://bitbucket.org/berkeleylab/metabat.git
      mkdir build 
      cd build 
      cmake .. 
      make
      make install
      cd ..
   fi
   echo "DONE checking and installing metabat!"
}

install_maxbin(){
   echo "Checking and installing maxbin"
   if [ -d "${TOOLS_FOLDER}/MaxBin-2.2.6" ]; then
      echo "maxbin already installed"
   else
      echo "Installing maxbin"
      cd $TOOLS_FOLDER
      #MaxBin
      wget https://sourceforge.net/projects/maxbin/files/MaxBin-2.2.6.tar.gz
      tar zxvf MaxBin-2.2.6.tar.gz
      cd MaxBin-2.2.6
      # Settings file update
      mv setting setting.org
      echo "[FragGeneScan]" ${TOOLS_FOLDER}/FragGeneScan1.31 >> setting
      echo "[Bowtie2]" ${TOOLS_FOLDER}/bowtie2-2.3.3.1-linux-x86_64 >> setting
      echo "[HMMER3]" ${TOOLS_FOLDER}/hmmer-3.2.1 >> setting
      echo "[IDBA_UD]" ${TOOLS_FOLDER}/idba-master/bin >> setting
      
      cd src
      make
      cd ../../
      rm MaxBin-2.2.6.tar.gz
   fi
   echo "DONE checking and installing maxbin!"
}

install_bowtie(){
   echo "Checking and installing bowtie"
   if [ -d "${TOOLS_FOLDER}/bowtie2-2.3.3.1-linux-x86_64" ]; then
      echo "bowtie already installed"
   else
      echo "Installing bowtie"
      cd $TOOLS_FOLDER
      wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3.1/bowtie2-2.3.3.1-linux-x86_64.zip/download
      unzip download
      rm download
      export PATH=$PATH"${TOOLS_FOLDER}/bowtie2-2.3.3.1-linux-x86_64:"
   fi
   echo "DONE checking and installing bowtie!"
}


install_fraggenescan(){
   echo "Checking and installing FragGeneScan"

   if [ -d "${TOOLS_FOLDER}/FragGeneScan1.31" ]; then
      echo "FragGeneScan already installed"
   else
      echo "Installing FragGeneScan"
      cd $TOOLS_FOLDER
      wget https://sourceforge.net/projects/fraggenescan/files/latest/download/FragGeneScan1.31.tar
      tar zxvf FragGeneScan1.31.tar
      cd FragGeneScan1.31/
      make
      make clean
      make fgs
      cd ..
      rm FragGeneScan1.31.tar
   fi

   echo "DONE checking and installing FragGeneScan!"
}

install_hmmer(){
   echo "Checking and installing hmmer"

   if [ -d "${TOOLS_FOLDER}/hmmer-3.2.1" ]; then
      echo "hmmer already installed"
   else
      echo "Installing hmmer"
      cd $TOOLS_FOLDER
      wget http://eddylab.org/software/hmmer/hmmer.tar.gz
      tar -xvvf hmmer.tar.gz
      cd hmmer-3.2.1/
      ./configure --prefix /data1/Active_Projects/Metagenomic_QC/tools/hmmer-3.2.1/
      make
      make install
      rm hmmer.tar.gz
   fi

   export PATH=${TOOLS_FOLDER}/hmmer-3.2.1/:$PATH

   echo "DONE checking and installing hmmer!"

}

#Installing idba_ud
install_idba(){
   echo "Checking and installing idba"
   if [ -d "${TOOLS_FOLDER}/idba-master" ]; then
      echo "idba already installed"
   else
      echo "Installing idba"
      cd $TOOLS_FOLDER
      wget https://github.com/loneknightpy/idba/archive/master.zip
      unzip master.zip
      cd idba-master/
      ./build.sh
      ./configure
      make
      cd ..
      rm master.zip
   fi
   echo "DONE checking and installing idba!"
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


coverage_and_bining_main
