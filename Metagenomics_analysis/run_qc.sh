#!/usr/bin/evn bash

### QC Analysis

FASTQC_VERSION=v0.11.8

qc_main(){
   check_and_install
   fastqc_stats
   trimmomatic_sickle_QC
   bbmap_QC1
   bbmap_QC2
   creating_stats
   copy_results
}

# Check and install missing packages for the QC pipeline
check_and_install(){
   echo "Checking and installing packages for Quality Control"
   install_fastqc
   install_trimmomatic
   install_sickle
   install_bbmap
   install_prinseq
   downloads_indexing
   echo "DONE checking and installing packages for Quality Control!"
}

#### 1. Generating comprehensive report and stats of data quality using fastqc and BBMAP
fastqc_stats(){
   echo "Running fastqc stats"

   echo "Running fastqc"
   #### Run fastqc on our data (non interactively)
   cd ${ANALYSIS_FOLDER}
   mkdir -p ${ANALYSIS_FOLDER}/QC/fastqc

   find ${RAWDATA_FOLDER} \
   -name "*.fastq.gz" | \
   xargs -n 1 \
   ${TOOLS_FOLDER}/FastQC/fastqc \
   -o ${ANALYSIS_FOLDER}/QC/fastqc
   #DONE


   echo "Creating rawdata stats"
   mkdir -p ${ANALYSIS_FOLDER}/QC/bbmap
   rawdatalist=$(ls -d $RAWDATA_FOLDER/*_1.fastq.gz | awk '{print $NF}')
   for s in $rawdatalist
   do
      $TOOLS_FOLDER/bbmap/reformat.sh \
      threads=16 \
      in=${s} \
      in2=${s%_1*}_2.fastq.gz \
      2>&1 >/dev/null | awk '{print "RAWREADS "$0}' | tee -a $ANALYSIS_FOLDER/QC/bbmap/$(basename ${s%_1*}stats.txt)
   done

   echo "DONE Running fastqc stats!"
}

#### 2. Quality Control Trimming using trimmomatic and Sickle
trimmomatic_sickle_QC(){

   echo "Running trimmomatic and sickle"
   #Trimming low quality, short length reads, adapters
   mkdir -p $ANALYSIS_FOLDER/QC/trimmomatic
   mkdir -p $ANALYSIS_FOLDER/QC/sickle
   rawdatalist=$(ls -d $RAWDATA_FOLDER/*_1.fastq.gz | awk '{print $NF}')

   for s in $rawdatalist
   do
      fname=$(basename $s | sed -e "s/_1.fastq.gz//")
      #Running trimmomatic
      echo "Running trimmomatic"
      java -jar $TOOLS_FOLDER/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
      ${s%_1*}_1.fastq.gz ${s%_1*}_2.fastq.gz \
      -threads 16 \
      -trimlog ${ANALYSIS_FOLDER}/QC/trimmomatic/${fname}.trimlog.txt \
      -phred33 \
      ${ANALYSIS_FOLDER}/QC/trimmomatic/${fname}.1.trimmoclean.fq.gz \
      ${ANALYSIS_FOLDER}/QC/trimmomatic/${fname}.1.u.trimmoclean.fq.gz \
      ${ANALYSIS_FOLDER}/QC/trimmomatic/${fname}.2.trimmoclean.fq.gz \
      ${ANALYSIS_FOLDER}/QC/trimmomatic/${fname}.2.u.trimmoclean.fq.gz \
      ILLUMINACLIP:$TOOLS_FOLDER/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:1:50:30 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60

      #Running sickle
      echo "Running sickle pe"
      $TOOLS_FOLDER/sickle-1.33/sickle pe \
      -n \
      -f ${ANALYSIS_FOLDER}/QC/trimmomatic/${fname}.1.trimmoclean.fq.gz \
      -r ${ANALYSIS_FOLDER}/QC/trimmomatic/${fname}.2.trimmoclean.fq.gz \
      -t sanger \
      -o ${ANALYSIS_FOLDER}/QC/sickle/${fname}.1.trimmoclean.sickleclean.fq \
      -p ${ANALYSIS_FOLDER}/QC/sickle/${fname}.2.trimmoclean.sickleclean.fq \
      -s ${ANALYSIS_FOLDER}/QC/sickle/${fname}.u.trimmoclean.sickleclean.fq \
      -q 20 \
      -l 60

      #QC for unpaired reads
      echo "Running sickle se"
      $TOOLS_FOLDER/sickle-1.33/sickle se \
      -n \
      -f ${ANALYSIS_FOLDER}/QC/trimmomatic/${fname}.1.u.trimmoclean.fq.gz \
      -o ${ANALYSIS_FOLDER}/QC/sickle/${fname}.1.u.trimmoclean.sickleclean.fq \
      -t sanger \
      -q 20 \
      -l 60

      echo "Running sickle se"
      $TOOLS_FOLDER/sickle-1.33/sickle se \
      -n \
      -f ${ANALYSIS_FOLDER}/QC/trimmomatic/${fname}.2.u.trimmoclean.fq.gz \
      -o ${ANALYSIS_FOLDER}/QC/sickle/${fname}.2.u.trimmoclean.sickleclean.fq \
      -t sanger \
      -q 20 \
      -l 60

      #Combining all unpaired files
      cat ${ANALYSIS_FOLDER}/QC/sickle/${fname}.1.u.trimmoclean.sickleclean.fq \
      ${ANALYSIS_FOLDER}/QC/sickle/${fname}.2.u.trimmoclean.sickleclean.fq \
      ${ANALYSIS_FOLDER}/QC/sickle/${fname}.u.trimmoclean.sickleclean.fq \
      > ${ANALYSIS_FOLDER}/QC/sickle/${fname}.unpaired.trimmoclean.sickleclean.fq

      rm ${ANALYSIS_FOLDER}/QC/sickle/${fname}.1.u.trimmoclean.sickleclean.fq ${ANALYSIS_FOLDER}/QC/sickle/${fname}.2.u.trimmoclean.sickleclean.fq ${ANALYSIS_FOLDER}/QC/sickle/${fname}.u.trimmoclean.sickleclean.fq
   done

   echo "DONE running trimmomatic and sickle!"
   #DONE
}

#### 3. Quality control removing phix adapters and sequencing artifacts using BBMAP
bbmap_QC1(){
   echo "Running bbmap"

   # paired data
   sicklelist=$(ls -d ${ANALYSIS_FOLDER}/QC/sickle/*1.trimmoclean.sickleclean.fq | awk '{print $NF}')
   echo "Running bbduk for paired"
   for s in $sicklelist
   do
      sname=$(basename ${s} | sed -e "s/1.trimmoclean.sickleclean.fq//")
      $TOOLS_FOLDER/bbmap/bbduk.sh \
      threads=8 \
      in=${s} \
      in2=${s%1*}2.trimmoclean.sickleclean.fq \
      k=31 \
      ref=${TOOLS_FOLDER}/bbmap/resources/sequencing_artifacts.fa.gz,${TOOLS_FOLDER}/bbmap/resources/phix_adapters.fa.gz \
      out1=$ANALYSIS_FOLDER/QC/bbmap/${sname}1.trimmoclean.sickleclean.bbdukclean.fq \
      out2=$ANALYSIS_FOLDER/QC/bbmap/${sname}2.trimmoclean.sickleclean.bbdukclean.fq \
      minlength=60
   done

   # unpaired data
   echo "Running bbduk for unpaired"
   sickleunplist=$(ls -d ${ANALYSIS_FOLDER}/QC/sickle/*unpaired.trimmoclean.sickleclean.fq| awk '{print $NF}')
   for s in $sickleunplist
   do
      sname=$(basename ${s} | sed -e "s/unpaired.trimmoclean.sickleclean.fq//")
      $TOOLS_FOLDER/bbmap/bbduk.sh \
      threads=8 \
      in=${s} \
      k=31 \
      ref=${TOOLS_FOLDER}/bbmap/resources/sequencing_artifacts.fa.gz,${TOOLS_FOLDER}/bbmap/resources/phix_adapters.fa.gz \
      out1=$ANALYSIS_FOLDER/QC/bbmap/${sname}unpaired.trimmoclean.sickleclean.bbdukclean.fq \
      minlength=60
   done

   echo "DONE running bbmap!"
}

#### 4. Removing the host contamination and generating the stats of the data using BBMAP
bbmap_QC2(){
   echo "Running bbmap 2"

   bbduklist=$(ls -d ${ANALYSIS_FOLDER}/QC/bbmap/*.1.trimmoclean.sickleclean.bbdukclean.fq | awk '{print $NF}')

   echo "Running bbwrap"
   for s in $bbduklist
   do
      $TOOLS_FOLDER/bbmap/bbwrap.sh \
      threads=16 \
      minid=0.95 \
      maxindel=3 \
      bwr=0.16 \
      bw=12 \
      quickmatch \
      fast \
      minhits=2 \
      qtrim=rl \
      trimq=20 \
      minlength=60 \
      in=${s},${s%1*}unpaired.trimmoclean.sickleclean.bbdukclean.fq \
      in2=${s%1*}2.trimmoclean.sickleclean.bbdukclean.fq,NULL \
      outu1=${s%1*}1.final.clean.fq \
      outu2=${s%1*}2.final.clean.fq \
      outu=${s%1*}u.clean.fq \
      path=${REFERENCE_FOLDER}/mouse/ 2>&1 >/dev/null | \
      awk '{print "HOST CONTAMINATION SEQUENCES "$0}' | \
      tee -a $ANALYSIS_FOLDER/QC/bbmap/$(basename ${s%1*}stats.txt)
   done
   #DONE

   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/bbmap/*.1.final.clean.fq | awk '{print $NF}')
   #in=pe1_1.fq,pe2_1.fq,se.fq in2=pe1_2.fq,pe2_2.fq
   echo "Running bbmerge"
   for s in $finallist
   do
      $TOOLS_FOLDER/bbmap/bbmerge.sh \
      threads=4 \
      in1=${s} \
      in2=${s%1*}2.final.clean.fq \
      out=${s%1*}merged.final.clean.fq \
      outu1=${s%1*}1.unmerged.final.clean.fq \
      outu2=${s%1*}2.unmerged.final.clean.fq \
      mininsert=60 \
      2>&1 >/dev/null | awk '{print "MERGED "$0}' | \
      tee -a $ANALYSIS_FOLDER/QC/bbmap/$(basename ${s%1*}stats.txt)
   done
   #DONE

   #creating interleaved file for diamond
   echo "Running reformat"
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/bbmap/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      ${TOOLS_FOLDER}/bbmap/reformat.sh \
      in1=${i} \
      in2=${i%1*}2.unmerged.final.clean.fq  \
      out=${i%1*}12.interleaved.final.clean.fa
   done

   for s in $finallist
   do
      cat ${s%1*}merged.final.clean.fq  ${s%1*}u.clean.fq  > ${s%1*}u.final.clean.fq
   done
   #DONE

   echo "Running reformat"
   for s in $finallist
   do
      $TOOLS_FOLDER/bbmap/reformat.sh \
      threads=16 \
      in=${s%1*}u.final.clean.fq \
      2>&1 >/dev/null | awk '{print "UNPAIRED "$0}' | tee -a $ANALYSIS_FOLDER/QC/bbmap/$(basename ${s%1*}stats.txt)
   done

   echo "Running reformat"
   for s in $finallist
   do
      $TOOLS_FOLDER/bbmap/reformat.sh \
      threads=16 \
      in1=${s%1*}1.unmerged.final.clean.fq \
      in2=${s%1*}2.unmerged.final.clean.fq  2>&1 >/dev/null | \
      awk '{print "PAIRED "$0}' | \
      tee -a  $ANALYSIS_FOLDER/QC/bbmap/$(basename ${s%1*}stats.txt)
   done

   echo "DONE running bbmap 2!"
}

#### Downloading and indexing the Host genome (human and mouse) using prinseq for cleaning and BBMAP for indexing the genome
#Creating human reference database
downloads_indexing(){
   echo "Downloading and indexing genome"

   if [ -d "${REFERENCE_FOLDER}/human" ]; then
      echo "human genome already exists"
   else
      ln -s $LINKPATH_DB/human ${REFERENCE_FOLDER}/human
   fi

   if [ -d "${REFERENCE_FOLDER}/mouse" ]; then
      echo "mouse genome already exists"
   else
      ln -s $LINKPATH_DB/mouse ${REFERENCE_FOLDER}/mouse
   fi

   echo "DONE downloading and indexing genome!"
}

creating_stats(){

   echo "Creating stats"

   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/bbmap/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for s in $finallist
   do
      grep 'RAWREADS' ${s%1*}stats.txt  | grep 'Input:' | awk '{print "RAWREADS COUNT""\t"$3/2}' | tee ${s%1*}finalstats.txt
      grep 'RAWREADS' ${s%1*}stats.txt  | grep 'Input:' | awk '{print "BASES RAWREADS "$5}' | tee -a ${s%1*}finalstats.txt
      grep 'HOST CONTAMINATION SEQUENCES' ${s%1*}stats.txt | grep "Reads Used:"  | awk '{printf $4" "}' | awk '{print "READS BIO "$1/2 + $2}' | tee -a ${s%1*}finalstats.txt
      egrep ^UNPAIRED ${s%1*}stats.txt  | grep 'Input:' | awk '{print $3}' | awk '{print "READS CLEAN_UNPAIRED "$1}' | tee -a ${s%1*}finalstats.txt
      egrep ^UNPAIRED ${s%1*}stats.txt  | grep 'Input:' | awk '{print "BASES CLEAN_UNPAIRED "$5}' | tee -a ${s%1*}finalstats.txt
      egrep ^PAIRED ${s%1*}stats.txt  | grep 'Input:' | awk '{print $3}' | awk '{print "READS CLEAN_PAIRED "$1}' | tee -a ${s%1*}finalstats.txt
      egrep ^PAIRED ${s%1*}stats.txt  | grep 'Input:' | awk '{print "BASES CLEAN_PAIRED "$5}' | tee -a ${s%1*}finalstats.txt
   done
   echo "DONE creating stats!"
}

copy_results(){
   echo "Copying clean files to the folder"
   mkdir -p ${ANALYSIS_FOLDER}/QC/final_QC_output/
   cp ${ANALYSIS_FOLDER}/QC/bbmap/*.unmerged.final.clean.fq ${ANALYSIS_FOLDER}/QC/final_QC_output/
   cp ${ANALYSIS_FOLDER}/QC/bbmap/*u.final.clean.fq ${ANALYSIS_FOLDER}/QC/final_QC_output/
   cp ${ANALYSIS_FOLDER}/QC/bbmap/*interleaved.final.clean.fa ${ANALYSIS_FOLDER}/QC/final_QC_output/
   cp ${ANALYSIS_FOLDER}/QC/bbmap/*1.final.clean.fq ${ANALYSIS_FOLDER}/QC/final_QC_output/
   cp ${ANALYSIS_FOLDER}/QC/bbmap/*2.final.clean.fq ${ANALYSIS_FOLDER}/QC/final_QC_output/

   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/bbmap/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for s in $finallist
   do
      cp ${s%1*}finalstats.txt ${ANALYSIS_FOLDER}/QC/final_QC_output/
   done

   echo "DONE copying clean files to the folder!"
}

# Check and install missing packages for the QC pipeline
# Install fastqc
install_fastqc(){
   echo "Checking and installing fastqc"
   if [ -d "${TOOLS_FOLDER}/FastQC" ]; then
      echo "fastqc already installed"
   else
      #Installing fastqc
      echo "Installing fastqc"
      cd $TOOLS_FOLDER
      wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_${FASTQC_VERSION}.zip
      unzip fastqc_${FASTQC_VERSION}.zip
      cd $TOOLS_FOLDER/FastQC
      chmod 755 fastqc
      cd ..
      rm fastqc_${FASTQC_VERSION}.zip -r
   fi
   echo "DONE checking and installing fastqc!"
}

# Install trimmomatic
install_trimmomatic(){
   echo "Checking and installing Trimmomatic"
   if [ -d "${TOOLS_FOLDER}/Trimmomatic-0.39" ]; then
      echo "Trimmomatic already installed"
   else
      echo "Installing Trimmomatic"
      #Installing trimmomatic
      cd $TOOLS_FOLDER
      wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
      unzip Trimmomatic-0.39.zip
      rm Trimmomatic-0.39.zip
   fi
   echo "DONE checking and installing trimmomatic!"
}

install_sickle(){
   echo "Checking and installing sickle"
   if [ -d "${TOOLS_FOLDER}/sickle-1.33" ]; then
      echo "sickle already installed"
   else
      echo "Installing sickle"
      #Installing sickle
      cd $TOOLS_FOLDER
      wget https://github.com/najoshi/sickle/archive/v1.33.tar.gz
      tar -xzvf v1.33.tar.gz
      cd sickle-1.33/
      make
      cd ..
      rm v1.33.tar.gz
   fi
   echo "DONE checking and installing sickle!"
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

#Installing prinseq
install_prinseq(){
   echo "Checking and installing prinseq"
   if [ -d "${TOOLS_FOLDER}/prinseq-lite-0.20.4" ]; then
      echo "prinseq already installed"
   else
      echo "Installing prinseq"
      cd $TOOLS_FOLDER
      wget https://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz
      tar xvf prinseq-lite-0.20.4.tar.gz
      rm prinseq-lite-0.20.4.tar.gz
   fi
   echo "DONE checking and installing prinseq!"
}



qc_main
