### Workflow for whole genome shotgun metagenomics analysis

#!/usr/bin/env bash

main(){
    create_folders
    set_variables
    fastqc_stats
    trimmomatic_sickle_QC
    bbmap_QC1
    downloads_indexing
    bbmap_QC2
    creating_stats
}


### Prerequisite
create_folders(){
# Creating folder structure
NAME=Metagenomic_QC
#DATE=$(date +%Y-%m-%d)
#ROOT_FOLDER_NAME=$DATE-$NAME
ROOT_FOLDER_NAME=$NAME
for FOLDER in analysis tools rawdata reference
do
    mkdir -p $ROOT_FOLDER_NAME/$FOLDER
done
}


# setting variable path
set_variables(){
TOOLS_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/tools
RAWDATA_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/rawdata
ANALYSIS_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/analysis
REFERENCE_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/reference
#DONE
}


#### 1. Generating comprehensive report and stats of data quality using fastqc and BBMAP
fastqc_stats(){
#### Run fastqc on our data (non interactively)
cd $ANALYSIS_FOLDER
mkdir -p $ANALYSIS_FOLDER/QC/fastqc
find $RAWDATA_FOLDER -name "*.fastq.gz" | xargs -n 1 $TOOLS_FOLDER/FastQC/fastqc -o $ANALYSIS_FOLDER/QC/fastqc
#DONE

#for rawdata
mkdir -p ${ANALYSIS_FOLDER}/QC/bbmap
rawdatalist=$(ll -d $RAWDATA_FOLDER/*R1* | awk '{print $NF}')
for s in $rawdatalist
  do
    $TOOLS_FOLDER/bbmap/reformat.sh \
    threads=16 \
    in=${s} \
    in2=${s%R1*}R2_2.fastq.gz \
    2>&1 >/dev/null | awk '{print "RAWREADS "$0}' | tee -a $ANALYSIS_FOLDER/QC/bbmap/${s%R1*}stats.txt
  done
}


#### 2. Quality Control Trimming using trimmomatic and Sickle
trimmomatic_sickle_QC(){
#Trimming low quality, short length reads, adapters
mkdir -p $ANALYSIS_FOLDER/QC/trimmomatic
mkdir -p $ANALYSIS_FOLDER/QC/sickle
rawdatalist=$(ll -d $RAWDATA_FOLDER/*R1* | awk '{print $NF}')

for s in $rawdatalist
  do
  fname=$(basename $s | sed -e "s/_R1_1.fastq.gz//")
  #Running trimmomatic
  java -jar $TOOLS_FOLDER/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
  ${s%_R*}_R1_1.fastq.gz ${s%_R*}_R2_2.fastq.gz \
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
  sickle pe \
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
  sickle se \
  -n \
  -f ${ANALYSIS_FOLDER}/QC/trimmomatic/${fname}.1.u.trimmoclean.fq.gz \
  -o ${ANALYSIS_FOLDER}/QC/sickle/${fname}.1.u.trimmoclean.sickleclean.fq \
  -t sanger \
  -q 20 \
  -l 60

  sickle se \
  -n \
  -f ${ANALYSIS_FOLDER}/QC/trimmomatic/${fname}.2.u.trimmoclean.fq.gz \
  -o ${ANALYSIS_FOLDER}/QC/sickle/${fname}.2.u.trimmoclean.sickleclean.fq \
  -t sanger \
  -q 20 \
  -l 60

  #Combining all unpaired files
  cat ${ANALYSIS_FOLDER}/QC/sickle/${fname}.1.u.trimmoclean.sickleclean.fq \
  ${ANALYSIS_FOLDER}/QC/sickle/${fname}.2.u.trimmoclean.sickleclean.fq \
  ${ANALYSIS_FOLD$TOOLS_FOLDER/FastQC/fastqcER}/QC/sickle/${fname}.u.trimmoclean.sickleclean.fq \
   > ${ANALYSIS_FOLDER}/QC/sickle/${fname}.unpaired.trimmoclean.sickleclean.fq

  rm ${ANALYSIS_FOLDER}/QC/sickle/${fname}.1.u.trimmoclean.sickleclean.fq ${ANALYSIS_FOLDER}/QC/sickle/${fname}.2.u.trimmoclean.sickleclean.fq ${ANALYSIS_FOLDER}/QC/sickle/${fname}.u.trimmoclean.sickleclean.fq
done
#DONE
}


#### 3. Quality control removing phix adapters and sequencing artifacts using BBMAP
bbmap_QC1(){
# paired data
sicklelist=$(ll -d ${ANALYSIS_FOLDER}/QC/sickle/*1.trimmoclean.sickleclean.fq | awk '{print $NF}')
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
      minl```bashength=60
    done

# unpaired data
sickleunplist=$(ll -d ${ANALYSIS_FOLDER}/QC/sickle/*unpaired.trimmoclean.sickleclean.fq| awk '{print $NF}')
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
      #DONE
    }

#### 3. Downloading and indexing the Host genome (human and mouse) using prinseq for cleaning and BBMAP for indexing the genome
#Creating human reference database
downloads_indexing(){
mkdir -p $REFERENCE_FOLDER/human
cd $REFERENCE_FOLDER/human
for i in {1..22} X Y MT
do
 wget ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/hs_ref_GRCh38.p12_chr$i.fa.gz
done

for i in {1..22} X Y MT
do
  gzip -dvc hs_ref_GRCh38.p12_chr$i.fa.gz >>$REFERENCE_FOLDER/human/hsref_GRCh38_p12.fa
done

#Clean
perl ${TOOLS_FOLDER}/prinseq-lite-0.20.4/prinseq-lite.pl \
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

#Creating mouse reference database
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
perl ${TOOLS_FOLDER}/prinseq-lite-0.20.4/prinseq-lite.pl \
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

#indexing the genome
$TOOLS_FOLDER/bbmap/bbmap.sh \
ref=${REFERENCE_FOLDER}/human/hsref_GRCh38_p12_clean.fasta \
path=${REFERENCE_FOLDER}/human

$TOOLS_FOLDER/bbmap/bbmap.sh \
ref=${REFERENCE_FOLDER}/mouse/mmref_GRCm38_p4_clean.fasta \
path=${REFERENCE_FOLDER}/mouse
}


#### 4. Removing the host contamination and generating the stats of the data using BBMAP
bbmap_QC2(){
bbduklist=$(ll -d ${ANALYSIS_FOLDER}/QC/bbmap/*.1.trimmoclean.sickleclean.bbdukclean.fq | awk '{print $NF}')

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
    minhit```bashs=2 \
    qtrim=rl \
    trimq=20 \
    minlength=60 \
    in=${s},${s%1*}unpaired.trimmoclean.sickleclean.bbdukclean.fq \
    in2=${s%1*}2.trimmoclean.sickleclean.bbdukclean.fq,NULL \
    outu1=${s%1*}1.final.clean.fq \
    outu2=${s%1*}2.final.clean.fq \
    outu=${s%1*}u.clean.fq \
    path=${REFERENCE_FOLDER}/mouse/ 2>&1 >/dev/null | awk '{print "HOST CONTAMINATION SEQUENCES "$0}' | tee -a $ANALYSIS_FOLDER/QC/bbmap/${s%1*}stats.txt
  done
#DONE

finallist=$(ll -d ${ANALYSIS_FOLDER}/QC/bbmap/*.1.final.clean.fq | awk '{print $NF}')
for s in $finallist
  do
    $TOOLS_FOLDER/bbmap/bbmerge.sh \
    threads=4 \
    in1=${s} \
    in2=${s%1*}2.final.clean.fq \
    out=${s%1*}merged.final.clean.fq \
    outu1=${s%1*}1.merged.final.clean.fq \
    outu2=${s%1*}2.merged.final.clean.fq \
    mininsert=60 \
    2>&1 >/dev/null | awk '{print "MERGED "$0}' | tee -a $ANALYSIS_FOLDER/QC/bbmap/${s%1*}stats.txt
  done
  #DONE

for s in $finallist
  do
    cat ${s%1*}merged.final.clean.fq  ${s%1*}u.clean.fq  > ${s%1*}u.final.clean.fq
  done
  #DONE

for s in $finallist
  do
    $TOOLS_FOLDER/bbmap/reformat.sh \
    threads=16 \
    in=${s%1*}u.final.clean.fq \
    2>&1 >/dev/null | awk '{print "UNPAIRED "$0}' | tee -a $ANALYSIS_FOLDER/QC/bbmap/${s%1*}stats.txt
  done


for s in $finallist
  do
    $TOOLS_FOLDER/bbmap/reformat.sh \
    threads=16 \
    in1=${s%1*}1.merged.final.clean.fq \
    in2=${s%1*}2.merged.final.clean.fq  2>&1 >/dev/null | awk '{print "PAIRED "$0}' | tee -a ${s%1*}stats.txt
  done
  #DONE
}


#Creating read count stats for cleaned and raw data
creating_stats(){
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
}

main
