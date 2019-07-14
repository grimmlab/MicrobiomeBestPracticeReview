# bacteria data. Most of the parameters are used as default and in case you
# would like to add or remove any command or parameter, you can edit it in
# here. Running this bash script is simply produce output which will then
# be used for plotting for more details on workflow please visit
# https://www.mothur.org/wiki/MiSeq_SOP

NAME=Test_Mothur_16S_workflow

mothur_16S_workflow(){
  install_tools
  download_testdata
  download_reference_database
  run_fastqc
  run_set_paths #don't comment this function
  run_mothur_preprocessing
  run_mothur_alignment
  run_mothur_post_aligment_quality_check
  run_mothur_classify_and_cluster_seq_to_OTUs
  run_mothur_phylogenetic_analysis
  run_mothur_downstream_analysis
  run_modify_phyliptre
  #run_R_plotting
}

#DATE=$(date +%Y-%m-%d)
ROOT_FOLDER_NAME=$NAME

for FOLDER in analysis tools rawdata reference
do
    mkdir -p $ROOT_FOLDER_NAME/$FOLDER
done

TOOLS_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/tools
RAWDATA_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/rawdata
ANALYSIS_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/analysis
REFERENCE_FOLDER=$(pwd)/$ROOT_FOLDER_NAME/reference
BIN_FOLDER=$(pwd)/bin

install_tools(){
# Download and compile [Boost](https://www.boost.org/users/download/)
cd $TOOLS_FOLDER
wget https://github.com/mothur/mothur/releases/download/v1.40.5/Mothur.linux_64.noReadLine.zip
unzip Mothur.linux_64.noReadLine.zip -d $TOOLS_FOLDER
rm Mothur.linux_64.noReadLine.zip
#sudo ln -s $TOOLS_FOLDER/mothur/mothur /usr/local/bin/mothur

#fastqc
cd $TOOLS_FOLDER
FASTQC_VERSION=v0.11.8
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_${FASTQC_VERSION}.zip
unzip fastqc_${FASTQC_VERSION}.zip
cd $TOOLS_FOLDER/FastQC
chmod 755 fastqc
rm fastqc_${FASTQC_VERSION}.zip
}


download_testdata(){
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
  #cd ${RAWDATA_FOLDER}
  #S=0
  #for FILE in ERR20420{42..46}
  #do
  	#wget -r -l2 -A.gz ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR204/
  #S=$(($S+1))
  #wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR204/${FILE}/S${S}_R1_001.fastq.gz
  #wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR204/${FILE}/S${S}_R2_001.fastq.gz
  #done

  cd ${RAWDATA_FOLDER}
  wget www.mothur.org/w/images/d/d6/MiSeqSOPData.zip
  unzip MiSeqSOPData.zip
  rm Mock_*
  rm HMP_MOCK.v35.fasta*
  rm -rf MiSeqSOPData.zip
  cd ${RAWDATA_FOLDER}/MiSeqSOPData
  rm Mock_S280_*
  rm HMP_MOCK.v35.fasta
}


download_reference_database(){
## Download the references and taxonomy
#This [blog](http://blog.mothur.org/2018/01/10/SILVA-v132-reference-files/) will help you to understand how reference can be formatted for mothur. Here we are downloading already formatted [SILVA](https://mothur.org/wiki/Silva_reference_files#Release_132) reference by mothur with following steps below

cd $REFERENCE_FOLDER
REFERENCE_VERSION=v132
mkdir -p ${REFERENCE_FOLDER}/silva
cd  ${REFERENCE_FOLDER}/silva
wget https://mothur.org/w/images/3/32/Silva.nr_${REFERENCE_VERSION}.tgz
tar -xzvf Silva.nr_${REFERENCE_VERSION}.tgz
rm Silva.nr_${REFERENCE_VERSION}.tgz


#RDP reference files can be download from [here](https://www.mothur.org/wiki/RDP_reference_files)
cd $REFERENCE_FOLDER
RPD_VERSION=16
wget https://www.mothur.org/w/images/c/c3/Trainset${RPD_VERSION}_022016.pds.tgz
tar -xzvf Trainset${RPD_VERSION}_022016.pds.tgz
rm Trainset${RPD_VERSION}_022016.pds.tgz


#Greengenes reference files can be downloaded from  [here](https://www.mothur.org/wiki/Greengenes-formatted_databases)
cd ${REFERENCE_FOLDER}
mkdir -p ${REFERENCE_FOLDER}/greengenes
cd  ${REFERENCE_FOLDER}/greengenes
wget http://www.mothur.org/w/images/6/68/Gg_13_8_99.taxonomy.tgz
tar -xzvf Gg_13_8_99.taxonomy.tgz
rm Gg_13_8_99.taxonomy.tgz
}


run_fastqc(){
  cd ${ANALYSIS_FOLDER}
  mkdir -p ${ANALYSIS_FOLDER}/QC
  find ${RAWDATA_FOLDER}/MiSeq_SOP/ -name "*.fastq" | xargs -n 1 ${TOOLS_FOLDER}/FastQC/fastqc -o QC/
  cd ${ANALYSIS_FOLDER}/QC
  #unzip *.zip
}

run_set_paths(){
  mkdir -p $ANALYSIS_FOLDER/mothur_output
  cd  $ANALYSIS_FOLDER/mothur_output
  PROCESSORS=8
  PROJECT_NAME="MiSeq_16S"
}

run_mothur_preprocessing(){
  ${TOOLS_FOLDER}/mothur/mothur \
  "#make.file(inputdir=${RAWDATA_FOLDER}/MiSeq_SOP/, \
  type=fastq, \
  prefix=${PROJECT_NAME})"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#make.contigs(file=$RAWDATA_FOLDER/MiSeq_SOP/${PROJECT_NAME}.files, \
  inputdir=${RAWDATA_FOLDER}/MiSeq_SOP, outputdir=${ANALYSIS_FOLDER}/mothur_output/, \
  processors=${PROCESSORS})"

  # Count number of reads after contigs
  ${TOOLS_FOLDER}/mothur/mothur \
  "#count.groups(group=${PROJECT_NAME}.contigs.groups)"

  # Summary of contigs created with above command
  ${TOOLS_FOLDER}/mothur/mothur \
  "#summary.seqs(fasta=${PROJECT_NAME}.trim.contigs.fasta, \
  processors=${PROCESSORS})"

  # require manual intervention if the statistics is not good
 #For V4 region
  ${TOOLS_FOLDER}/mothur/mothur \
  "#screen.seqs(fasta=${PROJECT_NAME}.trim.contigs.fasta, \
  group=${PROJECT_NAME}.contigs.groups, \
  summary=${PROJECT_NAME}.trim.contigs.summary, \
  maxambig=0, \
  maxlength=275, \
  maxhomop=8)"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#summary.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.fasta, \
  processors=${PROCESSORS})"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#unique.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.fasta, \
  format=name)"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#count.seqs(name=${PROJECT_NAME}.trim.contigs.good.names, \
  group=${PROJECT_NAME}.contigs.good.groups)"
}


run_mothur_alignment(){
  ${TOOLS_FOLDER}/mothur/mothur \
  "#align.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.fasta, \
  reference=${REFERENCE_FOLDER}/silva/silva.nr_v132.align, \
  processors=${PROCESSORS})"
}

run_mothur_post_aligment_quality_check(){
  ${TOOLS_FOLDER}/mothur/mothur \
  "#summary.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.align, \
  count=${PROJECT_NAME}.trim.contigs.good.count_table, \
  processors=${PROCESSORS})"

  #check for the positions, If its V1-V3 then start=6388, end=13861 and for V3-V4 region start=13861,end=25319
  #For V4 region start=11894, end=25319
  #Require manual intervention if data is not v4 region
  ${TOOLS_FOLDER}/mothur/mothur \
  "#screen.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.align, \
  count=${PROJECT_NAME}.trim.contigs.good.count_table, \
  summary=${PROJECT_NAME}.trim.contigs.good.unique.summary, \
  start=13862, \
  end=23444, \
  maxhomop=6, \
  processors=${PROCESSORS})"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#summary.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.align, \
  count=${PROJECT_NAME}.trim.contigs.good.good.count_table, \
  processors=${PROCESSORS})"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#filter.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.align,\
  processors=${PROCESSORS}, \
  vertical=T, \
  trump=.)"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#unique.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.fasta, \
  count=${PROJECT_NAME}.trim.contigs.good.good.count_table)"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#pre.cluster(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.fasta,\
  count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.count_table, \
  diffs=2, \
  processors=${PROCESSORS})"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#chimera.uchime(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.fasta, \
  count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.count_table, \
  dereplicate=t, \
  processors=${PROCESSORS})"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#remove.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.fasta, \
  count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.count_table, \
  accnos=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos)"


  ${TOOLS_FOLDER}/mothur/mothur \
  "#summary.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, \
  count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.count_table, \
  processors=${PROCESSORS})"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#classify.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,\
  count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, \
  reference=${REFERENCE_FOLDER}/silva/silva.nr_v132.align, \
  taxonomy=${REFERENCE_FOLDER}/silva/silva.nr_v132.tax, \
  cutoff=80, \
  processors=${PROCESSORS})"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#remove.lineage(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, \
  count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, \
  taxonomy=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, \
  taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#summary.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, \
  count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, \
  processors=${PROCESSORS})"
}

run_mothur_classify_and_cluster_seq_to_OTUs(){
### In this section sequences are clustered and classifies into OTUs
  $TOOLS_FOLDER/mothur/mothur \
  "#dist.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta)"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#cluster.split(column=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, \
  count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, \
  taxonomy=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, \
  splitmethod=classify, \
  taxlevel=4, \
  cutoff=0.03, \
  method=opti, \
  processors=${PROCESSORS})"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#make.shared(list=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, \
  count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, \
  label=0.03)"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#classify.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, \
  count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, \
  template=${REFERENCE_FOLDER}/greengenes/gg_13_8_99.fasta, \
  taxonomy=${REFERENCE_FOLDER}/greengenes/gg_13_8_99.gg.tax, \
  cutoff=80, \
  processors=${PROCESSORS})"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#classify.otu(list=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list,\
   taxonomy=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.gg.wang.taxonomy, \
   count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, \
   label=0.03, \
   cutoff=80, \
   basis=otu, \
   probs=F)"
}


run_mothur_phylogenetic_analysis(){
  ${TOOLS_FOLDER}/mothur/mothur \
  "#dist.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, \
  output=lt)"

   ${TOOLS_FOLDER}/mothur/mothur \
   "#clearcut(phylip=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.dist)"
}


run_mothur_downstream_analysis(){
  ${TOOLS_FOLDER}/mothur/mothur \
  "#rarefaction.single(shared=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#make.biom(shared=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared,  \
  label=0.03, \
  reftaxonomy=${REFERENCE_FOLDER}/greengenes/gg_13_8_99.gg.tax, \
  constaxonomy=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy)"
}

run_modify_phyliptre(){
  ${TOOLS_FOLDER}/mothur/mothur \
  "#get.oturep(phylip=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.dist,\
  list=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, \
  count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, \
  fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta)"

  python ${BIN_FOLDER}/modify_phylip.py \
  ${ANALYSIS_FOLDER}/mothur_output/${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.fasta

  ${TOOLS_FOLDER}/mothur/mothur \
  "#dist.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.otu_modified.fasta, \
  output=lt)"

  ${TOOLS_FOLDER}/mothur/mothur \
  "#clearcut(phylip=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.otu_modified.phylip.dist)"
}

#To run this function Phyloseq pacakge should be installed
run_R_plotting(){
    mkdir -p ${ANALYSIS_FOLDER}/plots
    Rscript ${BIN_FOLDER}/plot.R \
    -l ${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list \
    -tx ${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy \
    -t ${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.otu_modified.phylip.tre \
    -s ${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared \
    -m ${RAWDATA_FOLDER}/MiSeq_SOP/mouse.dpw.metadata
    -o ${ANALYSIS_FOLDER}/plots
}

mothur_16S_workflow
