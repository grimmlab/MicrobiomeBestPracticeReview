# bacteria data. Most of the parameters are used as default and in case you
# would like to add or remove any command or parameter, you can edit it in
# here. Running this bash script is simply produce output which will then
# be used for plotting for more details on workflow please visit
# https://www.mothur.org/wiki/MiSeq_SOP

mothur_16S_workflow(){
   check_and_install
   run_set_paths  #do not uncomment this function
   run_fastqc
   run_mothur_preprocessing
   run_mothur_alignment
   run_mothur_post_alignment_quality_check_1
   run_mothur_post_alignment_quality_check_uchime
   run_mothur_post_alignment_quality_check_2
   run_mothur_classify_and_cluster_seq_to_OTUs
   run_mothur_phylogenetic_analysis
   run_mothur_downstream_analysis
   run_modify_phyliptre
   run_R_plotting
   cd $MYCURRDIR
}


check_and_install(){
   echo "Running Check and install for mothur workflow"

   if [ -d "${TOOLS_FOLDER}/mothur" ]; then
      echo "mothur already installed"
   else
      cd $TOOLS_FOLDER
      wget https://github.com/mothur/mothur/releases/download/v.1.42.3/Mothur.linux_64_noReadline.zip
      unzip Mothur.linux_64_noReadline.zip -d $TOOLS_FOLDER
      rm Mothur.linux_64.noReadline.zip
   fi

   #fastqc
   cd $TOOLS_FOLDER
   if [ -d "${TOOLS_FOLDER}/FastQC" ]; then
      echo "FastQC already installed"
   else
      FASTQC_VERSION=v0.11.8
      wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_${FASTQC_VERSION}.zip
      unzip fastqc_${FASTQC_VERSION}.zip
      cd $TOOLS_FOLDER/FastQC
      chmod 755 fastqc
      cd ..
      rm fastqc_${FASTQC_VERSION}.zip
   fi

   echo "DONE running check and install for mothur workflow"
}

run_fastqc(){
   echo "Running fastqc"

   cd ${ANALYSIS_FOLDER}
   mkdir -p ${ANALYSIS_FOLDER}/QC
   #chmod 755 fastqc
   find ${RAWDATA_FOLDER}/ -name "*.fastq" | xargs -n 1 ${TOOLS_FOLDER}/FastQC/fastqc -o QC/
   cd ${ANALYSIS_FOLDER}/QC
   #unzip *.zip

   echo "DONE running fastqc"
}

run_set_paths(){
  MYCURRDIR=$PWD
  mkdir -p $ANALYSIS_FOLDER/mothur/mothur_output
  cd  $ANALYSIS_FOLDER/mothur/mothur_output
  PROCESSORS=8
  PROJECT_NAME="MiSeq_16S"
}

run_mothur_preprocessing(){
   echo "Running mothur preprocessing"

   ${TOOLS_FOLDER}/mothur/mothur \
   "#make.file(inputdir=${RAWDATA_FOLDER}/, \
   type=fastq, \
   prefix=${PROJECT_NAME})"

   ${TOOLS_FOLDER}/mothur/mothur \
   "#make.contigs(file=$RAWDATA_FOLDER/${PROJECT_NAME}.files, \
   inputdir=${RAWDATA_FOLDER}/, outputdir=${ANALYSIS_FOLDER}/mothur/mothur_output/, \
   processors=${PROCESSORS})"

    

   cd ${ANALYSIS_FOLDER}/mothur/mothur_output/

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

   echo "DONE running mothur preprocessing"
}


run_mothur_alignment(){
   echo "Running mothur alignment"

   ${TOOLS_FOLDER}/mothur/mothur \
   "#align.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.fasta, \
   reference=${REFERENCE_FOLDER}/silva/silva.nr_v132.align, \
   processors=${PROCESSORS})"

   echo "DONE running mothur alignment"
}

run_mothur_post_alignment_quality_check_1(){
   echo "Running mothur post alignment quality check"

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

}

run_mothur_post_alignment_quality_check_uchime(){

   echo "Running uchime command"
   ${TOOLS_FOLDER}/mothur/mothur \
   "#chimera.uchime(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.fasta, \
   count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.count_table, \
   dereplicate=t)"#, \
   #processors=${PROCESSORS})"
}

run_mothur_post_alignment_quality_check_2(){
   echo "Running mothur post alignment quality check"
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

   echo "DONE running mothur post alignment quality check"
}

run_mothur_classify_and_cluster_seq_to_OTUs(){
   echo "Running mothur classify and cluster seq to OTUs"

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

   echo "DONE running mothur classify and cluster seq to OTUs"
}


run_mothur_phylogenetic_analysis(){
   echo "Running mothur phylogenetic analysis"
	
   ${TOOLS_FOLDER}/mothur/mothur \
   "#dist.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, \
   output=lt)"

   ${TOOLS_FOLDER}/mothur/mothur \
   "#clearcut(phylip=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.dist)"

   echo "DONE running mothur phylogenetic analysis"
}


run_mothur_downstream_analysis(){
   echo "Running mothur downstream analysis"

   ${TOOLS_FOLDER}/mothur/mothur \
   "#rarefaction.single(shared=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared)"

   ${TOOLS_FOLDER}/mothur/mothur \
   "#make.biom(shared=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared,  \
   label=0.03, \
   reftaxonomy=${REFERENCE_FOLDER}/greengenes/gg_13_8_99.gg.tax, \
   constaxonomy=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy)"

   echo "DONE running mothur downstream analysis"
}

echo "Path check.................... $pwd"
run_modify_phyliptre(){
   echo "Running mothur modify phyliptre"

   ${TOOLS_FOLDER}/mothur/mothur \
   "#get.oturep(phylip=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.dist,\
   list=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, \
   count=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, \
   fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta)"

   python ${BIN_FOLDER}/modify_phylip.py \
   ${ANALYSIS_FOLDER}/mothur/mothur_output/${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.fasta

   ${TOOLS_FOLDER}/mothur/mothur \
   "#dist.seqs(fasta=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.otu_modified.fasta, \
   output=lt)"

   ${TOOLS_FOLDER}/mothur/mothur \
   "#clearcut(phylip=${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.otu_modified.phylip.dist)"

   echo "DONE running mothur modify phyliptre"
}

#To run this function Phyloseq pacakge should be installed
run_R_plotting(){
   echo "Running R plotting"
   
   mkdir -p ${ANALYSIS_FOLDER}/mothur/mothur_output/plots/
   Rscript ${BIN_FOLDER}/plots.R \
   -l ${ANALYSIS_FOLDER}/mothur/mothur_output/${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list \
   -tx ${ANALYSIS_FOLDER}/mothur/mothur_output/${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy \
   -t ${ANALYSIS_FOLDER}/mothur/mothur_output/${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.otu_modified.phylip.tre \
   -s ${ANALYSIS_FOLDER}/mothur/mothur_output/${PROJECT_NAME}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared \
   -m ${RAWDATA_FOLDER}/../metadata/mouse.dpw.metadata \
   -o ${ANALYSIS_FOLDER}/mothur/mothur_output/plots/

   echo "DONE running R plotting"
}

mothur_16S_workflow
