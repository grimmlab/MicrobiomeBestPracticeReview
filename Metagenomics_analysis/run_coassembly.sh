#!/usr/bin/evn bash

### Assembly Analysis


coassembly_main(){
   check_and_install
   run_coassembly
   #assembly_stats
   filter_on_assembly_length
}

# Check and install missing packages for the QC pipeline
check_and_install(){
   install_quast
   install_megahit
   install_spades
   #install_ray
   install_idba
}



# Run assembly 
run_coassembly(){
   echo "Running all coassemblies"

   echo "Running coassembly"
   mkdir -p ${ANALYSIS_FOLDER}/coassembly

   # Megahit 
   run_coassem_megahit

   # Spades
   run_coassem_spades

   # Ray
   # run_coassem_ray

   # IDBA
   #run_coassem_idba

   echo "DONE running all coassemblies!"
}

run_coassem_idba(){
   #Assembly with IDBA
   #first covert paired into single
   #1. Merged the paired end data and covert to fasta format.
   echo "Running idba/fq2fa"
   mkdir -p ${ANALYSIS_FOLDER}/coassembly/idba_ud
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      ${TOOLS_FOLDER}/idba-master/bin/fq2fa \
      --merge ${i} ${i%1*}2.unmerged.final.clean.fq \
      ${i%1*}12.final.clean.fa

      ${TOOLS_FOLDER}/idba-master/bin/idba_ud \
      -l ${i%1*}12.final.clean.fa \
      --mink 20 \
      --maxk 124 \
      --num_threads 16 \
      -o ${ANALYSIS_FOLDER}/coassembly/idba_ud/$(basename ${i%.1*})
   done
}

run_coassem_ray(){
   mkdir -p ${ANALYSIS_FOLDER}/coassembly/ray
   #Assembly with Ray
   echo "Running ray"
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      mpiexec -n 16 \
      ray \
      - k 31 \
      -p ${i} ${i%1*}2.unmerged.final.clean.fq \
      -s ${i%1*}u.final.clean.fq  \
      -o ${ANALYSIS_FOLDER}/coassembly/ray/ray_31_assembly/$(basename ${i%.1*})
   done
}

run_coassem_spades(){
   mkdir -p ${ANALYSIS_FOLDER}/coassembly/spades
   # Assembly with spades
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   echo "Running spades"
   #   spades.py --pe1-1 lib1_forward_1.fastq --pe1-2 lib1_reverse_1.fastq \
   #   --pe1-1 lib1_forward_2.fastq --pe1-2 lib1_reverse_2.fastq \
   #   -o spades_output
   cospd_1="${TOOLS_FOLDER}/SPAdes-3.13.0-Linux/bin/spades.py --meta "
   for i in ${finallist}
   do
      cospd="$cospd --pe1-1 ${i} --pe1-2 ${i%1*}2.unmerged.final.clean.fq --pe-s ${i%1*}u.final.clean.fq "
   done

   cospd_2=" -k auto -o ${ANALYSIS_FOLDER}/coassembly/spades  -t 16"

   $cospd_1 $cospd $cospd_2
}

run_coassem_megahit(){
   mkdir -p ${ANALYSIS_FOLDER}/coassembly/megahit
   # Assembly Megahit
   echo "Running megahit"
   R1s=`ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
   R2s=`ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.2.unmerged.final.clean.fq | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
   Rs=`ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.u.final.clean.fq | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`

   ${TOOLS_FOLDER}/MEGAHIT-1.2.1-beta-Linux-static/bin/megahit \
   --pe-1 ${R1s} \
   --pe-2 ${R2s} \
   --pe-s ${Rs} \
   -t 16 \
   --presets meta-large \
   --mem-flag 2 \
   -o ${ANALYSIS_FOLDER}/coassembly/megahit/


}

assembly_stats(){
   echo "Running stats"
   # creating stats of all the assemblers
   #Stats with quast
   mkdir -p ${ANALYSIS_FOLDER}/coassembly/coassembly_stats

   ${TOOLS_FOLDER}/quast-5.0.2/metaquast.py \
   --threads 14 \
   --gene-finding \
   -o $ANALYSIS_FOLDER/coassembly/coassembly_stats \
   -l megahit,SPAdes,IDBA-UD \
   $ANALYSIS_FOLDER/coassembly/megahit/megahit_final.contigs.fa \
   $ANALYSIS_FOLDER/coassembly/spades/contigs.fasta \
   ${ANALYSIS_FOLDER}/coassembly/idba_ud/contig.fa
   #${ANALYSIS_FOLDER}/coassembly/ray/$(basename ${i%.1*})/ray_31_assembly/  \

   echo "DONE running stats!"
}


filter_on_assembly_length(){
   echo "Running assembly length filter"

   #running megahit
   #spades
   ${TOOLS_FOLDER}/bbmap/reformat.sh \
   in=${ANALYSIS_FOLDER}/coassembly/spades/contigs.fasta \
   out=${ANALYSIS_FOLDER}/coassembly/spades/contigs_1000_filtered.fasta \
   minlength=1000

   #megahit
   ${TOOLS_FOLDER}/bbmap/reformat.sh \
   in=${ANALYSIS_FOLDER}/coassembly/megahit/megahit_final.contigs.fa \
   out=${ANALYSIS_FOLDER}/coassembly/megahit/megahit_final_1000_filtered.fasta \
   minlength=1000

   #idba_ud
   #${TOOLS_FOLDER}/bbmap/reformat.sh \
   #in=${ANALYSIS_FOLDER}/coassembly/idba_ud/contig.fa \
   #out=${ANALYSIS_FOLDER}/coassembly/idba_ud/contig_1000_filtered.fa \
   #minlength=1000

   echo "DONE running assembly length filter!"

}



#Installing quast 
install_quast(){
   echo "Checking and installing quast"
   if [ -d "${TOOLS_FOLDER}/quast-5.0.2" ]; then
      echo "quast already installed"
   else
      echo "Installing quast"
      cd $TOOLS_FOLDER
      wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz
      tar -xzf quast-5.0.2.tar.gz
      rm quast-5.0.2.tar.gz
   fi
   echo "DONE checking and installing quast!"
}


#Installing megahit
install_megahit(){
   echo "Checking and installing MEGAHIT"
   if [ -d "${TOOLS_FOLDER}/MEGAHIT-1.2.1-beta-Linux-static" ]; then
      echo "MEGAHIT already installed"
   else
      echo "Installing MEGAHIT"
      cd $TOOLS_FOLDER
      wget https://github.com/voutcn/megahit/releases/download/v1.2.1-beta/MEGAHIT-1.2.1-beta-Linux-static.tar.gz
      tar zvxf MEGAHIT-1.2.1-beta-Linux-static.tar.gz
      rm MEGAHIT-1.2.1-beta-Linux-static.tar.gz
   fi
   echo "DONE checking and installing MEGAHIT!"
}

#Installing spades
install_spades(){
   echo "Checking and installing SPAdes"
   if [ -d "${TOOLS_FOLDER}/SPAdes-3.13.0-Linux" ]; then
      echo "SPAdes already installed"
   else
      echo "Installing SPAdes"
      cd $TOOLS_FOLDER
      wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz
      tar zvxf SPAdes-3.13.0-Linux.tar.gz
      rm SPAdes-3.13.0-Linux.tar.gz
   fi
   echo "DONE checking and installing SPAdes!"
}


#Installing ray 
install_ray(){
   echo "Checking and installing Ray"
   if [ -d "${TOOLS_FOLDER}/Ray-2.3.1" ]; then
      echo "Ray already installed"
   else
      echo "Installing Ray"
      cd $TOOLS_FOLDER
      wget https://sourceforge.net/projects/denovoassembler/files/Ray-2.3.1.tar.bz2/download 
      mv download Ray-2.3.1.tar.bz2
      tar xjf Ray-2.3.1.tar.bz2
      cd Ray-2.3.1
      make
      cd ..
      rm Ray-2.3.1.tar.bz2
      #cp Ray /home/user/bin
   fi
   echo "DONE checking and installing Ray!"
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


coassembly_main
