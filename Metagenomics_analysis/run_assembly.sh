#!/usr/bin/evn bash

### Assembly Analysis


assembly_main(){
   check_and_install
   run_assembly
   assembly_stats
   filter_on_assembly_length
}

# Check and install missing packages for the QC pipeline
check_and_install(){
   install_quast
   install_megahit
   install_spades
   #install_ray
   install_bbmap
   install_idba
}



# Run assembly 
run_assembly(){
   echo "Running all assemblies"

   echo "Running assembly"
   mkdir -p ${ANALYSIS_FOLDER}/assembly

   # Megahit 
   run_assem_megahit

   # Spades
   run_assem_spades

   # Ray
   # run_assem_ray

   # IDBA
   run_assem_idba

   echo "DONE running all assemblies!"
}

run_assem_idba(){
   #Assembly with IDBA
   #first covert paired into single
   #idba_ud -b y /usr/bin/idba_ud -r reads12.fas --num_threads 14 -o ${ANALYSIS_FOLDER}/assembly/idba_ud_out/${i%.1*}
   #1. Merged the paired end data and covert to fasta format.
   echo "Running idba/fq2fa"
   mkdir -p ${ANALYSIS_FOLDER}/assembly/idba_ud
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      ${TOOLS_FOLDER}/idba-master/bin/fq2fa --merge ${i} ${i%1*}2.unmerged.final.clean.fq ${i%1*}12.final.clean.fa

      ${TOOLS_FOLDER}/idba-master/bin/idba_ud \
      -l ${i%1*}12.final.clean.fa \
      --mink 20 \
      --maxk 124 \
      --num_threads 16 \
      -o ${ANALYSIS_FOLDER}/assembly/idba_ud/$(basename ${i%.1*})
   done
}

run_assem_ray(){
   mkdir -p ${ANALYSIS_FOLDER}/assembly/ray
   #Assembly with Ray
   echo "Running ray"
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   for i in ${finallist}
   do
      mpiexec -n 16 \
      ray \
      -k 31 \
      -p ${i} ${i%1*}2.unmerged.final.clean.fq \
      -s ${i%1*}u.final.clean.fq  \
      -o ${ANALYSIS_FOLDER}/assembly/ray/ray_31_assembly/$(basename ${i%.1*})
   done
}

run_assem_spades(){
   mkdir -p ${ANALYSIS_FOLDER}/assembly/spades
   # Assembly with spades
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   echo "Running spades"
   for i in ${finallist}
   do
      ${TOOLS_FOLDER}/SPAdes-3.13.0-Linux/bin/spades.py \
      --meta \
      -1 ${i} \
      -2 ${i%1*}2.unmerged.final.clean.fq \
      -s ${i%1*}u.final.clean.fq \
      -k auto \
      -o ${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*}) \
      -t 24
   done
}

run_assem_megahit(){
   mkdir -p ${ANALYSIS_FOLDER}/assembly/megahit
   # Assembly Megahit
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   #running megahit
   echo "Running megahit"
   for i in ${finallist}
   do
     #To run a co-assembly using multiple samples, you group together the R1 and R2 reads as shown for three libraries (six read files):
     #megahit [options] {-1 pe1_R1.fq,pe2_R1.fq,pe3_R1.fq -2 pe1_R2.fq,pe2_R2.fq,pe3_R2.fq} -o outdir
      ${TOOLS_FOLDER}/MEGAHIT-1.2.1-beta-Linux-static/bin/megahit \
      -1 ${i} \
      -2 ${i%1*}2.unmerged.final.clean.fq \
      -r ${i%1*}u.final.clean.fq \
      -t 24 \
      --presets meta-large \
      --mem-flag 2 \
      -o ${ANALYSIS_FOLDER}/assembly/megahit/$(basename ${i%.1*})
   done

   cat ${ANALYSIS_FOLDER}/assembly/megahit/$(basename ${i%.1*})/final.contigs.fa | cut -d ' ' -f 1 > ${ANALYSIS_FOLDER}/assembly/megahit/$(basename ${i%.1*})/megahit_final.contigs.fa
}

assembly_stats(){
   echo "Running stats"
   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   # creating stats of all the assemblers
   #Stats with quast
   mkdir -p ${ANALYSIS_FOLDER}/assembly/assembly_stats
   for i in ${finallist}
   do
      ${TOOLS_FOLDER}/quast-5.0.2/metaquast.py \
      --threads 14 \
      --gene-finding \
      -o $ANALYSIS_FOLDER/assembly/assembly_stats/$(basename ${i%.1*}) \
      -l megahit,SPAdes,IDBA-UD \
      $ANALYSIS_FOLDER/assembly/megahit/$(basename ${i%.1*})/megahit_final.contigs.fa \
      $ANALYSIS_FOLDER/assembly/spades/$(basename ${i%.1*})/contigs.fasta \
      ${ANALYSIS_FOLDER}/assembly/idba_ud/$(basename ${i%.1*})/contig.fa
      #${ANALYSIS_FOLDER}/assembly/ray/$(basename ${i%.1*})/ray_31_assembly/  \
   done
   echo "DONE running stats!"
}


filter_on_assembly_length(){
   echo "Running assembly length filter"

   finallist=$(ls -d ${ANALYSIS_FOLDER}/QC/final_QC_output/*.1.unmerged.final.clean.fq | awk '{print $NF}')
   #running megahit
   for i in ${finallist}
   do
      #spades
      ${TOOLS_FOLDER}/bbmap/reformat.sh \
      in=${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs.fasta \
      out=${ANALYSIS_FOLDER}/assembly/spades/$(basename ${i%.1*})/contigs_1000_filtered.fasta \
      minlength=1000

      #megahit
      ${TOOLS_FOLDER}/bbmap/reformat.sh \
      in=${ANALYSIS_FOLDER}/assembly/megahit/$(basename ${i%.1*})/megahit_final.contigs.fa \
      out=${ANALYSIS_FOLDER}/assembly/megahit/$(basename ${i%.1*})/megahit_final_1000_filtered.fasta \
      minlength=1000

      #idba_ud
      ${TOOLS_FOLDER}/bbmap/reformat.sh \
      in=${ANALYSIS_FOLDER}/assembly/idba_ud/$(basename ${i%.1*})/contig.fa \
      out=${ANALYSIS_FOLDER}/assembly/idba_ud/$(basename ${i%.1*})/contig_1000_filtered.fa \
      minlength=1000
   done

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

assembly_main
