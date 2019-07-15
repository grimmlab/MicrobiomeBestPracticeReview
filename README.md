## Introduction

This review paper aims to provide a comprehensive workflow to perform amplicon and shotgun metagenomics analysis. There are two workflows provided. First workflow for amplicon,  using the standard mothur workflow, and along with it visualization methods for the mothur output. Second workflow for metagenomics, using a variety of tools openly available which have been stitched together to form a usable pipeline.

Both the workflows are controlled by bash scripts: amplicon_analysis.sh and metagenomics_analysis.sh. The bash scripts contain functions which call the respective underlying tools. Of-course, the tools have to exist in the system before using them, hence, a function called as "check_and_install" is written into each script which checks if the tools exists in a certain path or not.

Since the workflows utilize so many different tools, it requires quiet a bit of patience for the download and installation process. Please go through the steps below before you begin to use the workflows.

## I. Metagenomic Sequencing Analysis Workflow

## Prerequisites
Although the "check_and_install" function is designed to install and setup the required software on the fly, there are some basic prerequisites that need to be satisfied:

### OS
Any Linux based distro should work. Out test OS is:

Distributor ID: Ubuntu
Description:    Ubuntu 18.04.2 LTS
Release:        18.04
Codename:       bionic

which can be got by typing 'lsb_release -a' on a Ubuntu based system.

### Hardware:
It is no secret that the hardware highly influences the speed of the workflows. The most time consuming tasks are the ones that involve assemblies, reference based alignment. A modest configuration consists of 16+cores and 100 GB of RAM with 1TB of diskspace. A majority of the diskspace is occupied by reference databases like nr database, kraken database, etc. Our HW configuration consists of 20 core CPU with 125 GB.

### Software and packages
Some software should be installed by the user directly as the workflow depends on a lot of external software.
Without these the workflow will fail to run.

1. gcc, g++
2. java
3. python2.7 and python3: pip, libraries (sys, os, shutil, stat, re, time, tarfile, operator, math, Bio, argparse)
4. perl libraries (Bio)
5. R
6. git
7. metabat: Install instructions can be found under https://bitbucket.org/berkeleylab/metabat/src/master/README.md. Metabat should be visible in the system PATH.
8. checkM (checkm-genome): Install instructions can be found under https://github.com/Ecogenomics/CheckM/wiki/Installation.
After installation the checkM database needs to be built using https://data.ace.uq.edu.au/public/CheckM_databases/ and building by using `checkm data setRoot PATH_TO_DOWNLOADED_DATABASE`
**NOTE**: Make sure checkM is placed finally under `/usr/local/bin`

## Steps to run the Metagenomics workflow (metagenomics_analysis.sh)
1. Preparing databases:
```bash
sh prepare_databases.sh
```
Insert the `LINKPATH_DB=/xxx/.../references/ to 'metagenomics_analysis.sh'`
```bash
LINKPATH_DB=/xxx/.../references/
```

The databases form the core of the workflows. Unfortunately, the databases are huge and  take a long time to download and index. If these databases already exist in your system pleease modify the scripts with the correct paths. Otherwise choose the missing databases and run 'prepare_databases.sh' where the databases will be installed under the 'references' in the current directory. At the end of the preparation of databases a path will shown in the stdout which needs to be plug-in to the metabgenomics_analysis.sh script (to LINKPATH_DB variable).

The following databases are installed:

Human and Mouse reference genome:
ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/
ftp://ftp.ncbi.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/

Kraken database:
http://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz
(needs to be indexed with kraken2)

Megan database:
http://ab.inf.uni-tuebingen.de/data/software/megan6/download/

NR database: Non-redundant database can be found at ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
(needs to be index with diamond)

Metaphlan database:
https://bitbucket.org/biobakery/metaphlan2/downloads/mpa_v20_m200.tar
(needs to be built with bowtie2)

checkM database:
Needs to be manually installed (Please check prerequisites).

Humann2 database:
Downloaded using humann2_databases.

2. Add source path for raw data

In the 'metagenomics_analysis.sh' add the path for the raw-data. Please note that the workflow will make a local copy of the rawdata before proceeding further.

```bash
SRC_RAWDATA=/path_to_my_rawdata_samples/.../.../
```

3. Set name of workflow

Next choose an appropriate name for the analysis in the 'metagenomics_analysis.sh' script. All the sub-folders like tools, analysis, rawdata copy, etc will be created under this folder name.

```bash
NAME=MY_METAGENOMIC_ANALYSIS_EXP
```

4. Run the workflow

Finally, the workflow is ready to be run
```bash
sh metagenomics_analysis.sh
```
There are messages on the stdout showing the status and progress of the analysis.

The script consists of several sub-scripts and functions. Each sub-script has its own "check_and_install". The "check_and_install" checks for the tools required to run the respective script and installs them if they are missing.

**NOTE**: The installation of Megan is an interactive installation and requires the user to input Y/N and memory options(between ~8GB-16GB depending on sample size). We recommend to use default options. Megan will be installed in the user home directory.


## Step-by-Step Analysis

The metagenomics workflow is a time-consuming workflow. Hence, the bash scripts are kept as simple as possible. In order to perform only one type of analysis, you can always comment the remaining functions.

For example, the quality control (`run_qc`) can be run only once initially and then commented for any further analysis for reruns.

If the appropriate steps have already been run, then these can be commented and other steps can be run. This is of-course, not true for steps dependent on previous outputs.

## Brief Description of the Each Step
1. run_qc.sh
2. run_assembly.sh
3. run_coassembly.sh
4. run_reference_analysis.sh
5. run_comparative_analysis.sh
6. run_coverage_bining.sh
7. run_binrefinement.sh
8. run_bin_taxonomic_classification.sh
9. run_bin_functional_classification.sh

![Workflow](/data1/Active_Projects/paper_scripts/MicroReviewPaper/Figure4_for_README.tif)


1. Quality control `(run_qc.sh)`: This scripts is running series of steps with different tools to perform quality control. [FastQC](https://github.com/s-andrews/FastQC/) is used to generate comprehensive report of data quality on raw data. Followed by this is a series of steps including removal of adapters, low quality reads, sequencing_artifacts, phix_adapters and host contamination is performed using trimmomatic, sickle and bbmap.

2. Metagenomic Single Sample Assembly `(run_assembly.sh)`: In this step genomes from more than one species with nonuniform coverage are de novo assembled in order to characterize these metagenomes. Three assemblers (Megahit,SPAdes and IDBA) are integrated in this step. After the assembly the assembly stats are generated for user to decide which assembler worked best for their data. After the stats, assembly filter was perform to filter contigs with minimum 1000bp length.

3. Metagenomic Coassembly `(run_coassembly.sh)`: This step is similar to step 2. except that here the samples are assembled in group with Megahit and SPAdes.

4. Reference based analysis `(run_reference_analysis.sh)`: The use of reference based is bit complicated due to the fact that here we are dealing not with single genome but to the unknown number and distribution. There way to deal with this by using all the available prokaryotic genomes and align them to the reads or use marker gene approach. In this step, different state of art tools like kraken2 and metaphlan2, Diamond, Megan Humann2 are used to classify the output reads from quality control step. Diamond blastx is performed against NR database and meganized using megan databases. This output can be visualized on megan6 software for different options. Also, humann2 which is pipeline itself for taxonomic and functional Classification has also been integrated in this script.

5. Comparative analysis `(run_comparative_analysis)`: The assembled filtered contigs are then used taxonomic classification and functional annotation using kraken, diamond, megan and prokka. Prokka annotate the data by predicting genes using Prodigal  and then perform functional annotation on these genes. For homology search prokka used CDD PFAM, TIGRFAM  databases on prodigal translated protein output.


6. Coverage and binning `(run_coverage_and_bining)`: In this script spade contigs are used for further analysis due to good stats but user can change this in `run_coverage_and_bining` script in "ref" variable.
The indexed contigs sequences are mapped on its own reads to create sam and bam file using BBMAP. The generated bam file can be used in IGV/IGB viewer for checking the alignment stats. Depth file is generated from sorted bam file, which is then used for binning. The two binning tools used are Metabat and Maxbin. Binning, along with known species also attempt to recover uncultivated microbial populations which might be paying an important role in a particular sample.

7. Bin refinement `(run_binrefinement.sh)`: Different binning tools lead to different and uncertain binning results. This could be due to different algorithm behind these tools.In this step we are using two programs Binning_refiner and CheckM. Binning_refiner merges the results of different binning programs and significantly reduce the contamination level of genome bins and increase the total size of contamination-free and “good-quality” genome bins, Whereas CheckM, provides a set of tools for assessing the quality of genomes recovered from metagenomes. CheckM also provides tools for identifying genome bins that are likely candidates for merging based on marker set compatibility, similarity in genomic characteristics, and proximity within a reference genome tree.


8. Bin taxonomic classification `(run_bin_taxonomic_classification.sh)`: To assign taxonomic labels to the bins Kraken2 is used.

9. Bin functional classification `(run_bin_functional_classification.sh)`:
Bin are functionally annotated using prokka. Refer to point 6 for database detail.


## II. Amplicon Sequencing Analysis Workflow
Two major tool categories exist: (i) Operational Taxonomic Units (OTU) based and (ii) Amplicon Sequence Variant (ASV) based tools. OTU based methods cluster the reads based on a predefined identity threshold (commonly 97%) into operational taxonomic units. On the other hand, ASV based tool utilises a denoising approach on biological sequences in samples before the introduction of amplification and sequencing errors. In this review, we have included a stepwise systematic workflow for 16S rRNA using mothur and its visualization. For ASV based methods a very detailed workflow is explained for DADA2 [here](https://benjjneb.github.io/dada2/tutorial.html). In future, we aim to add additional workflows.

## Prerequisites
Although the "check_and_install" function is designed to install and setup the required software on the fly, there are some basic prerequisites that need to be satisfied:

### OS
Any Linux based distro should work. Out test OS is:

Distributor ID: Ubuntu
Description:    Ubuntu 18.04.2 LTS
Release:        18.04
Codename:       bionic

which can be got by typing 'lsb_release -a' on a Ubuntu based system.

### Hardware:
In mothur formation of the distance matrix for a large file requires

### Software and packages
Some software should be installed by the user directly as the workflow depends external software.
Without these the workflow will fail to run.

1. R  and its libraries (Phyloseq, ggplot2)
2. Python3: pip, libraries (argparse, re, Bio)

In this workflow the precompiled binaries for Mothur is used which will be automatically downloaded
