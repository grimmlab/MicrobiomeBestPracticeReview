## Introduction

This review paper aims to provide a comprehensive workflow to perform amplicon and shotgun metagenomics analysis. There are two workflows provided. First workflow for amplicon,  using the standard mothur workflow, and along with it visualization methods for the mothur output. Second workflow for metagenomics, using a variety of tools openly available which have been stitched together to form a usable pipeline.

Both the workflows are controlled by bash scripts: amplicon_analysis.sh and metagenomics_analysis.sh. The bash scripts contain functions which call the respective underlying tools. Of-course, the tools have to exist in the system before using them, hence, a function called as "check_and_install" is written into each script which checks if the tools exists in a certain path or not.

Since the workflows utilize so many different tools, it requires quiet a bit of patience for the download and installation process. Please go through the steps below before you begin to use the workflows.

## I Metagenomic Sequencing Analysis Workflow

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
3. python2.7 and python3: pip, libraries (sys, os, shutil, stat, re, time, tarfile, operator, math, )
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

**NOTE**: The installation of Megan is an interactive installation and requires the user to input Y/N and memory options. We recommend to use default options. Megan will be installed in the user home directory.


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

1. Quality control `(run_qc.sh)`: This scripts is running series of steps with different tools to perform quality control. [FastQC](https://github.com/s-andrews/FastQC/) is used to generate comprehensive report of data quality on raw data. Followed by this is a series of steps including removal of adapters, low quality reads, sequencing_artifacts, phix_adapters and host contamination is performed using trimmomatic, sickle and bbmap.

2. Metagenomic Assembly `(run_assembly.sh)`: The


## II Amplicon Sequencing Analysis Workflow
