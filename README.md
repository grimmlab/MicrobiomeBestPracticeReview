# __Introduction__

<p style='text-align: justify;'> This review paper aims to provide a comprehensive workflow to perform amplicon and shotgun metagenomics analysis. There are two workflows provided. First workflow for amplicon,  using the standard mothur workflow, and along with it visualization methods for the mothur output. Second workflow for metagenomics, using a variety of tools openly available which have been stitched together to form a usable pipeline.</p>

<p style='text-align: justify;'> Both the workflows are controlled by bash scripts: amplicon_analysis.sh and metagenomics_analysis.sh. The bash scripts contain functions which call the respective underlying tools. Of-course, the tools have to exist in the system before using them, hence, a function called as "check_and_install" is written into each script which checks if the tools exists in a certain path or not.</p>

<p style='text-align: justify;'> Since the workflows utilize so many different tools, it requires quiet a bit of patience for the download and installation process. Please go through the steps below before you begin to use the workflows.</p>

![Workflow](https://github.com/dominikgrimm/MicroReviewPaper/blob/master/Figure4_README.png)

## __I.__ __Metagenomic Sequencing Analysis Workflow__

## Prerequisites
Although the `check_and_install` function is designed to install and setup the required software on the fly, there are some basic prerequisites that need to be satisfied:

### OS
Any Linux based distro should work. Out test OS is:

Distributor ID: Ubuntu <br/>
Description:    Ubuntu 18.04.2 LTS <br/>
Release:        18.04 <br/>
Codename:       bionic <br/>

'lsb_release -a' on a Ubuntu based system.


### Hardware:
<p style='text-align: justify;'> It is no secret that the hardware highly influences the speed of the workflows. The most time consuming tasks are the ones that involve assemblies, reference based alignment. A modest configuration consists of 16+cores and 100 GB of RAM with 1TB of diskspace. A majority of the diskspace is occupied by reference databases like nr database, kraken database, etc. Our HW configuration consists of 20 core CPU with 125 GB.</p>

### Software and packages
<p style='text-align: justify;'> Some software should be installed by the user directly as the workflow depends on a lot of external software.
Without these the workflow will fail to run. </p>

- gcc, g++
- java
- python2.7 and python3: pip, libraries (sys, os, shutil, stat, re, time, tarfile, operator, math, Bio, argparse)
- perl libraries (Bio)
- R
- git
- <p style='text-align: justify;'>  metabat: Install instructions can be found under https://bitbucket.org/berkeleylab/metabat/src/master/README.md. Metabat should be visible in the system PATH.</p>
- <p style='text-align: justify;'> checkM (checkm-genome): Install instructions can be found under https://github.com/Ecogenomics/CheckM/wiki/Installation.</br>
After installation the checkM database needs to be built using https://data.ace.uq.edu.au/public/CheckM_databases/ and building by using `checkm data setRoot PATH_TO_DOWNLOADED_DATABASE`</p>
**NOTE**: Make sure checkM is placed finally under `/usr/local/bin`

## Steps to run the Metagenomics workflow (metagenomics_analysis.sh)
  **1. Preparing databases:**
  ```bash
  sh prepare_databases.sh
  ```
  Insert the `LINKPATH_DB=/xxx/.../references/ to 'metagenomics_analysis.sh'`
  ```bash
  LINKPATH_DB=/xxx/.../references/
  ```
  <p style='text-align: justify;'> The databases form the core of the workflows. Unfortunately, the databases are huge and  take a long time to download and index. If these databases already exist in your system pleease modify the scripts with the correct paths. Otherwise choose the missing databases and run 'prepare_databases.sh' where the databases will be installed under the 'references' in the current directory. At the end of the preparation of databases a path will shown in the stdout which needs to be plug-in to the metabgenomics_analysis.sh script (to LINKPATH_DB variable). </p>

  The following databases are installed:
  - Human and Mouse reference genome:
  ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/Assembled_chromosomes/seq/
  ftp://ftp.ncbi.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/

  - Kraken database:
  http://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz
  (needs to be indexed with kraken2)

  - Megan database:
  http://ab.inf.uni-tuebingen.de/data/software/megan6/download/

  - NR database: Non-redundant database can be found at ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
  (needs to be index with diamond)

  - Metaphlan database:
  https://bitbucket.org/biobakery/metaphlan2/downloads/mpa_v20_m200.tar
  (needs to be built with bowtie2)

  - checkM database:
  Needs to be manually installed (Please check prerequisites).

  - Humann2 database:
  Downloaded using humann2_databases.

**2. Add source path for raw data**

  <p style='text-align: justify;'> In the 'metagenomics_analysis.sh' add the path for the raw-data. Please note that the workflow will make a local copy of the rawdata before proceeding further.</p>

  ```bash
  SRC_RAWDATA=/path_to_my_rawdata_samples/.../.../
  ```
  <p style='text-align: justify;'> **NOTE**: The sample data must always be paired end and compressed in the\*.fastq.gz format. Also the names of the pair must end with \*\_1.fastq.gz and \*\_2.fastq.gz. Example: Sample\_1.fastq.gz and Sample\_2.fastq.gz. </p>


**3. Set name of workflow**

<p style='text-align: justify;'>  Next choose an appropriate name for the analysis in the 'metagenomics_analysis.sh' script. All the sub-folders like tools, analysis, rawdata copy, etc will be created under this folder name. </p>

  ```bash
  NAME=MY_METAGENOMIC_ANALYSIS_EXP
  ```

**4. Run the workflow**

  Finally, the workflow is ready to be run
  ```bash
  sh metagenomics_analysis.sh
  ```
  There are messages on the stdout showing the status and progress of the analysis.

  <p style='text-align: justify;'> The script consists of several sub-scripts and functions. Each sub-script has its own "check_and_install". The "check_and_install" checks for the tools required to run the respective script and installs them if they are missing.</p>

  <p style='text-align: justify;'> **NOTE**: The installation of Megan is an interactive installation and requires the user to input Y/N and memory options(between ~8GB-16GB depending on sample size). We recommend to use default options. Megan will be installed in the user home directory.</p>


## Step-by-Step Analysis

<p style='text-align: justify;'>The metagenomics workflow is a time-consuming workflow. Hence, the bash scripts are kept as simple as possible. In order to perform only one type of analysis, you can always comment the remaining functions.</br>

For example, the quality control (`run_qc`) can be run only once initially and then commented for any further analysis for reruns. </br>

If the appropriate steps have already been run, then these can be commented and other steps can be run. This is of-course, not true for steps dependent on previous outputs. </p>

## Brief Description of the Each Step
  **1. run_qc.sh**</br>
  **2. run_assembly.sh**</br>
  **3. run_coassembly.sh**</br>
  **4. run_reference_analysis.sh**</br>
  **5. run_comparative_analysis.sh**</br>
  **6. run_coverage_bining.sh**</br>
  **7. run_binrefinement.sh**</br>
  **8. run_bin_taxonomic_classification.sh**</br>
  **9. run_bin_functional_classification.sh**</br>

**1. Quality control `(run_qc.sh)`**: <p style='text-align: justify;'> This scripts is running series of steps with different tools to perform quality control. [FastQC](https://github.com/s-andrews/FastQC/) is used to generate comprehensive report of data quality on raw data. Followed by this is a series of steps including removal of adapters, low quality reads, sequencing artifacts, phix adapters and host contamination is performed using trimmomatic, sickle and bbmap.  
**NOTE:** Its very important to review the QC result and change the parameters in the script based e.g read length and read quality etc. </p>

**2. Metagenomic Single Sample Assembly `(run_assembly.sh)`**: <p style='text-align: justify;'>In this step genomes from more than one species with nonuniform coverage are de novo assembled in order to characterize these metagenomes. Three assemblers (Megahit,SPAdes and IDBA) are integrated in this step. After the assembly, assembly stats is generated for user to decide which assembler worked best on their data. After the stats, assembly filter is performed to filter contigs with minimum 1000bp length. </p>

**3. Metagenomic Coassembly `(run_coassembly.sh)`**: <p style='text-align: justify;'> This step is similar to step 2. except that here the samples are assembled in group with Megahit and SPAdes. </p>

**4. Reference based analysis `(run_reference_analysis.sh)`:** <p style='text-align: justify;'>The use of reference based is bit complicated due to the fact that here we are dealing not with single genome but to the unknown number and distribution. There way to deal with this by using all the available prokaryotic genomes and align them to the reads or use marker gene approach. In this step, different state of art tools like kraken2 and metaphlan2, Diamond, Megan and Humann2 are used to classify the output reads from quality control step. Diamond blastx is performed against NR database and meganized using megan databases. This output can be visualized on megan6 software for different options. Also, humann2 which is a reference based pipeline itself for taxonomic and functional Classification has also been integrated in this script. </p>

**5. Comparative analysis `(run_comparative_analysis)`**: <p style='text-align: justify;'>The assembled filtered contigs are then used for taxonomic classification and functional annotation using kraken, diamond, megan and prokka. Prokka annotate the data by predicting genes using Prodigal and then perform functional annotation on these genes. For homology search prokka uses CDD PFAM, TIGRFAM  databases on prodigal translated protein output. </p>


**6. Coverage and binning `(run_coverage_and_bining)`**: <p style='text-align: justify;'> In this script spade contigs are used for further analysis due to good stats but user can change this in `run_coverage_and_bining` script in "ref" variable. </br>
The indexed contigs sequences are backmapped on its own reads to create sam and bam file using BBMAP. The generated bam file can be used in IGV/IGB viewer for checking the alignment statistics. Depth file is generated from sorted bam file, which is then used for binning. The two binning tools used are Metabat and Maxbin. Along with known species, binning attempt to recover uncultivated microbial populations which might be paying an important role in a particular sample.</p>

**7. Bin refinement `(run_binrefinement.sh)`**: <p style='text-align: justify;'> Different binning tools lead to different and uncertain binning results. This could be due to different algorithm behind these tools.In this step we are using two programs Binning_refiner and CheckM. Binning_refiner merges the results of different binning programs and significantly reduce the contamination level of genome bins and increase the total size of contamination-free and “good-quality” genome bins. CheckM, provides a set of tools for assessing the quality of genomes recovered from metagenomes. CheckM also provides tools for identifying genome bins that are likely candidates for merging based on marker set compatibility, similarity in genomic characteristics, and proximity within a reference genome tree.</p>

**8. Bin taxonomic classification `(run_bin_taxonomic_classification.sh)`**: To assign taxonomic labels to the bins Kraken2 is used.</p>

**9. Bin functional classification `(run_bin_functional_classification.sh)`**:
Bin are functionally annotated using prokka. Refer to point 6 for database detail.


## __II. Amplicon Sequencing Analysis Workflow__
<p style='text-align: justify;'> Two major tool categories exist: (i) Operational Taxonomic Units (OTU) based and (ii) Amplicon Sequence Variant (ASV) based tools. OTU based methods cluster the reads based on a predefined identity threshold (commonly 97%) into operational taxonomic units. On the other hand, ASV based tool utilizes a denoising approach on biological sequences in samples before the introduction of amplification and sequencing errors. In this review, we have included a stepwise systematic workflow for 16S rRNA using mothur and DADA2 and its visualization.</p>

## Prerequisites
<p style='text-align: justify;'> Although the "check_and_install" function is designed to install and setup the required software on the fly, there are some basic prerequisites that need to be satisfied:</p>

### OS
Any Linux based distro should work. Out test OS is:

Distributor ID: Ubuntu </br>
Description:    Ubuntu 18.04.2 LTS</br>
Release:        18.04</br>
Codename:       bionic</br>

type 'lsb_release -a' on a Ubuntu based system.

### Hardware:
A modest configuration consists of 4+cores and 16-32GB of RAM.

### Software and packages
Some software should be installed by the user directly as the workflow depends external software.</br>
Without these the workflow will fail to run.

- R  and its libraries (Phyloseq, ggplot2, DADA2, argparser)
- Python3: pip, libraries (argparse, re, Bio)

In this workflow the precompiled binaries for Mothur is used which will be automatically downloaded.

## Example data
This example data is a 16S rRNA gene region of gut samples collected longitudinally from a mouse post-weaning. The fastq files are generated on Miseq platform with 2X250 for V4 region.


## Steps to run the Amplicon workflow (amplicon_analysis.sh)

**1. Add source path for raw data**:
In the `amplicon_analysis.sh` add the path for the raw-data. Please note that the workflow will make a local copy of the rawdata before proceeding further.

```bash
SRC_RAWDATA=/path_to_my_rawdata_samples/.../.../
```
**NOTE**: The sample data must always be paired end and compressed in the\*.fastq.gz format.

**2. Set name of workflow**:

Next choose an appropriate name for the analysis in the 'amplicon_analysis.sh' script. All the sub-folders like tools, analysis, rawdata copy, etc will be created under this folder name.

```bash
NAME=MY_AMPLICON_ANALYSIS_EXP
```
**3. Download reference databases**:

The required reference databases will automatically be downloaded in reference folder.
The following databases are installed:

- Silva reference database:
https://mothur.org/w/images/3/32/Silva.nr_v132.tgz
https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz

- Greengenes database:
http://www.mothur.org/w/images/6/68/Gg_13_8_99.taxonomy.tgz

- RDP database
wget https://www.mothur.org/w/images/c/c3/Trainset16_022016.pds.tgz

**4. Run the workflow**:
The `amplicon_analysis.sh` consist of two main work flows which are called automatically. The example data used is from  mothur and dada2 workflow is [miseq SOP](http://www.mothur.org/wiki/MiSeq_SOP)

Finally, the workflow is ready to be run!
```bash
sh amplicon_analysis.sh
```
