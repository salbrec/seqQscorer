# Machine Learning Quality Assessment of NGS Data

seqQscorer is a python implementation that handles quality statistics or report summaries (quality features) as input to calculate a probability of an input NGS sample to be of low quality. This probability is calculated with pre-trained classification models. The quality features are derived from FastQ and BAM files as shown in the Figure below and described in detail in our [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/768713v2). 

The following Figure describes the workflow implemented to receive and preprocess NGS data from ENCODE and applying a grid-search to find the optimal classification model. The optimization depends on the experimental context (species or assay) and the quality features that are provided by the user. Already computed classification models are not available in the github. However, the software conatins settings for an over all well-performing genereic model and multiple more specialized model, that can be trained with the ENCODE data or new data. The first time a model is needed, it is trained and serialized into the folder `models`. Afterwards it is not necessary to compute it again. Note, seqQscorer trains models on the preprocessed ENCODE data (in utils) as described in the article.

Additionally the script `deriveFeatureSets.py` is available in this repository. It allows the user to preprocess all quality feature sets that are needed to properly run seqQscorer and provides the results in an already readable way for seqQscorer. An example for its usage is provided below. 

<img src="figures/workflow.png" width="850">

## Software installation

Especially the preprocessing requires several bioinformatic tools and software packages. The easiest and fastest way to get ready for seqQscorer is pulling the docker and running the scripts inside the docker. The following descriptions explain how to get started with docker. However, it is also possible to install everything manually. For both it is of course necessary to clone this repository:

```
git clone https://github.com/salbrec/seqQscorer.git
```

We also listed all commands (bottom of this README) that were ran starting from an Anaconda installation on a Linux OS. 

However, for the start with docker, please open a Linux terminal and run the following commands to first install docker, then pulling the image, and finally running the image.

```
sudo apt-get install docker
sudo apt-get install docker.io 

sudo docker login

sudo docker pull salbrec/seqqdocker
sudo docker run -i -t -v "/home/:/home/" salbrec/seqqdocker /bin/bash
```
#### Check out further installation guides for running Docker with Windows 10 or creating your own conda environment (on the bottom of this README)

## Getting the first information from the software within your installation

Change directory into the seqQscorer repository and run the usage information:

```
python deriveFeatureSets.py --help
python seqQscorer.py --help
```

The expected output looks like this, for `python deriveFeatureSets.py --help`:

```

usage: deriveFeatureSets.py [-h] --fastq1 FASTQ1 [--fastq2 FASTQ2] --btidx
                            BTIDX [--outdir OUTDIR] [--cores CORES]
                            [--fastqc {1,2}] [--assembly {GRCh38,GRCm38}]
                            [--gtf GTF]

seqQscorer Preprocessing - derives feature sets needed leveraged by seqQscorer

optional arguments:
  -h, --help            show this help message and exit
  --fastq1 FASTQ1, -f1 FASTQ1
                        Input fastq file. Either the fastq file for a single-
                        end sample or the fastq file for read 1 of a paired-
                        end sample.
  --fastq2 FASTQ2, -f2 FASTQ2
                        In case of a paired-end sample, the fastq file for
                        read 2.
  --btidx BTIDX, -ix BTIDX
                        Filename prefix for Bowtie2 Genome Index (minus
                        trailing .X.bt2).
  --outdir OUTDIR, -o OUTDIR
                        Output directory. Default: "./feature_sets/"
  --cores CORES, -c CORES
                        Defines the number of processors (CPUs) to be used by
                        bowtie2 and samtools. (decreases runtime)
  --fastqc {1,2}, -f {1,2}
                        The fastq on which FastQC is applied on. Can
                        optionally be selected for paired-end samples.
  --assembly {GRCh38,GRCm38}, -a {GRCh38,GRCm38}
                        Species assembly needed to define the gene structure /
                        annotation used by the bioconductor functions. (has to
                        be consistent with the species used in for Bowtie2)
  --gtf GTF, -g GTF     File path for a gtf file to be used to get the LOC and
                        TSS features. (--assembly will be ignored then)



```

Expected output for `python seqQscorer.py --help`:

```

usage: seqQscorer.py [-h] --indir INDIR [--species {generic,human,mouse}]
                     [--assay {generic,ChIP-seq,DNase-seq,RNA-seq}]
                     [--runtype {generic,single-ended,paired-ended}]
                     [--model MODEL] [--noRAW] [--noMAP] [--noLOC] [--noTSS]
                     [--noFS] [--bestCalib] [--peaktype {narrow,broad}]
                     [--probOut PROBOUT] [--compOut COMPOUT]
                     [--inputOut INPUTOUT]

seqQscorer - A machine learning application for quality assessment of NGS data

optional arguments:
  -h, --help            show this help message and exit
  --indir INDIR, -i INDIR
                        Input directory containing the feature set files
  --species {generic,human,mouse}, -s {generic,human,mouse}
                        Species specifying the model used.
  --assay {generic,ChIP-seq,DNase-seq,RNA-seq}, -a {generic,ChIP-seq,DNase-seq,RNA-seq}
                        Assay specifying the model used.
  --runtype {generic,single-ended,paired-ended}, -r {generic,single-ended,paired-ended}
                        Run-Type specifying the model used.
  --model MODEL, -m MODEL
                        Path to a serialized model, trained on own data. If
                        used, the parameters --species, --assay, and --runtype
                        have no impact on the classification model.
  --noRAW               Ignore all RAW features.
  --noMAP               Ignore all MAP features.
  --noLOC               Ignore all LOC features.
  --noTSS               Ignore all TSS features.
  --noFS                Switch off feature selection. (has only an impact if
                        the best performance was achieved with chi2 or RFE)
  --bestCalib           Classifier setting is used that achieved the lowest
                        brier score, hence the best calibration of the
                        probabilities.
  --peaktype {narrow,broad}, -pt {narrow,broad}
                        Optionally specify the peak-type for ChIP-seq data.
  --probOut PROBOUT, -po PROBOUT
                        To specify an output file for the probabilities.
                        Output will be tab-separated.
  --compOut COMPOUT, -co COMPOUT
                        To specify an out file for the comprehensive output.
                        Output will be kind of tab-separated.
  --inputOut INPUTOUT, -io INPUTOUT
                        To specify an out file that will contain the parsed
                        input. Output will be tab-separated.


```


## Preprocessing for fastq files

You'll need a bowtie2 genome index as input for seqQscorer. If you don't have the genome index, check out the README in [genome index](utils/genome_index/) for a reference for the download from the official Bowtie2 webpage.
If you would like to use a small example file right away within the docker, there are one here: `/var/examples/single/ENCFF165NJF.fastq.gz`. Even smaller examples are here for the paired-end test: `/var/examples/paired/`. Note, these are examples just for testing the installation, he files were reduced randomly and especially the paired-end example has only ~100k reads.

All seqQscorer feature sets can be derived by using the provided python script applied on an input fastq file or a pair of fastq files in case of paired-end sequencing. 

```
python deriveFeatureSets.py --fastq1 /var/examples/single/ENCFF165NJF.fastq.gz --btidx ./utils/genome_index/GRCh38_noalt_as/GRCh38_noalt_as --assembly GRCh38
```
The results will be in the default output folder `./feature_sets/`, use `--outdir` to specify the destination of the feature sets.

The following run represents a paired-end example using the genome index for *Mus musculus*. The parameter `--cores` allows the usage of multiple CPUs to accelarate especially the mapping.

```
python deriveFeatureSets.py --fastq1 /var/examples/paired/ENCFF310LVJ.fastq.gz --fastq2 /var/examples/paired/ENCFF410LTA_r2.fastq.gz
--cores 4 --btidx ./utils/genome_index/mm10/mm10 --assembly GRCm38 --outdir ./mouse_pe/
```

### Preprocessing with own gene structure files

The preprocessing procedure runs for the genome assemblies GRCh38 and GRCm38 which were used within the study. In case you would like to run everything with your own data for another species or an older human or mouse genome assembly, it is possible to use gtf files for the Bioconductor packages that derive the LOC and TSS features. 

Since we ran into compatibility problems while implementing this option, we recommend to use an **NCBI** genome index for Bowtie2, downloaded from this website: [http://bowtie-bio.sourceforge.net/bowtie2/index.shtml](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). 

The gtf file should be downloaded from Ensembl via their FTP Download website: [https://www.ensembl.org/info/data/ftp/index.html](https://www.ensembl.org/info/data/ftp/index.html).

With the index and gtf present, it is straight forward to preprocess fastq files for other organisms. An example for *Rattus norvegicus*:

```
python deriveFeatureSets.py --fastq1 /var/examples/ENCFF165NJF.fastq.gz --btidx ./utils/genome_index/Rnor_6.0/Rnor_6.0 --outdir ./rat_data/ --gtf ./utils/gene_structure/Rattus_norvegicus.Rnor_6.0.101.gtf -c 4
```

## Applying seqQscorer on preprocessed data

After deriving the feature sets the application of seqQscorer can be as simple as the following line. Check out the parameters that allow you to specify the feature sets used and especially the classification model that is applied.

```
python seqQscorer.py --indir ./feature_set_examples/
```

By default the generic classification model is (trained and) used to calculate the quality probabilities. According to our analyses, its performance is comparable to the more specialized models. However, depending on the data that was available for our investigations, some specialized model are available. You can specify the model with the paramters `--species`, `--assay`, and `--runtype`. seqQscorer will then automatically select the model that achieved the highest auROC (area under ROC curve). The specialized models, available for all feature set combinations out of RAW, MAP, LOC, and TSS are: (human, ChIP-seq, single-ended), (mouse, ChIP-seq, single-ended), (human, ChIP-seq, paired-ended), (human, DNase-seq, paired-ended), (mouse, DNase-seq, paired-ended), (human, RNA-seq, single-ended).

By default, seqQscorer uses all feature sets, but there are also parameters to be used if certain feature sets should be ignored.

seqQscorer prints some interesting information to the console. With options such as `--probOut` and `--compOut` it is possible to save your results to a given file name.

From our preprocessed grid search, we already defined algorithms and parameter settings that perfromed well on the different datasets. Some of them achieved the best performance when applied together with a feature selection method. However, if you want to switch off the feature selection, use `--noFS`. 

By default the model is selected that achieved the highest predictive performance, thus the highest auROC. You can also tell seqQscorer to use the model that achieved the best calibration with respect to the probabilities. This is done by `--bestCalib`. Note that sometimes the calibration (expressed by the Brier-loss) could not be drastically improved by another model that differs from the model that achieved the highest auROC.

## Further installation guides

### Installation with Docker Desktop on Windows

To install Docker Desktop follow the instructions on their website:
https://docs.docker.com/docker-for-windows/install/

Use git from powershell to clone seqQscorer
```
git clone https://github.com/salbrec/seqQscorer.git

```
To get the image and activate it is similar to Linux. 
However it is advisable to only link the SeqQscorer folder.
Docker mentions, that binding Windows Volumes can lead to performance drops and suggest to bind mounted folders from the linux filesystem in wsl rather than a folder on the windows filesystem.
Both works fine and can be accessed via powershell from the windows side or from the bash from the Linux/WSL side.

Below is an example from powershell, for linux just add sudo in front.
```
docker pull salbrec/testing
docker run -i -t --name seqQscorer -v "C:/Users/User/seqQscorer:/seqQscorer" salbrec/seqqdocker 
```
Now you can just change to the newqscorer folder and start usign the software!
```
(SeqQscorer) root@ xxx : cd seqQscorer
```
In this example the SeqQscorer folder that is on windows is to find in the root of the docker image.
The docker image is named SeqQscorer and can be invoked by this name in the future.
You can copy files from the windows side (like your fastq's) and compute from the docker side.

Docker advises to use WSL, the mounted Linux System for Windows. If you want to use this, the installation and handling would be similar to the normal Linux installation, but the Installation of Docker Desktop for windows also needs to be done.

### Installation with ANACONDA  

First, install anaconda in case you do not have it in your linux machine. We recommend to use the one that is suggested here. For the installation of Anaconda run the following two lines in your terminal.

```
wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
bash Anaconda3-2020.11-Linux-x86_64.sh

```
Accept licence and installation requirements with "return" and "yes", but follow the instructions, you might like to change the directory for anaconda. After installation it is necessary to initialize conda with:
```
source ~/.bashrc
```

You can use our yml file `conda_env.yml` to create the conda environment with this call: `conda env create -f conda_env.yml`.
Afterwards the R packages are installed by running thses lines within R.

```
# Within R run the following lines to install the R packages needed
install.packages("BiocManager")
BiocManager::install("ChIPpeakAnno")
BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
```



