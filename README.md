# Machine Learning Quality Assessment of NGS Data

seqQscorer is a python implementation that handles quality statistics or report summaries (quality features) as input to calculate a probability of an input sample to be of low quality. This probability is calculated with pre-trained classification models. The quality features are derived from FastQ and BAM files as shown in the Figure below and described in detail in our [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/768713v2). 

The following Figure describes as well the workflow implemented to receive and preprocess NGS data from ENCODE and applying a grid-search to find the optimal classification model depending on the experimental context (species or assay) and the quality features that are provided by the user. The best Random Forest models are publicly available in this repository.

From a comprehensive feature analysis, done for this study, we derived exploratory statistics. These statistics can be used as [guidelines](statistical_guidelines/) for scientists doing quality control of NGS data. 

Furthermore, a description is provided for the application of seqQscorer to an example, also available in this repository. This example includes 8 FastQ files from two paired-end ChIP-seq experiments from ENCODE: ENCSR931HNY and ENCSR568PGX. Scripts and an explanation for the quality feature preprocessing can be found [here](example_dataset/).

<img src="figures/workflow.png" width="800">

## Applying the tool to preprocessed data

The tool runs within a Linux terminal and can be applied either to a single file or to a directory containing multiple files. These are the options of seqQscorer, followed by some examples:

Display the help text:
	
	--help		Print this help text and exit
	
These options describe the experimental setting of the data 
that is provided by the user:
	
	--spec		defines the organism of the given data
			"human" or "mouse" are available, by default it is "None"
	
	--assay		defines the assay of the given data
			"ChIP-seq", "DNase-seq", and "RNA-seq" are available, by default it is "None"
	
	--rt		defines the run-type used in the given data
			"single-ended" or "paired-ended", by default it is "None"
	
These options specify the directories containing the statistics 
or report summaries describing the data:

	--raw		path to the file or directory for the FastQC report(s)
	
	--map		path to the file or directory for the Bowtie2 mapping statistic(s)
	
	--loc		path to the file or directory for the location features (ChIPseeker)
	
	--tss		path to the file or directory for the TSS features (ChIPpeakAnno)

Additional Options:

	--out		path that specifies a file for the output of the predictions

#### Example 1: Only with the FastQC report summaries using the model for paired-ended ChIP-seq
```
python seqQscorer.py --spec human --assay ChIP-seq --rt paired-ended --raw ./example_dataset/FastQC/read1_summaries/
```

#### Example 2: Only with the Bowtie2 mapping statistics
```
python seqQscorer.py --spec human --assay ChIP-seq  --rt paired-ended --map ./example_dataset/bowtie/mapping_statistics/
```

#### Example 3: Involving all feature sets
```
python seqQscorer.py --spec human --assay ChIP-seq --rt paired-ended --raw ./example_dataset/FastQC/read1_summaries/ --map ./example_dataset/bowtie/mapping_statistics/ --loc ./example_dataset/reads_annotation/statistics/ --tss ./example_dataset/tss_annotation/statistics/
```

#### Example 4: Applying seqQscorer on only one file
```
python seqQscorer.py --spec human --assay ChIP-seq --rt paired-ended --raw ./example_dataset/FastQC/read1_summaries/ENCFF165NJF.txt --map ./example_dataset/bowtie/mapping_statistics/ENCFF165NJF.txt --loc ./example_dataset/reads_annotation/statistics/ENCFF165NJF.tsv --tss ./example_dataset/tss_annotation/statistics/ENCFF165NJF.tsv
```

#### Default Settings

It is necessary to provide at least one feature type (RAW, MAP, LOC, or TSS), however for the options regarding the experimental setting, it is not necessary to provide this information, though it is recommended. By default, the assay, species, and run-type are set to `None`. Without specifying these options the classification model is loaded that was trained on the whole dataset containing all different types of data. 

Especially when seqQscorer is applied to a dataset containing files from different assays and species, it could be of interest to use a more generic model. Even though it was not possible to evaluate the performance on data from other assays or species, it is still possible to use seqQscorer especially when using FastQC and / or the mapping statistics, as those feature sets were less varying between the assays.

##### Required packages and the version used for testing:

- python (v3.7)
- numpy (v1.17.3)
- pandas (v0.25.3)
- sklearn (v0.21.3)

Other versions might be compatible as well. In case it is not running with your python installation, the easiest is to create a conda envirnment with the following steps and running the tool within this environment. Installation commands are executed on a Linux terminal and require not more than 15 minutes.

##### Installation with ANACONDA  

First, install anaconda in case you do not have it in your linux machine. We highly recommend to install the most recent one.

```
wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh

```
Accept licence and installation requirements with "return" and "yes", but follow the instructions, you might like to change the directory for anaconda. After installation it is necessary to initialize conda with:
```
source ~/.bashrc
```

##### Create a conda environment `seqscore` with anaconda:

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n seqscore python=3.7 anaconda
```
Finally activate the environment before running the algorithm:

`conda activate seqscore`

The implementation was tested on:
- Ubuntu 16.04.6 LTS (Xenial Xerus)
- Ubuntu 18.04.3 LTS (Bionic Beaver)
- CentOS Linux 7 (Core)


