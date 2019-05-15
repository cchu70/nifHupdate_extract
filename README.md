# nifHUpdate
nifH Update is a python pipeline that updates an existing databases of labeled nifH gene sequences. 
Eventually, this pipeline will be generalized such that any database can be updated. 

nifHUpdate requires a configuration file with the file paths to relevant files. Dependencies and the configuration requirements are explained in detail below.

## minimap Method
The edirect method had several shortcomings, namely the time it took to retrieve relevant sequences from NCBI Nucleotide database, and it's reliance on proper annotation of the queried term. Instead, this pipeline takes fasta files of whole genome sequences from some existing database. To extract relevant genomes sequences, we use minimap to perform a preliminary alignment with your existing database. 

# Tutorial
The following steps will guide you through how to set up and run this pipeline

1) Activate your environment if desired. I used ddocent_env (http://www.ddocent.com/bioconda/), which contains many of the tools necessary.
```
source activate ddocent_env
```
2) Download additional dependencies if necessary (requirements.txt). 
```
pip install cd-hit-auxtools
```
3) Make relevant file-of-file-names for your local database you are searching.

Here is an example of a fofn file `Klebsiella_genomes.fofn`
```
data/home/user/nuccore/570_1.fasta
data/home/user/nuccore/570_2.fasta
data/home/user/nuccore/570_3.fasta
```
4) Create Configuration File. Below are each label, and what they are used for. Those with a * are mandatory. The others have default values. 
- PREFIX * will be used to name all intermediate and output files and directories.
- DBFILE * file path to your existing database in fasta format.
- NUCCORE * file path to a text file containing file paths to all relevant (unzipped) fasta files of your local database you are searching from, even if it is a single fasta file
- MIN_MINIMAP_ALIGNLEN Default = 200 bp. Minimap will drop any alignments that are below this length.
- MIN_BLASTN_ALIGNLEN Default = 200 bp. Blastn will drop any alignments that are below this length.
- PIDENT_CUTOFF Default = 75%, which is the cutoff used for degree of similarity among organisms from the same family. This is used to capture homologs. 
- '#' Used to denote a comment. 

These labels can be specified in any order. The program will halt if there are any formatting errors or unnecessary labels. Here is an example configuration file called `my_config.txt`
```
# May 4, 2019 klebsiella genome test run
PREFIX  AAA
DBFILE    /Users/nifH_extract/oldDatabase.fasta
NUCCORE    /Users/nifH_extract/Klebsiella_genomes.fofn
MIN_MINIMAP_ALIGNLEN 300
```

4) Run the following command in your desired directory.

```
$ path/to/nifHUpdate_Lib/nifHUpdate_caller.py my_config.txt my_log_file_name.txt
```
`my_config.txt` is the configuration file you made. `my_log_file_name.txt` does not have to already exist. It will be created in your current working directory.

You can restart the pipeline at any point adding the `-s` option and specifying what stage to start on. For example:

```
$ path/to/nifHUpdate_Lib/nifHUpdate_caller.py my_config.txt my_log_file_name.txt -s blastn
```
View what stages are available below. Running the command without the `-s` option will default to rerunning the pipeline from the beginning. It is recommended that if you want to save the result of each run, instead of rerunning with the same prefix, make a new configuration file with a different prefix label. 

### Stages
- minimap
  This stage does the preliminary mapping of your old database onto your subject database. It outputs `.paf` files, which provides information about the alignments
- minimap filter
  This stage parses the `.paf` and pulls out the relevant fasta records from your subject database. This reduced set is then passed to the next stage.
- blastn
  Blastn does a more thorough alignment to pull out the exact sequences with a high enough similarity (specified in your configuration file). This produces the results in output format 6, with the addition of other details used to extract the relevant sequences. 
- filter_best_alignments
  This extracts the aligned regions that are above predetermined thresholds and puts them into new fasta files. It also formats the fasta file ID and description to include which cluster that sequence aligned to. 
- cluster
  From the renamed records, this stage simply sorts all the sequences that aligned to the same cluster in its own fasta file.
- deduplicate
  Removes copies of the same sequence in each of the cluster files

### Packages
nifHUpdate requires the following packages, and is highly recommended to run in an environment
- Python 
- Biopython
- cd-hit-auxtools 
- blastn


