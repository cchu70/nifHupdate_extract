# nifHUpdate
nifH Update is a python pipeline that updates an existing databases of labeled nifH gene sequences. 
Eventually, this pipeline will be generalized such that any database can be updated using either of the following two methods. 

Both methods requires a configuration file with the file paths to relevant files. Dependencies and the configuration requirements are explained in detail below.

# Sample run
The following steps will guide you through how to set up and run this pipeline

1) Activate your environment
2) Download dependencies (requirements.txt). 
3) Create Configuration File (details below)
4) Run the following command in your desired directory.

```
$ path/to/nifHUpdate_Lib/nifHUpdate_caller.py my_config.txt my_log_file_name.txt
```
`my_config.txt` is the condiguration file you made. `my_log_file_name.txt` does not have to already exist. It will be created in your current working directory. If your configuration file and environment are set up properly, the pipeline will output your new database. 

You can restart the pipeline at any point adding the `-s` option and specifying what stage to start on. For example:

```
$ path/to/nifHUpdate_Lib/nifHUpdate_caller.py my_config.txt my_log_file_name.txt -s blastn
```
View what stages are available below.

## edirect Method
This uses Entrez-direct (or edirect) to search NCBI databases using a query term.
### Stages
- esearch
- efetch
- blastn
- filter_best_alignments
- cluster
- deduplicate

### Packages
nifHUpdate requires the following packages, and is highly recommended to run in an environment
- Python 
- Biopython
- edirect (installation instructions here: https://www.ncbi.nlm.nih.gov/books/NBK179288/)
- cd-hit-auxtools 
- blastn
I used ddocent_env (http://www.ddocent.com//bioconda/), which contains many of the tools necessary.

### Configuration File
The following describes the purpose of each label in the config file
- PREFIX will be used to name all intermediate and output files and directories
- DBFILE file path to your existing database in fasta format
- DBNAME name of the database files created for blastn
- 

Here is an example Configuration File

## minimap Method
The edirect method had several shortcomings, namely the time it took to retrieve relevant sequences from NCBI Nucleotide database, and it's reliance on proper annotation of the queried term. 
