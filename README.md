# nifHUpdate
nifH Update is a python pipeline that updates existing databases of gene sequences. 
It is a generalized such that any database can be updated using either of the two methods. 

## edirect Method
This uses Entrez-direct (or edirect) to search NCBI databases using key terms. 
### Packages
nifHUpdate requires the following packages, and is highly recommended to run in an environment
- Python 
- Bio
- cd-hit-auxtools
ddocent_env contains many of the tools necessary. I used ddocent_env (http://www.ddocent.com//bioconda/)
