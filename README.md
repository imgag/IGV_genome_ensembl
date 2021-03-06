# IGV_genome_ensembl
A IGV reference genome annotated with Ensembl gene/transcript info

The tool takes a gff3 file with Ensembl annotations and converts it into a genePred file. Then it uses the HGNC ids in the gff3 file to annotate the genes/transcripts with the correct names (from the HGNC file). After that the genePred file is modified to fit the requirements of IGV. In the last step the gene file in the reference genome file is replaced.  

To convert the gff3 file to genePred the tool `gff3ToGenePred` from http://hgdownload.soe.ucsc.edu/admin/exe/ is used (which is automatically downloaded on runtime).

## Requirements
- Python 3
- Linux x64

## Usage
To run the python script itself:
```
python gff_to_genepred_converter.py [-h] gff_file hgnc_file genome_file output
```  
To get extended help:  
```
python gff_to_genepred_converter.py -h
``` 
 
To download all required data and run the tool:
```
make convert_GRCh37
``` 
or:
```
make convert_GRCh38
``` 
This generates a .genome file called `GRCh37_ensembl.genome`/`GRCh38_ensembl.genome` which can be imported into IGV via `Genomes` -> `Load Genome from file...`

