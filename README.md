# Note: This project is discontinued! See the [GSvar documentation](https://github.com/imgag/ngs-bits/blob/master/doc/GSvar/igv_integration.md#igv-gene-trackgenome) for details how we build a IGV genome based on Ensembl now.

# IGV genome based on Ensembl
A IGV reference genome annotated with Ensembl gene/transcript info.

## Requirements
- Python 3
- Linux x64
- `bgzip` and `tabix` for JSON format

## New JSON genome format
Since IGV version 2.11.0 a new genome file format based on JSON is supported (https://github.com/igvteam/igv/wiki/JSON-Genome-Format).
The script `generate_igv_genome.py` creates a genome file in this format. Additional it downloads all linked data sources to the local storage and updates the gene names of included gff3 with the current HGNC symbol.

### Usage
To run the python script itself:
```
python3 generate_igv_genome.py template.json hgnc_complete_set.tsv output_genome.json
```  
To get extended help:  
```
python3 generate_igv_genome.py -h
``` 
To generate a genome file for GRCh38 with ensembl gene annotation:
```
make create_json_GRCh38
```

## Old .genome format
The tool takes a gff3 file with Ensembl annotations and converts it into a genePred file. Then it uses the HGNC ids in the gff3 file to annotate the genes/transcripts with the correct names (from the HGNC file). After that the genePred file is modified to fit the requirements of IGV. In the last step the gene file in the reference genome file is replaced.  

To convert the gff3 file to genePred the tool `gff3ToGenePred` from http://hgdownload.soe.ucsc.edu/admin/exe/ is used (which is automatically downloaded on runtime).


### Usage
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

