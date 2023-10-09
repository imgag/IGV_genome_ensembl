CURDATE=$(shell date +"%Y-%m-%d")

help:
	@cat Makefile
	

hgnc_complete_set.tsv:
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt > hgnc_complete_set.tsv

create_json_GRCh38: hgnc_complete_set.tsv
	mkdir -p $(CURDATE)_GRCh38
	python3 generate_igv_genome.py GRCh38_template.json hgnc_complete_set.tsv $(CURDATE)_GRCh38/GRCh38_ensembl.json
	
	
### For old IGV genome format
Homo_sapiens.GRCh37.87.gff3:
	wget ftp://ftp.ensembl.org/pub/grch37/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
	gunzip Homo_sapiens.GRCh37.87.gff3.gz

Homo_sapiens.GRCh38.110.gff3:
	wget https://ftp.ensembl.org/pub/release-110/gff3/homo_sapiens/Homo_sapiens.GRCh38.110.gff3.gz
	gunzip Homo_sapiens.GRCh38.110.gff3.gz

1kg_v37.genome:	
	wget http://igv.broadinstitute.org/genomes/1kg_v37.genome

hg38.genome:
	wget https://s3.amazonaws.com/igv.org.genomes/hg38/hg38.genome

convert_GRCh37: Homo_sapiens.GRCh37.87.gff3 hgnc_complete_set.tsv 1kg_v37.genome
	python gff_to_genepred_converter.py Homo_sapiens.GRCh37.87.gff3 hgnc_complete_set.tsv 1kg_v37.genome GRCh37_ensembl.genome

convert_GRCh38: Homo_sapiens.GRCh38.104.gff3 hgnc_complete_set.tsv hg38.genome
	python gff_to_genepred_converter.py Homo_sapiens.GRCh38.110.gff3 hgnc_complete_set.tsv hg38.genome GRCh38_ensembl.genome


