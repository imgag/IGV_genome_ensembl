help:
	@cat Makefile
	
Homo_sapiens.GRCh37.87.gff3:
	wget ftp://ftp.ensembl.org/pub/grch37/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
	gunzip Homo_sapiens.GRCh37.87.gff3.gz

hgnc_complete_set.txt:
	wget -O - ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc_complete_set.txt.gz | gunzip > hgnc_complete_set.txt
	
1kg_v37.genome:	
	wget http://igv.broadinstitute.org/genomes/1kg_v37.genome
	
convert: Homo_sapiens.GRCh37.87.gff3 hgnc_complete_set.txt 1kg_v37.genome
	python gff_to_genepred_converter.py Homo_sapiens.GRCh37.87.gff3 hgnc_complete_set.txt 1kg_v37.genome 1kg_v37_ensembl.genome

all: convert