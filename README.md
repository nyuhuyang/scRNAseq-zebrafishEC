# scRNAseq-zebrafishEC
single cell RNA-Seq analysis of zebrafish endothelial cells. Data from Todd R. Evans lab (tre2003@med.cornell.edu) and postdoc Yahui Lan (yal2011@med.cornell.edu)

Creating a Reference Package with cellranger mkref
ftp://ftp.ensembl.org/pub/release-91/fasta/danio_rerio/cdna/
ftp://ftp.ensembl.org/pub/release-91/gtf/danio_rerio/

cellranger mkgtf Danio_rerio.GRCz10.91.gtf Danio_rerio.GRCz10.91.filtered.gtf --attribute=gene_biotype:protein_coding

cellranger mkref --genome=Danio_rerio.GRCz10 --fasta=Danio_rerio.GRCz10.cdna.all.fa --genes=Danio_rerio.GRCz10.91.filtered.gtf

cellranger mkref --genome=Danio_rerio.GRCz10 --fasta=Danio_rerio.GRCz10.cdna.all.fa --genes=Danio_rerio.GRCz10.91.filtered.gtf


cellranger mkref --genome=Danio_rerio.GRCz10 --fasta=Danio_rerio.GRCz10.cdna.all.fa --genes=Danio_rerio.GRCz10.91.chr.gtf
