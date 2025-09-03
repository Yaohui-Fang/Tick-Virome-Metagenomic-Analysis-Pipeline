#! /bin/bash
#############################################################
# Title: Pipeline of Metagenomic Sequencing
# Author: Yao-hui Fang
# Date: 3/16/2024
# Version: 1.3
# Description: QC,  Filtration, Assembly, Quant and Annotation
# Requre file list: NGS data: *.fq.gz
#############################################################

##1.QC##
rm -rf 1.QC && mkdir 1.QC && cd 1.QC

fastqc ../*.fq.gz

cd ..

##2.Filtration##
rm -rf 2.Filtration && mkdir 2.Filtration && cd 2.Filtration

for i in `ls ../*.R1.fq.gz`;do bt=$(basename $i .R1.fq.gz); bowtie2 -x /data/Tickdb/Tickdb -1 ../${bt}.R1.fq.gz -2 ../${bt}.R2.fq.gz -S ${bt}.sam;done

for i in `ls *.sam`;do bt=$(basename ${i} .sam);awk '{if($4==0)print;}' ${bt}.sam | cut -f 1 | uniq > ${bt}.txt;done

for i in `ls *.sam`;do bt=$(basename ${i} .sam);seqtk subseq ../${bt}.R1.fq.gz ${bt}.txt > ${bt}.R1.hostfree.fq;seqtk subseq ${bt}.R2.fq.gz ${bt}.txt > ${bt}.R2.hostfree.fq;done

rm *sam

cd ..

##3.Assembly##
rm -rf 3.Assembly && mkdir 3.Assembly && cd 3.Assembly

for i in `ls ../2.Filtration/*.R1.hostfree.fq`;do bt=$(basename ${i} .R1.hostfree.fq);Trinity --seqType fq --max_memory 300G --left ../2.Filtration/${bt}.R1.hostfree.fq  --right ../2.Filtration/${bt}.R2.hostfree.fq --CPU 48 --output ${bt}_trinity_out;mv ${bt}_trinity_out/Trinity.fasta ${bt}_Trinity.fasta;done

rm -rf *trinity_out/

cd ..

##4.Quant##
rm -rf 4.Quant && mkdir 4.Quant && cd 4.Quant

for i in `ls ../3.Assembly/*_Trinity.fasta`;do bt=$(basename ${i} _Trinity.fasta);align_and_estimate_abundance.pl --transcripts ../3.Assembly/${bt}_Trinity.fasta --seqType fq --left ../2.Filtration/${bt}.R1.hostfree.fq --right ../2.Filtration/${bt}.R2.hostfree.fq --est_method RSEM --aln_method bowtie2 --trinity_mode --prep_reference --sample_name ${bt} --output_dir rsem_quant --thread_count 48;done

cd ..

##5.Annotation##
rm -rf 5.Annotation && mkdir 5.Annotation && cd 5.Annotation

for i in `ls ../3.Assembly/*_Trinity.fasta`;do bt=$(basename ${i} _Trinity.fasta);blastn -query ../3.Assembly/${bt}_Trinity.fasta -db /data/db/blastdb/nt -out ${bt}_blastn_results.tsv -outfmt "6 qseqid sseqid pident length evalue bitscore stitle staxids" -evalue 0.00001 -num_threads 48 -max_target_seqs 1;done

for i in `ls ../3.Assembly/*_Trinity.fasta`;do bt=$(basename ${i} _Trinity.fasta);diamond blastx -q ../3.Assembly/${bt}_Trinity.fasta -d /data/db/blastdb/nr.dmnd -o ${bt}_blastx_results.tsv -e 0.00001 -k 1 -p 48 --outfmt "6 qseqid sseqid pident length evalue bitscore stitle staxids";done

echo -e "qseqid\tsseqid_blastn\tpident_blastn\tlength_blastn\tevalue_blastn\tbitscore_blastn\tstitle_blastn\tstaxids_blastn\tsseqid_blastx\tpident_blastx\tlength_blastx\tevalue_blastx\tbitscore_blastx\tstitle_blastx\tstaxids_blastx\tTPM" > header.txt

for i in `ls *_blastn_results.tsv`;do bt=$(basename ${i} _blastn_results.tsv);

awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8}' ${bt}_blastn_results.tsv | sort -k1,1 > ${bt}_blastn_clean.tsv;

awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8}' ${bt}_blastx_results.tsv | sort -k1,1 > ${bt}_blastx_clean.tsv;

awk 'BEGIN{OFS="\t"} NR > 1 {print $1, $6}' ../4.Quant/${bt}.genes.results | sort -k1,1 > ${bt}_tpm_clean.tsv;

join -t $'\t' -a1 -a2 -e "NA" -o auto ${bt}_blastn_clean.tsv ${bt}_blastx_clean.tsv > ${bt}_blastn_blastx_merged.tsv;

join -t $'\t' -a1 -e "NA" -o auto ${bt}_blastn_blastx_merged.tsv ${bt}_tpm_clean.tsv > ${bt}_merged_final.tsv;

cat header.txt ${bt}_merged_final.tsv > ${bt}_merged_annotation_TPM.tsv;

done

rm *blastn_clean.tsv *blastx_clean.tsv *tpm_clean.tsv *blastn_blastx_merged.tsv *header.txt *merged_final.tsv

cd ..

##6.Virus##
rm -rf 6.Virus && mkdir 6.Virus && cd 6.Virus

cut -f8,15 ../5.Annotation/merged_annotation_TPM.tsv | tail -n +2 | tr '\t' '\n' | grep -v '^$' | sort | uniq > all_staxids.txt

taxonkit lineage all_staxids.txt > staxid_lineage.txt

grep 'Viruses' staxid_lineage.txt | cut -f1 > virus_taxids.txt

awk 'BEGIN{FS=OFS="\t"} NR==FNR{ v[$1]; next } FNR==1{ print; next } ($8 in v || $15 in v){ print }' virus_taxids.txt ../5.Annotation/merged_annotation_TPM.tsv > virus_annotation_TPM.tsv

rm all_staxids.txt staxid_lineage.txt virus_taxids.txt

cd ..
