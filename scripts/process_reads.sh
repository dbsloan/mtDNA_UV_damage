#!/bin/sh

for file in *.fastq
do

cutadapt -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG --discard-untrimmed -o ${file%fastq}trimmed.fastq $file

bowtie2 --phred33 -x sacCer3_genome.fasta -U ${file%fastq}trimmed.fastq -S ${file%fastq}trimmed.sam >>${file%fastq}bowtie.log.txt 2>>${file%fastq}bowtie.err.txt
done 


for file in *.trimmed.sam
do
python3 custom_dedup.py $file  ${file%sam}dedup.sam
done




for file in *.sam
do

samtools view -h -q30 $file >  ${file%sam}unique.sam

grep -E '^@|chrM' ${file%sam}unique.sam > ${file%sam}unique.mt.sam
grep -v -E 'chrM' ${file%sam}unique.sam > ${file%sam}unique.nuc.sam
done


for file in *trimmed.unique.mt.sam
do
python3 xr_seq_region_sum_RPKM_211006.py  $file ${file%mt.sam}nuc.sam sacCer3_chrM.coords.txt chrM ${file%sam}coordlike.txt ${file%sam}nuc_vs_org_RPKM.txt ${file%sam}region_cov_sum_RPKM.txt ${file%sam}gene_cov_sum.txt ${file%.trimmed.unique.mt.sam} ${file%sam}filter.sam 12071326
done


for file in *trimmed.dedup.unique.mt.sam
do
python3 xr_seq_region_sum_RPKM_211006.py  $file ${file%mt.sam}nuc.sam sacCer3_chrM.coords.txt chrM ${file%sam}coordlike.txt ${file%sam}nuc_vs_org_RPKM.txt ${file%sam}region_cov_sum_RPKM.txt ${file%sam}gene_cov_sum.txt ${file%.trimmed.dedup.unique.mt.sam} ${file%sam}filter.sam 12071326
done


for file in *trimmed.unique.mt.filter.sam
do
python3 xr_seq_read_nuc_comp.210816.py $file chrM ${file%.trimmed.unique.mt.filter.sam} ${file%sam}nuc_comp.txt ${file%sam}di_py_freq.txt

done

for file in *trimmed.dedup.unique.mt.filter.sam
do
python3 xr_seq_read_nuc_comp.210816.py $file chrM ${file%.trimmed.dedup.unique.mt.filter.sam} ${file%sam}nuc_comp.txt ${file%sam}di_py_freq.txt

done

for file in *trimmed.dedup.unique.nuc.sam
do
samtools view $file  |cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort -n| uniq -c > ${file%sam}length.txt
c=${file%.trimmed.dedup.unique*sam}
awk -v c="$c" '{print $0 "\t" c }' ${file%sam}length.txt > ${file%sam}length.name.txt
done

cat *trimmed.dedup.unique.nuc.length.name.txt > dedup.nuc.length.dist.txt




for file in *trimmed.dedup.unique.mt.filter.sam
do
samtools view $file  |cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort -n| uniq -c > ${file%sam}length.txt
c=${file%.trimmed.dedup.unique*sam}
awk -v c="$c" '{print $0 "\t" c }' ${file%sam}length.txt > ${file%sam}length.name.txt
done

cat *trimmed.dedup.unique.mt.filter.length.name.txt > dedup.mt.length.dist.txt



for file in *trimmed.unique.nuc.sam
do
samtools view $file  |cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort -n| uniq -c > ${file%sam}length.txt
c=${file%.trimmed.unique*sam}
awk -v c="$c" '{print $0 "\t" c }' ${file%sam}length.txt > ${file%sam}length.name.txt
done
cat *trimmed.unique.nuc.length.name.txt > nuc.length.dist.txt


for file in *trimmed.unique.mt.filter.sam
do
samtools view $file  |cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort -n| uniq -c > ${file%sam}length.txt
c=${file%.trimmed.unique*sam}
awk -v c="$c" '{print $0 "\t" c }' ${file%sam}length.txt > ${file%sam}length.name.txt
done

cat *trimmed.unique.mt.filter.length.name.txt > mt.length.dist.txt






