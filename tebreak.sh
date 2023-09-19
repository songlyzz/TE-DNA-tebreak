#!/bin/bash
#PBS -N tebreak-NTD-test
#PBS -l nodes=1:ppn=6
#PBS -l walltime=1000:00:00
#PBS -j oe
#PBS -q fat
data_dic="/histor/sun/huangjunsong/raw_data/TE_NTD_test"
work_dic="/histor/sun/huangjunsong/analysis_data/TE_pipeline_test"
ref_dic="/histor/sun/huangjunsong/reference"
software_dic="/histor/sun/huangjunsong/software"
ref_annotate="/histor/sun/huangjunsong/software/annovar/humandb/GRCh38/"
THREADS=16

mkdir ${work_dic}/LINE_results/results
mkdir ${work_dic}/LINE_results/annovar
mkdir ${data_dic}/map_data/bam_stat

source activate tebreak
for SAMPLE in A6163_A207 B1921_B60 B1921_B37 B1315_B113
do
#mapping
bwa mem -Y -M -t ${THREADS} ${ref_dic}/ensembl.Homo.GRCh38.fa ${data_dic}/raw_data/${SAMPLE}_1.clean.fq.gz ${data_dic}/raw_data/${SAMPLE}_2.clean.fq.gz > ${data_dic}/map_data/${SAMPLE}.sam
samtools sort -@ ${THREADS} -O BAM -o ${data_dic}/map_data/sorted_${SAMPLE}.bam ${data_dic}/map_data/${SAMPLE}.sam
echo "=========Mappppping ==================="

#质量评估
qualimap bamqc -bam sorted_${SAMPLE}.bam -outformat PDF:HTML -outdir ${data_dic}/map_data/bam_stat -nt ${THREADS} --java-mem-size=10G
echo "=========qualimap ==================="

#remove duplicates
#sambamba markdup --remove-duplicates -t ${THREADS} --tmpdir=${work_dic}/map_data/tmp/ ${work_dic}/map_data/sorted_${SAMPLE}.bam ${work_dic}/map_data/${SAMPLE}.rmdup.bam
#sambamba index ${work_dic}/map_data/${SAMPLE}.rmdup.bam
picard MarkDuplicates I=${data_dic}/map_data/sorted_${SAMPLE}.bam O=${data_dic}/map_data/${SAMPLE}.rmdup.bam M=${data_dic}/map_data/${SAMPLE}.rmdup.metrics.out
samtools index ${data_dic}/map_data/${SAMPLE}.rmdup.bam
echo "========Duplicate Removing ==========="

#obtain soft-clipped and discordant reads
cd ${data_dic}/map_data/
${software_dic}/tebreak/scripts/reduce_bam.py -b ${SAMPLE}.rmdup.bam #default = <input>.reduced.bam
echo "==============filter bam ==============="
#tebreak
tebreak --max_ins_reads 1000 -p ${THREADS} \
-b ${data_dic}/map_data/${SAMPLE}.rmdup.reduced.bam -r ${ref_dic}/ensembl.Homo.GRCh38.fa -o ${work_dic}/LINE_results/${SAMPLE}.tebreak.table.txt \
-d ${software_dic}/tebreak/lib/hg38.te.disctgt.txt -m ${software_dic}/tebreak/lib/hg38.centromere_telomere.bed -i ${software_dic}/tebreak/lib/teref.human.fa
echo "========tebreeeeeeeeeeaking ==========="

#filter
cd ${work_dic}/LINE_results/
${software_dic}/tebreak/scripts/general_filter.py -t ${work_dic}/LINE_results/${SAMPLE}.tebreak.table.txt -i ${software_dic}/tebreak/lib/teref.human.fa -r ${ref_dic}/ensembl.Homo.GRCh38.fa --numsplit 4 --numdiscord 4
${software_dic}/tebreak/scripts/annotate.py -t ${work_dic}/LINE_results/${SAMPLE}.tebreak.table.filter.txt -x ${software_dic}/tebreak/lib/nonref.collection.hg38.bed.gz -n KnownNonRef --nonref > ${work_dic}/LINE_results/results/${SAMPLE}.tebreak.txt
echo "=========ffffffilterrrrrr =============="

#annovar 
less -S ${work_dic}/LINE_results/results/${SAMPLE}.tebreak.txt | cut -f 2,3,4 | awk '{print $0"\\t""0""\\t""0"}' | grep -v Chromosome > ${work_dic}/LINE_results/annovar/${SAMPLE}.tebreak.annotate.input
echo "=========LESS ======================"

perl ${software_dic}/annovar/annotate_variation.pl -out ${work_dic}/LINE_results/annovar/${SAMPLE} -build GRCh38 ${work_dic}/LINE_results/annovar/${SAMPLE}.tebreak.annotate.input ${ref_annotate}
echo "=========ANNOVAR ======================"

cat ${work_dic}/LINE_results/annovar/${SAMPLE}.variant_function | cut -f 1 | awk -F "\\t" '{a[$1]++;sum++}END{for(i in a){print i"\\t"a[i]"\\t"a[i]/sum}}' | sort -k2 -nr | less > ${work_dic}/LINE_results/annovar/${SAMPLE}.proportion.txt
echo "=========CAT ======================"

done
rm ${data_dic}/map_data/*.sam
conda deactivate
