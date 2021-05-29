work="/public/home/fyy/rna_seq/hcmv/110"
db="/public/home/fyy/db"

for i in 516 517 518 519 520 521 522 523 524 525 526 527 528 529 530 
do
	mkdir SRR6326${i}
	
	source activate python27

	kneaddata  -i $work/SRR6326${i}_1.fastq.gz \
	-i $work/SRR6326${i}_2.fastq.gz \
	-o $work/SRR6326${i} \
	-v -t 15 --remove-intermediate-output \
	--trimmomatic /public/home/fyy/miniconda3/envs/python27/share/trimmomatic/ \
	--trimmomatic-options 'ILLUMINACLIP:TruSeq3-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50' \
	--bowtie2-options '--very-sensitive --dovetail' \
	--reference-db $db/human/Homo_sapiens \
	-db $db/rRNA/SILVA_128_LSUParc_SSUParc_ribosomal_RNA

	conda deactivate
	source activate rnaseq
	
	hisat2 -x /public/home/fyy/biosoft/virus/virus -p 10 \
	-1 $work/SRR6326${i}/SRR6326${i}_1_kneaddata_paired_1.fastq \
	-2 $work/SRR6326${i}/SRR6326${i}_1_kneaddata_paired_2.fastq \
	-S $work/SRR6326${i}/SRR6326${i}.sam

	samtools view -S $work/SRR6326${i}/SRR6326${i}.sam -b > $work/SRR6326${i}/SRR6326${i}.bam
	samtools sort $work/SRR6326${i}/SRR6326${i}.bam -o $work/SRR6326${i}/SRR6326${i}_sorted.bam
	samtools index $work/SRR6326${i}/SRR6326${i}_sorted.bam

	featureCounts -t exon -g gene_name \
	-a /public/home/fyy/biosoft/virus/gtf/virus.gtf \
	-o $work/SRR6326${i}/SRR6326${i}_counts_virus.txt \
	$work/SRR6326${i}/SRR6326${i}_sorted.bam
	
	conda deactivate
	cut -f 1,7 $work/SRR6326${i}/SRR6326${i}_counts_virus.txt |grep -v '^#' >$work/SRR6326${i}/SRR6326${i}_allvirus.txt
	grep 'NC_006273.2' $work/SRR6326${i}/SRR6326${i}.sam > $work/SRR6326${i}/SRR6326${i}_6237.txt
	sed -i '1d' $work/SRR6326${i}/SRR6326${i}_6237.txt
	cut -f 3,10 $work/SRR6326${i}/SRR6326${i}_6237.txt >$work/SRR6326${i}/SRR6326${i}_hcmv.txt
	awk '{ printf ">%s\n%s\n",$1,$2 }' $work/SRR6326${i}/SRR6326${i}_hcmv.txt > $work/SRR6326${i}/SRR6326${i}_hcmv.fasta
	
	if [ -s $work/SRR6326${i}/SRR6326${i}_6237.txt ]
	then
	source activate ciriquant
	bwa mem /public/home/fyy/db/hcmv/bwa_hcmv/bwahcmv $work/SRR6326${i}/SRR6326${i}_hcmv.fasta > $work/SRR6326${i}/SRR6326${i}_hcmv.sam
	conda deactivate
        source activate rnaseq

        samtools view -S $work/SRR6326${i}/SRR6326${i}_hcmv.sam -b > $work/SRR6326${i}/SRR6326${i}_hcmv.bam
        samtools sort $work/SRR6326${i}/SRR6326${i}_hcmv.bam -o $work/SRR6326${i}/SRR6326${i}_hcmv_sorted.bam
        samtools index $work/SRR6326${i}/SRR6326${i}_hcmv_sorted.bam

        featureCounts -t exon -g gene_name \
        -a /public/home/fyy/db/hcmv/bwa_hcmv/hcmv.gtf \
        -o $work/SRR6326${i}/SRR6326${i}_hcmv_virus_last.txt \
        $work/SRR6326${i}/SRR6326${i}_hcmv_sorted.bam

        conda deactivate

	cut -f 1,7 $work/SRR6326${i}/SRR6326${i}_hcmv_virus_last.txt |grep -v '^#' >$work/SRR6326${i}/SRR6326${i}_last.txt	

	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_Homo_sapiens_bowtie2_paired_contam_1.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_Homo_sapiens_bowtie2_paired_contam_2.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_Homo_sapiens_bowtie2_unmatched_1_contam.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_Homo_sapiens_bowtie2_unmatched_2_contam.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_unmatched_1_contam.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_unmatched_2_contam.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_1.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_2.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_unmatched_1.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_unmatched_2.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_paired_1.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_paired_2.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}.bam
	rm -f $work/SRR6326${i}/SRR6326${i}_sorted.bam
	rm -f $work/SRR6326${i}/SRR6326${i}_sorted.bam
	rm -f $work/SRR6326${i}_1.fastq.gz
	rm -f $work/SRR6326${i}_2.fastq.gz
	

	else
	
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_Homo_sapiens_bowtie2_paired_contam_1.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_Homo_sapiens_bowtie2_paired_contam_2.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_Homo_sapiens_bowtie2_unmatched_1_contam.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_Homo_sapiens_bowtie2_unmatched_2_contam.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_unmatched_1_contam.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_unmatched_2_contam.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_1.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_2.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_unmatched_1.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_unmatched_2.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_paired_1.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}_1_kneaddata_paired_2.fastq
	rm -f $work/SRR6326${i}/SRR6326${i}.bam
	rm -f $work/SRR6326${i}/SRR6326${i}_sorted.bam
	rm -f $work/SRR6326${i}/SRR6326${i}_sorted.bam
	rm -f $work/SRR6326${i}_1.fastq.gz
	rm -f $work/SRR6326${i}_2.fastq.gz
	fi
done
