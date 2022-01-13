#!/bin/bash
echo "Started RNA Seq script"


rnaSeq=/projects/micb405/analysis/Group1/rnaSeq
refGenome=/projects/micb405/analysis/Group1/rnaSeq/mm9Genome.fasta
textFiles=/projects/micb405/analysis/Group1/textFiles/SraSpecificRnaSeq.txt
gtfFile=/projects/micb405/analysis/Group1/rnaSeq/referenceGTFmm9/gencode.vM1.annotation.gtf
starFiles=/projects/micb405/analysis/Group1/rnaSeq/starSRRFiles

#downloading each SRR file from the server and output into the dir chipSeq

while read sraNum; do
      echo "Downloading \"$sraNum\"..."
      fasterq-dump --split-3 "$sraNum" -O "$rnaSeq"
done < $textFiles


# Index refGenome
STAR \
        --runMode genomeGenerate \
        --genomeDir /projects/micb405/analysis/Group1/rnaSeq/STARIndex \
        --genomeFastaFiles "$refGenome" \
        --sjdbGTFfile "$gtfFile" \
        --sjdbOverhang 100 \
        --runThreadN 16



# ALign files
for file in "$rnaSeq"/*.gz
do
        echo /projects/micb405/analysis/Group1/rnaSeq/SRR1536410 | cut -d / -f7
        path=`echo "$file" | cut -d . -f1`
        basename=`echo "$path" | cut -d / -f7`
        STAR \
                --genomeDir /projects/micb405/analysis/Group1/rnaSeq/STARIndex/ \
                --readFilesIn "$file" \
                --outFileNamePrefix "$path" \
                --runThreadN 8 \
                --limitBAMsortRAM 60000000000 \
                --outSAMattrRGline ID:"$basename".fastq SM:"$basename".fastq \
                --outBAMsortingThreadN 8 \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --outSAMstrandField intronMotif \
                --readFilesCommand zcat \
                --chimSegmentMin 20 \
                --genomeLoad NoSharedMemory

        mv *.out* /projects/micb405/analysis/Group1/rnaSeq/starSRRFiles  
        mv *_STARtmp /projects/micb405/analysis/Group1/rnaSeq/starSRRFiles 
        
        echo file
        echo "$file"
        echo "$basename"
done


# Performing htseq-count
for file in "$starFiles"/*.sortedByCoord.out.bam
do
        path=`echo "$file" | cut -d . -f1`
        basename=`echo "$path" | cut -d / -f8`
        samtools index "$file"

        htseq-count \
                -f bam \
                -r pos \
                --stranded=no \
                "$starFiles"/"$basename".sortedByCoord.out.bam \
                "$gtfFile" \
                > "$basename".htseq.out
done
