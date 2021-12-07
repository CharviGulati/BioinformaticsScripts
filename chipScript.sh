#!/bin/bash
echo "Started script"



textFiles=/projects/micb405/analysis/Group1/textFiles/SraSpecificChipSeq.txt
chipSeq=/projects/micb405/analysis/Group1/chipSeq
RefGenome=/projects/micb405/analysis/Group1/referenceGenome/mm9Genome.fasta
gtfFile=/projects/micb405/analysis/Group1/chipSeq/referenceGTFmm9/gencode.vM1.annotation.gtf.bed

macsFiles=/projects/micb405/analysis/Group1/chipSeq/files/macsFiles/




#referenceGenome is indexed already

#downloading each SRR file from the server and output into the dir chipSeq

while read sraNum; do
      echo "Downloading \"$sraNum\"..."
      fasterq-dump --split-3 "$sraNum" -O "$chipSeq"
done < $textFiles



# a simple loop that converts each SRR file from a FASTQ -> BAM while filtering our dups
# then removes unneeded files

for file in "$chipSeq"/*.fastq
do
#    echo "$file"
#    echo $file
#    SRR1521828.fastq
    echo "Aligning SRR to Reference Genome..."
    bwa aln "$RefGenome" "$file" > aligned.sai
    bwa samse "$RefGenome" aligned.sai "$file" > file.sam
    echo "Converting SAM to BAM..."
    samtools view -h -b file.sam > file.bam
    echo "Taging mates..."
    samtools fixmate -m file.bam file.fixmate.bam
    echo "Position sort for markdup..."
    samtools sort file.fixmate.bam -o fixmateSorted.bam
    echo "Marking duplicates..."
    samtools markdup -S fixmateSorted.bam markdup.bam
    echo "Sorting again..."
    samtools sort markdup.bam -o fileSortedFinal.bam
    echo "indexing BAM file..."
    samtools index fileSortedFinal.bam

    basename=`echo "$file" | cut -d . -f1`
    mv fileSortedFinal.bam $basename.sorted.bam 

    echo "Removing temporary files..."
    rm file.sam
    rm file.bam
    rm fileSorted.bam
    rm file.fixmate.bam
    rm fixmateSorted.bam
    rm markdup.bam
done


# using macs2 to call peaks of the finals sorted bam file for each SRR keeping in mind that the false discovery rate is 0.001
for file in "$chipSeq"/*.bam;
do
        echo "Using MACS2 for peak calling"
        basename=`echo "$file" | cut -d . -f1`
        macs2 callpeak -t "$file" -g mm -n $basename.bdg -B -q 0.001
done



for file in "$chipSeq"/*treat_pileup.bdg;
do
        echo "File is being sorted and compressed into BigWig"
        bedtools sort -i "$file" > ${file%.*}.sorted.bdg
        bedGraphToBigWig "$file" mm9.chrom.sizes "$file".sorted.bw
done







# echo "STARTING SCRIPT FOR PEAK INFO" > infoLog2.txt

for file in "$macsFiles"/*.bdg_peaks.narrowPeak;
do
        path=`echo "$file" | cut -d . -f1`
        basename=`echo "$path" | cut -d / -f10`

        echo "$path"            #/projects/micb405/analysis/Group1/chipSeq/files/macsFiles//SRR1521828
        echo "$basename"        #SRR1521834

        echo "CHiP Seq data on "$basename" " >> infoLog2.txt
        
        cp "$basename".bdg_peaks.narrowPeak "$basename".narrowPeak.bed
        #presort data by chromosome and then by start position
        sort -k1,1 -k2,2n "$basename".narrowPeak.bed > "$basename".sorted.bed
        # What genes (if any) fall directly under the peaks, and how large is the overlap?
        # For very large B files, invoke a “sweeping” algorithm that requires position-sorted 
        # input. When using -sorted, memory usage remains low even for very large files.
       
        bedtools intersect -sorted -a "$basename".sorted.bed -b "$gtfFile" > "$basename".narrowPeak.genes.bed
        # the overlap size is the last column

        
        echo "How many total overlaps are there" >> infoLog2.txt
        wc -l "$basename".narrowPeak.genes.bed >> infoLog2.txt
        
        echo "How many distinct genes are overlapped?" >> infoLog2.txt
        cut -f 10 "$basename".narrowPeak.genes.bed | sort -u | wc -l >> infoLog2.txt

        echo "---------------------" >> infoLog2.txt

done







echo "STARTING SCRIPT FOR PEAK INFO" > infoLog.txt
echo "CHiP Seq data on SRR1521829 and SRR1521831 for ST-HSC" >> infoLog.txt
cp SRR1521829.bdg_peaks.narrowPeak SRR1521829.narrowPeak.bed
cp SRR1521831.bdg_peaks.narrowPeak SRR1521831.narrowPeak.bed
        
#presort data by chromosome and then by start position
sort -k1,1 -k2,2n SRR1521829.narrowPeak.bed > SRR1521829.sorted.bed
sort -k1,1 -k2,2n SRR1521831.narrowPeak.bed > SRR1521831.sorted.bed
# What genes (if any) fall directly under the peaks, and how large is the overlap?
# For very large B files, invoke a “sweeping” algorithm that requires position-sorted 
# input. When using -sorted, memory usage remains low even for very large files.
       
bedtools intersect -sorted -a SRR1521829.sorted.bed -b SRR1521831.sorted.bed > STHSC.narrowPeak.bed
# the overlap size is the last column

echo "How many total overlaps are there" >> infoLog.txt
wc -l STHSC.narrowPeak.bed >> infoLog.txt
        
echo "How many distinct genes are overlapped?" >> infoLog.txt
cut -f 10 STHSC.narrowPeak.bed | sort -u | wc -l >> infoLog.txt

echo "------------------------------------------------------------------------------" >> infoLog.txt



echo "CHiP Seq data on SRR1521832 and SRR1521833 for MPP" >> infoLog.txt
cp SRR1521832.bdg_peaks.narrowPeak SRR1521832.narrowPeak.bed
cp SRR1521833.bdg_peaks.narrowPeak SRR1521833.narrowPeak.bed
        
#presort data by chromosome and then by start position
sort -k1,1 -k2,2n SRR1521832.narrowPeak.bed > SRR1521832.sorted.bed
sort -k1,1 -k2,2n SRR1521833.narrowPeak.bed > SRR1521833.sorted.bed
# What genes (if any) fall directly under the peaks, and how large is the overlap?
# For very large B files, invoke a “sweeping” algorithm that requires position-sorted 
# input. When using -sorted, memory usage remains low even for very large files.
       
bedtools intersect -sorted -a SRR1521832.sorted.bed -b SRR1521833.sorted.bed > MPP.narrowPeak.bed
# the overlap size is the last column

echo "How many total overlaps are there" >> infoLog.txt
wc -l MPP.narrowPeak.bed >> infoLog.txt
        
echo "How many distinct genes are overlapped?" >> infoLog.txt
cut -f 10 MPP.narrowPeak.bed | sort -u | wc -l >> infoLog.txt

echo "------------------------------------------------------------------------------" >> infoLog.txt




echo "CHiP Seq data on SRR1521834 and SRR1521836 for CMP" >> infoLog.txt
cp SRR1521834.bdg_peaks.narrowPeak SRR1521834.narrowPeak.bed
cp SRR1521836.bdg_peaks.narrowPeak SRR1521836.narrowPeak.bed
        
#presort data by chromosome and then by start position
sort -k1,1 -k2,2n SRR1521834.narrowPeak.bed > SRR1521834.sorted.bed
sort -k1,1 -k2,2n SRR1521836.narrowPeak.bed > SRR1521836.sorted.bed
# What genes (if any) fall directly under the peaks, and how large is the overlap?
# For very large B files, invoke a “sweeping” algorithm that requires position-sorted 
# input. When using -sorted, memory usage remains low even for very large files.
       
bedtools intersect -sorted -a SRR1521834.sorted.bed -b SRR1521836.sorted.bed > CMP.narrowPeak.bed
# the overlap size is the last column

echo "How many total overlaps are there" >> infoLog.txt
wc -l CMP.narrowPeak.bed >> infoLog.txt
        
echo "How many distinct genes are overlapped?" >> infoLog.txt
cut -f 10 CMP.narrowPeak.bed | sort -u | wc -l >> infoLog.txt

echo "------------------------------------------------------------------------------" >> infoLog.txt




echo "CHiP Seq data on SRR1521837 and SRR1521838 for GMP " >> infoLog.txt
cp SRR1521837.bdg_peaks.narrowPeak SRR1521837.narrowPeak.bed
cp SRR1521838.bdg_peaks.narrowPeak SRR1521838.narrowPeak.bed
        
#presort data by chromosome and then by start position
sort -k1,1 -k2,2n SRR1521837.narrowPeak.bed > SRR1521837.sorted.bed
sort -k1,1 -k2,2n SRR1521838.narrowPeak.bed > SRR1521838.sorted.bed
# What genes (if any) fall directly under the peaks, and how large is the overlap?
# For very large B files, invoke a “sweeping” algorithm that requires position-sorted 
# input. When using -sorted, memory usage remains low even for very large files.
       
bedtools intersect -sorted -a SRR1521837.sorted.bed -b SRR1521838.sorted.bed > GMP.narrowPeak.bed
# the overlap size is the last column

echo "How many total overlaps are there" >> infoLog.txt
wc -l GMP.narrowPeak.bed >> infoLog.txt
        
echo "How many distinct genes are overlapped?" >> infoLog.txt
cut -f 10 GMP.narrowPeak.bed | sort -u | wc -l >> infoLog.txt

echo "------------------------------------------------------------------------------" >> infoLog.txt




echo "CHiP Seq data on SRR1521851 and SRR1521854 for BMM" >> infoLog.txt
cp SRR1521851.bdg_peaks.narrowPeak SRR1521851.narrowPeak.bed
cp SRR1521854.bdg_peaks.narrowPeak SRR1521854.narrowPeak.bed
        
#presort data by chromosome and then by start position
sort -k1,1 -k2,2n SRR1521851.narrowPeak.bed > SRR1521851.sorted.bed
sort -k1,1 -k2,2n SRR1521854.narrowPeak.bed > SRR1521854.sorted.bed
# What genes (if any) fall directly under the peaks, and how large is the overlap?
# For very large B files, invoke a “sweeping” algorithm that requires position-sorted 
# input. When using -sorted, memory usage remains low even for very large files.
       
bedtools intersect -sorted -a SRR1521851.sorted.bed -b SRR1521854.sorted.bed > BMM.narrowPeak.bed
# the overlap size is the last column

echo "How many total overlaps are there" >> infoLog.txt
wc -l BMM.narrowPeak.bed >> infoLog.txt
        
echo "How many distinct genes are overlapped?" >> infoLog.txt
cut -f 10 BMM.narrowPeak.bed | sort -u | wc -l >> infoLog.txt

echo "------------------------------------------------------------------------------" >> infoLog.txt





echo "CHiP Seq data on SRR1521845 and SRR1521846 " >> infoLog.txt
cp SRR1521845.bdg_peaks.narrowPeak SRR1521845.narrowPeak.bed
cp SRR1521846.bdg_peaks.narrowPeak SRR1521846.narrowPeak.bed
        
#presort data by chromosome and then by start position
sort -k1,1 -k2,2n SRR1521845.narrowPeak.bed > SRR1521845.sorted.bed
sort -k1,1 -k2,2n SRR1521846.narrowPeak.bed > SRR1521846.sorted.bed
# What genes (if any) fall directly under the peaks, and how large is the overlap?
# For very large B files, invoke a “sweeping” algorithm that requires position-sorted 
# input. When using -sorted, memory usage remains low even for very large files.
       
bedtools intersect -sorted -a SRR1521845.sorted.bed -b SRR1521846.sorted.bed > G.narrowPeak.bed
# the overlap size is the last column

echo "How many total overlaps are there" >> infoLog.txt
wc -l G.narrowPeak.bed >> infoLog.txt
        
echo "How many distinct genes are overlapped?" >> infoLog.txt
cut -f 10 G.narrowPeak.bed | sort -u | wc -l >> infoLog.txt

echo "------------------------------------------------------------------------------" >> infoLog.txt






echo "CHiP Seq data on SRR1521848 and SRR1521850 " >> infoLog.txt
cp SRR1521848.bdg_peaks.narrowPeak SRR1521848.narrowPeak.bed
cp SRR1521850.bdg_peaks.narrowPeak SRR1521850.narrowPeak.bed
        
#presort data by chromosome and then by start position
sort -k1,1 -k2,2n SRR1521848.narrowPeak.bed > SRR1521848.sorted.bed
sort -k1,1 -k2,2n SRR1521850.narrowPeak.bed > SRR1521850.sorted.bed
# What genes (if any) fall directly under the peaks, and how large is the overlap?
# For very large B files, invoke a “sweeping” algorithm that requires position-sorted 
# input. When using -sorted, memory usage remains low even for very large files.
       
bedtools intersect -sorted -a SRR1521848.sorted.bed -b SRR1521850.sorted.bed > M.narrowPeak.bed
# the overlap size is the last column

echo "How many total overlaps are there" >> infoLog.txt
wc -l M.narrowPeak.bed >> infoLog.txt
        
echo "How many distinct genes are overlapped?" >> infoLog.txt
cut -f 10 M.narrowPeak.bed | sort -u | wc -l >> infoLog.txt

echo "------------------------------------------------------------------------------" >> infoLog.txt



echo "bedtools jaccard LTHSC and LTHSC"
bedtools jaccard -a LTHSC.narrowPeak.bed -b LTHSC.narrowPeak.bed
echo "bedtools jaccard LTHSC and STHSC"
bedtools jaccard -a LTHSC.narrowPeak.bed -b STHSC.narrowPeak.bed
echo "bedtools jaccard LTHSC and MPP"
bedtools jaccard -a LTHSC.narrowPeak.bed -b MPP.narrowPeak.bed
echo "bedtools jaccard LTHSC and CMP"
bedtools jaccard -a LTHSC.narrowPeak.bed -b CMP.narrowPeak.bed
echo "bedtools jaccard LTHSC and GMP"
bedtools jaccard -a LTHSC.narrowPeak.bed -b GMP.narrowPeak.bed
echo "bedtools jaccard LTHSC and BMM"
bedtools jaccard -a LTHSC.narrowPeak.bed -b BMM.narrowPeak.bed
echo "bedtools jaccard LTHSC and G"
bedtools jaccard -a LTHSC.narrowPeak.bed -b G.narrowPeak.bed
echo "bedtools jaccard LTHSC and M"
bedtools jaccard -a LTHSC.narrowPeak.bed -b M.narrowPeak.bed


echo "bedtools jaccard STHSC and STHSC"
bedtools jaccard -a STHSC.narrowPeak.bed -b STHSC.narrowPeak.bed
echo "bedtools jaccard STHSC and MPP"
bedtools jaccard -a STHSC.narrowPeak.bed -b MPP.narrowPeak.bed
echo "bedtools jaccard STHSC and CMP"
bedtools jaccard -a STHSC.narrowPeak.bed -b CMP.narrowPeak.bed
echo "bedtools jaccard STHSC and GMP"
bedtools jaccard -a STHSC.narrowPeak.bed -b GMP.narrowPeak.bed
echo "bedtools jaccard STHSC and BMM"
bedtools jaccard -a STHSC.narrowPeak.bed -b BMM.narrowPeak.bed
echo "bedtools jaccard STHSC and G"
bedtools jaccard -a STHSC.narrowPeak.bed -b G.narrowPeak.bed
echo "bedtools jaccard STHSC and M"
bedtools jaccard -a STHSC.narrowPeak.bed -b M.narrowPeak.bed

echo "bedtools jaccard MPP and MPP"
bedtools jaccard -a MPP.narrowPeak.bed -b MPP.narrowPeak.bed
echo "bedtools jaccard MPP and CMP"
bedtools jaccard -a MPP.narrowPeak.bed -b CMP.narrowPeak.bed
echo "bedtools jaccard MPP and GMP"
bedtools jaccard -a MPP.narrowPeak.bed -b GMP.narrowPeak.bed
echo "bedtools jaccard MPP and BMM"
bedtools jaccard -a MPP.narrowPeak.bed -b BMM.narrowPeak.bed
echo "bedtools jaccard MPP and G"
bedtools jaccard -a MPP.narrowPeak.bed -b G.narrowPeak.bed
echo "bedtools jaccard MPP and M"
bedtools jaccard -a MPP.narrowPeak.bed -b M.narrowPeak.bed

echo "bedtools jaccard CMP and CMP"
bedtools jaccard -a CMP.narrowPeak.bed -b CMP.narrowPeak.bed
echo "bedtools jaccard CMP and GMP"
bedtools jaccard -a CMP.narrowPeak.bed -b GMP.narrowPeak.bed
echo "bedtools jaccard CMP and BMM"
bedtools jaccard -a CMP.narrowPeak.bed -b BMM.narrowPeak.bed
echo "bedtools jaccard CMP and G"
bedtools jaccard -a CMP.narrowPeak.bed -b G.narrowPeak.bed
echo "bedtools jaccard CMP and M"
bedtools jaccard -a CMP.narrowPeak.bed -b M.narrowPeak.bed

echo "bedtools jaccard GMP and GMP"
bedtools jaccard -a GMP.narrowPeak.bed -b GMP.narrowPeak.bed
echo "bedtools jaccard GMP and BMM"
bedtools jaccard -a GMP.narrowPeak.bed -b BMM.narrowPeak.bed
echo "bedtools jaccard GMP and G"
bedtools jaccard -a GMP.narrowPeak.bed -b G.narrowPeak.bed
echo "bedtools jaccard GMP and M"
bedtools jaccard -a GMP.narrowPeak.bed -b M.narrowPeak.bed

echo "bedtools jaccard BMM and BMM"
bedtools jaccard -a BMM.narrowPeak.bed -b BMM.narrowPeak.bed
echo "bedtools jaccard BMM and G"
bedtools jaccard -a BMM.narrowPeak.bed -b G.narrowPeak.bed
echo "bedtools jaccard BMM and M"
bedtools jaccard -a BMM.narrowPeak.bed -b M.narrowPeak.bed

echo "bedtools jaccard G and G"
bedtools jaccard -a G.narrowPeak.bed -b G.narrowPeak.bed
echo "bedtools jaccard G and M"
bedtools jaccard -a G.narrowPeak.bed -b M.narrowPeak.bed

echo "bedtools jaccard M and M"
bedtools jaccard -a M.narrowPeak.bed -b M.narrowPeak.bed