# Loading module and conda environment that has all the packages

module load anaconda3
source activate cuttag

# some common variable I am goint to use during the analysis
index="/data/Blizard-CGH-MadapuraLab/pankaj/genomes/Bowtie2Index/human/genome"
chrlen="/data/Blizard-CGH-MadapuraLab/pankaj/genomes/hg38chr"
ecoli="/data/Blizard-CGH-MadapuraLab/pankaj/genomes/Bowtie2Index/ecoli_K_12_MG1655/genome"
seacr="/data/Blizard-CGH-MadapuraLab/pankaj/seacr_script/SEACR_1.3.sh"

# repo about fastp: https://github.com/OpenGene/fastp
# or conda install -c bioconda fastp

# here my filename.txt file contains path of all the fastq file
# for loop iteratres over all the files and generates the trimmed fq reads
for file in $(cat filename.txt | grep '1.fq.gz')  # common mistake, doublecheck the filename.txt file's name
do
    input=${file%_[12].fq.gz}
    output=$(dirname $file | xargs basename )_trim
    fastp -w $NSLOTS -j ${output}.json -h ${output}.html -i ${input}_1.fq.gz -I ${input}_2.fq.gz -o ${output}_1.fq.gz -O ${output}_2.fq.gz
done


# In this for loop, I do bowtie2 alignment with hg38 genome 
# I also generate the bedgraph files required for the seacr peak calling 
 
for file in *trim_1.fq.gz
do
    output=${file%_[12].fq.gz}
    
    # For alignment, end to end, discard everything else except properly paired reads
    bowtie2 -t -p $NSLOTS --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x $index -1 ${output}_1.fq.gz -2 ${output}_2.fq.gz | samtools view -bSh - -o ${output}_unsorted.bam
    samtools sort -o ${output}_sorted.bam ${output}_unsorted.bam
    rm ${output}_unsorted.bam
    
    # as the Cut and tag /cut and run enzymes are purified from the e.coli, the author the those papers argue the residual gDNA in the enzyme can be used
    # for normalization but I find these not very reliable
    bowtie2 -t -p $NSLOTS --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x $ecoli -1 ${output}_1.fq.gz -2 ${output}_2.fq.gz | samtools view -bSh - -o ${output}_ecoli_unsorted.bam
    samtools sort -o ${output}_ecoli_sorted.bam ${output}_ecoli_unsorted.bam
    rm ${output}_ecoli_unsorted.bam
    
    # index file
    samtools index ${output}_sorted.bam ${output}_sorted.bam.bai

    # for bigwig file generation
    bamCoverage -b ${output}_sorted.bam -o ${output}_cpm.bigwig --binSize 10 --normalizeUsing CPM

    # this is to generate the final bedgraph file that can would be used for seacr peak calling
    samtools view -uh -f 0x2 ${output}_sorted.bam | samtools sort -n -m 5G -@ 5 - | samtools fixmate - - | bedtools bamtobed -i stdin -bedpe | bedtools sort -i stdin | awk '$6-$2 < 1000 {print $1"\t"$2"\t"$6 }' | bedtools genomecov -bg -i stdin -g $chrlen > ${output}.bedgraph
done

# MACS2
# Unlike, Sonication driven fragmentation of DNA, enzymatic fragmentation are expected to have exact start and end position,
# for --keep-dup I consider keeping it 'all' or 10

macs2 callpeak --broad -t $file -c $control -f BAMPE -g hs --keep-dup 10 --name $name --outdir $outdir