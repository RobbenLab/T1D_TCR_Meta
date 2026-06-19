#!/bin/bash

#The purpose of this script is to take in the RNAseq data and output the data on them. The things we need:
#  RNA libary statistics (n reads, filtered, etc....)
#  HLA typing (outputs a genotype.tsv file, we probably want to get that into a )
#  TCR gene usage (MixCR we will get fully assembled which probably won't include all aligned genes but unique full lengths)

INPUT=$1 #This needs to be a list of the sample names and files for processing
FDIR=$2 #This is where the data is located
OUTDIR=$3 #The output directory (should already exist)
TMPDIR=$4 #Temp directory to hold files that will be deleted by the end
NUM_CPU=60

#mkdir $OUTDIR
#if [[ -d "$OUTDIR" ]]; then
#    : #Do nothing
#else
#    echo "Directory $OUTDIR does not exist."
#    exit 1
#fi

mkdir $TMPDIR
if [[ -d "$TMPDIR" ]]; then
  : #Do nothing
else
    echo "Directory $TMPDIR does not exist."
    exit 1
fi

#mkdir "$OUTDIR/HLA"
mkdir "$OUTDIR/MixCR/"
mkdir "$OUTDIR/Stat/"
mkdir "$OUTDIR/Counts/"

Paired="false"
while getopts "p" flag; do
  case "${flag}" in
    p) 
      Paired="true";; 
  esac
done
f="$INPUTDIR/*"


#Save each column of the input file to a variable so that we can for loop through
mapfile -t names < <(awk 'BEGIN{FS="\t"} {print $1}' "$INPUT")
mapfile -t r1 < <(awk 'BEGIN{FS="\t"} {print $2}' "$INPUT")
mapfile -t r2 < <(awk 'BEGIN{FS="\t"} {print $3}' "$INPUT")

input_length=${#names[@]}

if $Paired; then
  echo "Paired end dataset"
fi

#Create start genome index to start (takes about 
#STAR --runThreadN 120 \
#  --runMode genomeGenerate \
#  --genomeDir /data/Databases/Genomes/GRCh38_p14/STAR/ \
#  --genomeFastaFiles /data/Databases/Genomes/GRCh38_p14/GCF_000001405.40_GRCh38.p14_genomic.fna \
#  --sjdbGTFfile /data/Databases/Genomes/GRCh38_p14/genomic.gff \
#  --sjdbOverhang 99




###########################################################
#####       Start Loop ####################################
###########################################################


echo "Performing operations on $input_length files"
times=()

for i in "${!names[@]}"; do
  
  start_time=$SECONDS
  #Extract variables
  name=${names[$i]}
  file1=$FDIR/${r1[$i]}
  file2=$FDIR/${r2[$i]}
  
  echo "###############################################"
  echo "######### Index: $i, Sample: $name ###########"
  echo "###############################################"
  
#  echo $file1
#  echo $file2
  
  #####################################################
  ###########  Statistics  ############################
  #####################################################
  
  
  echo "Performing Trimming...."
  #Perform trimming on the file and save the trimmed to the   
  if $Paired; then
    cutadapt --report=minimal -q 15 --cores $NUM_CPU -m 10 --max-average-error-rate 0.3 -o $TMPDIR/Trimmed_${name}_r1.fastq.gz -p $TMPDIR/Trimmed_${name}_r2.fastq.gz $file1 $file2 1> "$OUTDIR/Stat/$name.txt"
    trim_file1="$TMPDIR/Trimmed_${name}_r1.fastq.gz"
    trim_file2="$TMPDIR/Trimmed_${name}_r2.fastq.gz"
  else
    cutadapt --report=minimal -q 15 --cores $NUM_CPU -m 10 --max-average-error-rate 0.3 -o $TMPDIR/Trimmed_$name.fastq.gz $file1 1> "$OUTDIR/Stat/$name.txt"
    trim_file1="$TMPDIR/Trimmed_${name}.fastq.gz"
  fi
  echo "Finished Trimming...."
  
  ######################################################
  ##########   Star Alignment  #########################
  ######################################################
  
  echo "Performing Star alignment...."
  
  if $Paired; then
    STAR --genomeDir /data/Databases/Genomes/GRCh38_p14/STAR/ \
      --runThreadN $NUM_CPU \
      --readFilesIn $trim_file1 $trim_file2 \
      --outFileNamePrefix $TMPDIR/Counts/ \
      --sjdbGTFfile /data/Databases/Genomes/GRCh38_p14/genomic.gtf \
      --outSAMtype BAM Unsorted \
      --readFilesCommand zcat
      
  else
    STAR --genomeDir /data/Databases/Genomes/GRCh38_p14/STAR/ \
      --runThreadN $NUM_CPU \
      --readFilesIn $trim_file1 \
      --outFileNamePrefix $TMPDIR/Counts/ \
      --sjdbGTFfile /data/Databases/Genomes/GRCh38_p14/genomic.gtf \
      --outSAMtype BAM Unsorted \
      --readFilesCommand zcat
      
  
  fi
  
  echo "Starting transcript counts..."
  
  #Run the counts
  samtools sort -n $TMPDIR/Counts/Aligned.out.bam -o $TMPDIR/Counts/sorted_by_name.bam
  htseq-count -q -t transcript -i gene_id -f bam -r name -s no $TMPDIR/Counts/sorted_by_name.bam /data/Databases/Genomes/GRCh38_p14/genomic.gtf > $OUTDIR/Counts/${name}_counts.tsv 
  
  echo "Finished STAR align and count....."
  
  
  
  ######################################################
  ############  HLA Typing #############################
  ######################################################

  mkdir -p "$OUTDIR/HLA/$name/"
  echo "Performing HLA Typing...."
  if $Paired; then
    #This command produces ~ 15 MB for 215 MB 
    run-t1k -1 $trim_file1 -2 $trim_file2 --preset hla -t $NUM_CPU -f  ~/Projects/T1D_RNA/HLA/hlaidx_rna_seq.fa --od "$OUTDIR/HLA/$name/"
  else
    run-t1k -u $trim_file1 --preset hla -t $NUM_CPU -f ~/Projects/T1D_RNA/HLA/hlaidx_rna_seq.fa --od "$OUTDIR/HLA/$name/"
  fi
  echo "Finishing HLA Typing...."

  #####################################################
  ##########   MixCR ##################################
  #####################################################

  mkdir -p "$OUTDIR/MixCR/$name/"
  
  echo "Performing MixCR...."
  #This will get full but we can expect only like 180 full sequences from 4 mil reads we will compare to alignment of TCR genes
  if $Paired; then
    mixcr analyze rna-seq -t $NUM_CPU -f \
      --species hsa \
      $trim_file1 \
      $trim_file2 \
      "$OUTDIR/MixCR/$name/" >/dev/null 2>&1
  else
    mixcr analyze rna-seq -t $NUM_CPU -f \
        --species hsa \
        $trim_file1 \
        "$OUTDIR/MixCR/$name/" >/dev/null 2>&1
  fi
  echo "Finishing MixCR"
  
  #End stuff
  end_time=$SECONDS 
  duration=$((end_time - start_time))
  times+=$duration
  
  
  
done


printf "%s\n" "${times[@]}" > "$OUTDIR/times.txt"
#Now delete temporary files
rm -rf $TMPDIR
