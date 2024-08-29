#Basecaller: 6.4.2
#Flowcell: R10.4.1
# https://github.com/calacademy-research/minibar

#!/bin/bash

conda activate longread_umi_custom

#wget https://raw.githubusercontent.com/SorenKarst/longread_umi/develop/scripts/install_conda.sh

#bash ./install_conda.sh develop

#/data/thesisguest1/Tools/usearch

############# SETTING VARIABLES #############

#These should be fixed as standard
#set working directory
wrk=/data/thesisguest1/work/UMI

#tools
minibar=/data/thesisguest1/Tools/minibar.py

#set the following parameters
#what is the name of your project folder where everything will run?

run_ID=RUN19

prj=$wrk/data_"$run_ID"

sample_id="$run_ID"_guppy_v6.4.2_sup_BC

raw_data_loc=$prj/raw_fastq_guppy_v6.4.2_sup

cd $prj


for bc in {22..22}; do 
  mkdir $prj/$sample_id"$bc"
  mkdir $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/
done


#only do barcode05

bc=22
cd $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/
nohup $minibar $prj/Replicate_ID.txt $prj/$sample_id"$bc".fastq -F -l 600 &


############# FILTER FILES for HIV #############
#https://www.samformat.info/sam-format-flag
conda activate longread_umi_custom

for rep_id in {1..8};do
minimap2 -t 15 -ax map-ont /data/thesisguest1/Pipe_3SEP20/db/blast/HXB2.fasta $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/sample_ID"$rep_id".fastq | samtools fastq -n -f 0x4 - > $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/sample_ID"$rep_id".nonHIV.fastq ;
minimap2 -t 15 -ax map-ont /data/thesisguest1/Pipe_3SEP20/db/blast/HXB2.fasta $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/sample_ID"$rep_id".fastq  | samtools fastq -n -F 0x4 - > $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/sample_ID"$rep_id".HIV.fastq ;
done

############# Check reads in files #############
conda activate longread_umi

for rep_id in {1..8}; do echo "$rep_id";
NanoComp --fastq $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/sample_ID"$rep_id".HIV.fastq $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/sample_ID"$rep_id".nonHIV.fastq -n ID"$rep_id"_HIV ID"$rep_id"_nonHIV --outdir compare-ID"$rep_id"_HIV
done

conda activate longread_umi

for rep_id in {1..8}; do echo "$rep_id";
NanoPlot -t 20 --fastq $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/sample_ID"$rep_id".fastq --maxlength 10000 -o $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/summary_ID"$rep_id"
done


NanoComp --fastq sample_ID1.fastq sample_ID2.fastq sample_ID3.fastq sample_ID4.fastq -n ID1 ID2 ID3 ID4 --outdir compare-IDs

