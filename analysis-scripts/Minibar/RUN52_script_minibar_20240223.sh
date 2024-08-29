#!/bin/bash

source activate longread_umi_HIV

#These should be fixed as standard
#set working directory
wrk=/data/thesisguest1/work/UMI

#set the following parameters
#what is the name of your project folder where everything will run?


run_ID=RUN52

prj=$wrk/data_"$run_ID"

sample_id="$run_ID"_guppy_v6.5.7_sup_BC

raw_data_loc=$prj/raw_fastq_guppy_v6.5.7_sup

cd $prj

#primers:
#outer tagging:
pin_inner_f=AAGTAGTGTGTGCCCGTCTGTTGTGTGAC
pin_inner_r=GGAAAGTCCCCAGCGGAAAGTCCCTTGTAG

hiv_f=$pin_inner_f
hiv_r=$pin_inner_r

ont_f=CAAGCAGAAGACGGCATACGAGAT
ont_r=AATGATACGGCGACCACCGAGATC


source activate minibar || conda activate minibar
minibar=/data/thesisguest1/Tools/minibar.py

#
source activate minibar

#for bc in {07..12}; do 
#  mkdir $prj/$sample_id"$bc"
#  mkdir $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/
#  cd $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/
#  $minibar $prj/Replicate_ID.txt $prj/$sample_id"$bc".fastq -F -l 600
#  find ./ -type f -name "*_*_*" -exec rm -f '{}' ';'
#done

source activate longread_umi_HIV || conda activate longread_umi_HIV
for bc in {07..12}; do 
for rep_id in {1..8};do
minimap2 -t 15 -ax map-ont /data/thesisguest1/Pipe_3SEP20/db/blast/HXB2.fasta $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/sample_ID"$rep_id".fastq | samtools fastq -n -f 0x4 - > $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/BC"$bc"_sample_ID"$rep_id".nonHIV.fastq ;
minimap2 -t 15 -ax map-ont /data/thesisguest1/Pipe_3SEP20/db/blast/HXB2.fasta $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/sample_ID"$rep_id".fastq  | samtools fastq -n -F 0x4 - > $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/BC"$bc"_sample_ID"$rep_id".HIV.fastq ;
done
done


for rep_id in {1..8}; do echo "$rep_id";
NanoComp --fastq $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/sample_ID"$rep_id".HIV.fastq $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/sample_ID"$rep_id".nonHIV.fastq -n ID"$rep_id"_HIV ID"$rep_id"_nonHIV --outdir compare-ID"$rep_id"_HIV
done
source activate ONT || conda activate ONT

for bc in {07..12}; do 
for rep_id in {1..6};do
cd $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"
NanoComp --fastq $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/BC"$bc"_sample_ID"$rep_id".HIV.fastq $prj/$sample_id"$bc"/minibar_demultiplex_$sample_id"$bc"/BC"$bc"_sample_ID"$rep_id".nonHIV.fastq -n BC"$bc"-ID"$rep_id"_HIV BC"$bc"-ID"$rep_id"_nonHIV --outdir compare-BC"$bc"-ID"$rep_id"_HIV -t 20
done
done