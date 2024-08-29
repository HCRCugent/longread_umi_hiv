#!/bin/bash

source activate longread_umi_HIV

#These should be fixed as standard
#set working directory
wrk=/data/thesisguest1/work/UMI

#set the following parameters
#what is the name of your project folder where everything will run?


run_ID=RUN50

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

#
for bc in {01..06}; do 
  mkdir $prj/$sample_id"$bc"
done
##
#
### Move data from guppy pass output to single in work directory fastq files 
for bc in {01..06}; do 
  cat $raw_data_loc/barcode"$bc"/*.gz > $prj/$sample_id"$bc".fastq.gz
  gunzip $prj/$sample_id"$bc".fastq.gz
done
### Because fastq files were not zipped: 
#for bc in {01..03}; do 
  #cat $raw_data_loc/barcode"$bc"/*.fastq > $prj/$sample_id"$bc".fastq
#done
#
#
########### FILTER FILES for HIV #############
##https://www.samformat.info/sam-format-flag
for bc in {01..06}; do 
  minimap2 -t 30 -ax map-ont /data/thesisguest1/Pipe_3SEP20/db/blast/HXB2.fasta $prj/$sample_id"$bc".fastq | samtools fastq -n -f 0x4 - > $prj/$sample_id"$bc"/$sample_id"$bc".nonHIV.fastq 
  minimap2 -t 30 -ax map-ont /data/thesisguest1/Pipe_3SEP20/db/blast/HXB2.fasta $prj/$sample_id"$bc".fastq | samtools fastq -n -F 0x4 - > $prj/$sample_id"$bc"/$sample_id"$bc".HIV.fastq 
done


# updated medaka model () & racon 3x + medaka 2x

for bc in {01..06}; do 
longread_umi nanopore_pipeline -d $prj/$sample_id"$bc"/$sample_id"$bc".HIV.fastq -o $prj/$sample_id"$bc"/umi_out_all_length_$sample_id"$bc"/ -v 30 -s 200 -e 200 -m 100 -M 10000 -f $ont_f -F $hiv_f  -r $ont_r -R $hiv_r -c 3 -p 2 -q r1041_e82_400bps_sup_v4.2.0 -U 'r103_min_high_g360' -t 30 -T 8 > $prj/$sample_id"$bc"/pipeline_sample_all_sup_nohup.log
done
